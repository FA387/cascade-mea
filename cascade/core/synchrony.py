"""
core/synchrony.py
-----------------
Spike-train synchrony metrics, matching the set computed by the Axion
Neural Metric Tool but applied uniformly to every manufacturer CASCADE supports.

Three families:

1. Cross-correlogram metrics  (Axion Neural Metric Tool v4.1, §2.3 / §4.4)
   Area under the well-wide pooled inter-electrode cross-correlation, the same
   area for the auto-correlation-normalised correlogram, and the full width at
   half height of each. Computed in the frequency domain, pooled over all unique
   electrode pairs.

2. Kreuz SPIKE-distance      (Kreuz et al. 2013, J Neurophysiol 109:1457-72)
   Parameter-free, time-scale independent. The "improved" dissimilarity profile
   of Eq 19, not the original Eq 13, which gave spuriously high values for
   event-like firing. Axion reports 1 - D_S so that 1 is perfect synchrony.

3. Synchrony Index           (Paiva, Park & Principe 2010, Neural Comput Appl
                              19:405-419, §2.3 — the Schreiber et al. measure)
   Gaussian-kernel correlation in [0, 1], 1 = perfect synchrony. Evaluated with
   the closed form of Eq 11, so no discretisation of the spike trains is needed.

Cross-manufacturer note
-----------------------
All timing parameters here are ABSOLUTE durations in seconds, never derived from
the sampling rate or the mean ISI. The avalanche pipeline deliberately infers its
bin width from the data (`core.avalanche.infer_dt`); synchrony deliberately does
not. If the correlogram bin scaled with sampling rate, "area under the normalised
cross-correlation at a 20 ms window" would mean a different physical quantity on a
12.5 kHz Axion plate than on a 19.75 kHz 3Brain plate, and the numbers would not
be comparable across manufacturers. That comparability is the point of CASCADE, so
it is preserved here by construction.
"""

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

# PySpike (Mulansky & Kreuz) is the reference implementation of the SPIKE-distance,
# written by an author of the 2013 paper. When it is installed CASCADE defers to it
# rather than to the local implementation below.
#
# Why: a direct comparison of the two showed CASCADE's own port agreeing to ~1e-3 on
# Poisson trains but differing by up to 2.7e-2 on structured ones. The published
# paper fixes the dissimilarity profile (Eqs 16-19, which both implementations agree
# on — PySpike's `dist_at_t` is algebraically identical to Eq 19) but leaves the
# treatment of the first and last interval underspecified beyond "add auxiliary
# spikes at t=0 and t=T", and PySpike's interval bookkeeping differs there. Rather
# than ship a near-miss of a published measure, CASCADE uses the reference when it
# can and records which path produced the number.
try:                                            # pragma: no cover - optional dep
    import pyspike as _pyspike
    HAVE_PYSPIKE = True
except Exception:                               # pragma: no cover
    _pyspike = None
    HAVE_PYSPIKE = False


# ── binning helper ─────────────────────────────────────────────────────────────
def _binary_matrix(spike_times: Dict[int, np.ndarray],
                   duration_s: float,
                   bin_s: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Bin every channel's spike times into a counts matrix.

    Parameters
    ----------
    spike_times : {channel_id: array of spike times in seconds}
    duration_s  : recording duration in seconds
    bin_s       : bin width in seconds (absolute — see module docstring)

    Returns
    -------
    X      : (n_channels, n_bins) float array of spike counts
    ch_ids : 1-D array of channel ids, row order of X
    """
    ch_ids = np.array(sorted(spike_times.keys()))
    n_bins = max(1, int(np.ceil(duration_s / bin_s)))
    edges = np.arange(n_bins + 1, dtype=np.float64) * bin_s

    X = np.zeros((ch_ids.size, n_bins), dtype=np.float64)
    for r, ch in enumerate(ch_ids):
        t = np.asarray(spike_times[ch], dtype=np.float64)
        t = t[(t >= 0.0) & (t < n_bins * bin_s)]
        if t.size:
            X[r] = np.histogram(t, bins=edges)[0]
    return X, ch_ids


# ── pooled cross-correlogram ───────────────────────────────────────────────────
def pooled_cross_correlogram(spike_times: Dict[int, np.ndarray],
                             duration_s: float,
                             bin_s: float = 0.001,
                             max_lag_s: float = 0.5,
                             ) -> Optional[dict]:
    """
    Well-wide pooled inter-electrode cross-correlogram, plus the
    auto-correlation-normalised version.

    Pooled over all unique pairs i < j. Computed in the frequency domain using

        sum_{i<j} conj(F_i) * F_j = ( |sum_i F_i|^2 - sum_i |F_i|^2 ) / 2

    which costs O(N * T log T) rather than O(N^2 * T log T) — for a 341-channel
    recording that is ~58,000 pairwise FFT products avoided. This is the same
    frequency-domain strategy the Axion guide cites (Halliday, Rosenberg, Breeze &
    Conway 2006).

    The normalised correlogram divides by sqrt(A_ii(0) * A_jj(0)) summed over the
    same pairs, i.e. it scales the cross-correlations by the geometric mean of the
    zero-lag auto-correlations. Per the Axion guide this "reduces the effects of
    high mean firing rate or high spiking regularity on individual electrodes or
    units which can artificially increase correlations".

    Returns
    -------
    dict with
      lags_s     : (2*L+1,) lag axis in seconds, symmetric around zero
      cc         : pooled cross-correlogram (counts of coincidences per lag bin)
      cc_norm    : auto-correlation-normalised correlogram
      n_pairs    : number of unique channel pairs pooled
      bin_s      : the bin width used
    None if fewer than 2 channels have spikes.
    """
    if not spike_times:
        return None

    X, ch_ids = _binary_matrix(spike_times, duration_s, bin_s)
    # drop silent channels — they contribute nothing and would divide by zero
    keep = X.sum(axis=1) > 0
    X, ch_ids = X[keep], ch_ids[keep]
    n_ch, n_bins = X.shape
    if n_ch < 2:
        return None

    max_lag_bins = int(round(max_lag_s / bin_s))
    max_lag_bins = max(1, min(max_lag_bins, n_bins - 1))

    # zero-pad to avoid circular wrap-around in the FFT correlation
    n_fft = int(2 ** np.ceil(np.log2(2 * n_bins)))

    # Accumulate sum_i F_i and sum_i |F_i|^2 one channel at a time rather than
    # materialising the whole (n_ch, n_fft/2+1) complex spectrum. Holding every
    # channel's FFT would cost O(N*T): for a 341-channel, 300 s recording at 1 ms
    # bins that is ~1.4 GB, which would dominate the pipeline's whole memory
    # footprint. Streaming makes it O(T) regardless of channel count.
    n_freq = n_fft // 2 + 1
    sum_F = np.zeros(n_freq, dtype=np.complex128)
    auto_spec = np.zeros(n_freq, dtype=np.float64)
    for r in range(n_ch):
        Fi = np.fft.rfft(X[r], n=n_fft)
        sum_F += Fi
        auto_spec += np.abs(Fi) ** 2

    # pooled cross-spectrum over unique pairs (see identity above)
    cross_spec = (np.abs(sum_F) ** 2 - auto_spec) / 2.0

    cc_full = np.fft.irfft(cross_spec, n=n_fft)
    ac_full = np.fft.irfft(auto_spec, n=n_fft)

    # assemble a symmetric lag axis; the pooled cross-correlogram is symmetric by
    # construction because summing over unordered pairs includes both directions
    lags = np.arange(-max_lag_bins, max_lag_bins + 1)
    cc = np.concatenate([cc_full[-max_lag_bins:], cc_full[:max_lag_bins + 1]])

    n_pairs = n_ch * (n_ch - 1) // 2

    # normalisation: geometric mean of zero-lag auto-correlations, summed over the
    # same unique pairs. sum_{i<j} sqrt(A_ii(0) A_jj(0)) = ( (sum_i s_i)^2 - sum_i s_i^2 ) / 2
    # with s_i = sqrt(A_ii(0)).
    zero_lag_auto = (X ** 2).sum(axis=1)          # A_ii(0) per channel
    s = np.sqrt(zero_lag_auto)
    denom = (s.sum() ** 2 - (s ** 2).sum()) / 2.0
    cc_norm = cc / denom if denom > 0 else np.zeros_like(cc)

    return dict(lags_s=lags * bin_s,
                cc=cc,
                cc_norm=cc_norm,
                ac=ac_full[:max_lag_bins + 1],
                n_pairs=int(n_pairs),
                n_channels=int(n_ch),
                bin_s=float(bin_s))


def area_under_correlogram(lags_s: np.ndarray,
                           cc: np.ndarray,
                           window_s: float) -> float:
    """
    Area under the correlogram within +/- window_s of zero lag.

    This is Axion's "Area Under Cross-Correlation" / "Area Under Normalized
    Cross-Correlation", the synchrony window being the "window of time around zero
    that is used to compute the area under the cross-correlation".

    Integrated with the trapezoid rule over the lag axis, so the value has units of
    correlogram-amplitude x seconds and does not depend on the bin width chosen.
    """
    if lags_s.size == 0:
        return float("nan")
    m = np.abs(lags_s) <= window_s
    if m.sum() < 2:
        return float("nan")
    return float(np.trapezoid(cc[m], lags_s[m]))


def fwhh_is_resolvable(lags_s: np.ndarray, cc: np.ndarray, bin_s: float,
                       fwhh: float) -> Tuple[bool, str]:
    """
    Decide whether a measured FWHH is a real width or an artefact of the settings.

    Two ways it can fail, both of which still produce a finite-looking number:

    * **bin-limited** — the peak is narrower than the correlogram bin, so the
      measured width just tracks `--sync_bin`. Halve the bin until FWHH stops
      changing.
    * **lag-truncated** — the correlogram never falls to half height inside the
      computed lag range, so the peak is wider than `--sync_max_lag` and has been
      cut off. Widen the lag.

    Returns (resolvable, reason). `reason` is "" when the value is trustworthy.
    The FWHH itself is never altered — this only reports on it.
    """
    if not np.isfinite(fwhh):
        return False, "no half-height crossing found within the lag range"
    if fwhh <= 2.5 * bin_s:
        return False, (f"peak is narrower than ~2 bins (FWHH {fwhh*1000:.2f} ms at a "
                       f"{bin_s*1000:.2f} ms bin) — value is bin-limited, "
                       f"reduce --sync_bin")
    span = float(np.max(lags_s)) if lags_s.size else 0.0
    if fwhh >= 0.9 * span:
        return False, (f"peak fills the lag range (FWHH {fwhh*1000:.1f} ms vs "
                       f"±{span*1000:.0f} ms) — value is lag-truncated, "
                       f"increase --sync_max_lag")
    return True, ""


def full_width_half_height(lags_s: np.ndarray, cc: np.ndarray) -> float:
    """
    Full width at half height of the correlogram, in seconds.

    Axion: "Distance along the x-axis (phase lag) from left half height to right
    half height (probability) of the cross-correlogram. Higher full widths indicate
    a wider correlogram (less synchrony) whereas lower full widths indicate a taller
    correlogram (greater synchrony)."

    Height is measured from the correlogram's baseline (its minimum over the
    returned lag range) up to its peak, and the crossings are linearly interpolated
    between bins so the result is not quantised to the bin width. NaN if the
    correlogram never falls to half height within the computed lag range, which is
    itself informative — it means the lag window was too narrow to contain the peak.
    """
    if cc.size < 3:
        return float("nan")
    base = float(np.min(cc))
    peak_i = int(np.argmax(cc))
    peak = float(cc[peak_i])
    if not np.isfinite(peak) or peak <= base:
        return float("nan")
    half = base + 0.5 * (peak - base)

    # walk left from the peak to the first bin at or below half height
    li = peak_i
    while li > 0 and cc[li] > half:
        li -= 1
    if cc[li] > half:
        return float("nan")
    # walk right
    ri = peak_i
    n = cc.size
    while ri < n - 1 and cc[ri] > half:
        ri += 1
    if cc[ri] > half:
        return float("nan")

    def _interp(i_out, i_in):
        """Linear crossing between an outside bin and an inside bin."""
        y0, y1 = cc[i_out], cc[i_in]
        if y1 == y0:
            return lags_s[i_in]
        f = (half - y0) / (y1 - y0)
        return lags_s[i_out] + f * (lags_s[i_in] - lags_s[i_out])

    left = _interp(li, min(li + 1, peak_i)) if li != peak_i else lags_s[peak_i]
    right = _interp(ri, max(ri - 1, peak_i)) if ri != peak_i else lags_s[peak_i]
    return float(abs(right - left))


# ── Kreuz SPIKE-distance ───────────────────────────────────────────────────────
def _corner_quantities(t1: np.ndarray, t2: np.ndarray, T: float):
    """
    Build the piecewise-linear SPIKE-distance profile S(t) for one pair of trains.

    Implements Kreuz et al. (2013) Eqs 3-5, 9-10, 16-17 and 19 — the *improved*
    profile. Auxiliary spikes are added at t=0 and t=T as the paper requires, which
    is what resolves the first/last ISI ambiguity.

    Returns (grid, S) where S is evaluated on `grid`, the sorted union of all spike
    times from both trains plus the endpoints. S(t) is piecewise linear between
    consecutive grid points, so trapezoidal integration over this grid is exact and
    no time sampling is needed.
    """
    # auxiliary leading/trailing spikes (paper, §2.1)
    a1 = np.unique(np.concatenate(([0.0], t1, [T])))
    a2 = np.unique(np.concatenate(([0.0], t2, [T])))

    grid = np.unique(np.concatenate((a1, a2)))

    # S(t) is piecewise linear but DISCONTINUOUS at spike times: crossing a spike
    # changes which corner spikes are in play, so the profile generally jumps. The
    # integral must therefore be taken interval by interval, using the left and
    # right limits within each open interval, rather than one value per grid point
    # (which would reuse the right-limit at both ends and bias the result upward).
    # Within an open interval the four corner spikes are constant and S is linear,
    # so evaluating at the two endpoints with those fixed corners is exact.
    n_iv = grid.size - 1
    S_left = np.zeros(n_iv, dtype=np.float64)    # S at the start of each interval
    S_right = np.zeros(n_iv, dtype=np.float64)   # S at the end of each interval

    for k in range(n_iv):
        t_a, t_b = grid[k], grid[k + 1]
        if t_b <= t_a:
            continue
        mid = 0.5 * (t_a + t_b)

        # corner spikes for the OPEN interval, Eqs 3-4 evaluated at the midpoint
        i1 = np.searchsorted(a1, mid, side="right") - 1
        i2 = np.searchsorted(a2, mid, side="right") - 1
        i1 = min(max(i1, 0), a1.size - 1)
        i2 = min(max(i2, 0), a2.size - 1)
        tP1, tP2 = a1[i1], a2[i2]
        tF1 = a1[min(i1 + 1, a1.size - 1)]
        tF2 = a2[min(i2 + 1, a2.size - 1)]

        xISI1 = tF1 - tP1          # Eq 5
        xISI2 = tF2 - tP2
        if xISI1 <= 0 or xISI2 <= 0:
            continue

        # Eq 16: each corner spike against its NEAREST counterpart in the other
        # train. This is the 2013 improvement; comparing P-to-P and F-to-F instead
        # is the original definition and reintroduces the spurious high values.
        dP1 = np.min(np.abs(a2 - tP1))
        dF1 = np.min(np.abs(a2 - tF1))
        dP2 = np.min(np.abs(a1 - tP2))
        dF2 = np.min(np.abs(a1 - tF2))

        mean_isi = 0.5 * (xISI1 + xISI2)
        if mean_isi <= 0:
            continue
        norm = 2.0 * mean_isi ** 2

        for which, t in ((0, t_a), (1, t_b)):
            xP1, xF1 = t - tP1, tF1 - t     # Eqs 9-10
            xP2, xF2 = t - tP2, tF2 - t
            # Eq 17, per train
            S1 = (dP1 * xF1 + dF1 * xP1) / xISI1
            S2 = (dP2 * xF2 + dF2 * xP2) / xISI2
            # Eq 19: cross-train ISI weighting (the firing-rate correction). Using
            # Eq 18 here instead would drop that correction.
            val = (S1 * xISI2 + S2 * xISI1) / norm
            if which == 0:
                S_left[k] = val
            else:
                S_right[k] = val

    return grid, S_left, S_right


def spike_distance_backend() -> str:
    """Which SPIKE-distance implementation is active: "pyspike" or "cascade"."""
    return "pyspike" if HAVE_PYSPIKE else "cascade"


def spike_distance(t1: np.ndarray, t2: np.ndarray, T: float) -> float:
    """
    Bivariate Kreuz SPIKE-distance D_S, Eq 1: the time average of S(t) over [0, T].

    Parameter-free and time-scale invariant. Bounded in [0, 1]; zero only for
    identical spike trains.

    Uses PySpike when installed (the reference implementation, see module header);
    falls back to the local port otherwise. `spike_distance_backend()` reports which
    path is active, and it is recorded in the synchrony output.
    """
    if T <= 0:
        return float("nan")
    if HAVE_PYSPIKE:
        a = np.asarray(t1, dtype=np.float64)
        b = np.asarray(t2, dtype=np.float64)
        if a.size == 0 or b.size == 0:
            return float("nan")
        return float(_pyspike.spike_distance(
            _pyspike.SpikeTrain(a, [0.0, T]),
            _pyspike.SpikeTrain(b, [0.0, T])))
    grid, S_left, S_right = _corner_quantities(np.asarray(t1, float),
                                               np.asarray(t2, float), T)
    if grid.size < 2:
        return float("nan")
    # trapezoid per interval, honouring the discontinuities at spike times
    widths = np.diff(grid)
    return float(np.sum(0.5 * (S_left + S_right) * widths) / T)


def averaged_spike_distance(spike_times: Dict[int, np.ndarray],
                            duration_s: float,
                            max_channels: Optional[int] = None,
                            rng_seed: int = 0) -> dict:
    """
    Averaged bivariate SPIKE-distance over all unique channel pairs (Eq 2), and the
    Kreuz Metric that Axion reports.

    Cost is O(N^2) pairs, each O(M log M) in the spike counts, so on dense arrays
    this dominates the whole pipeline. `max_channels` subsamples channels (seeded,
    so it is reproducible) and the result records what was used.

    Returns
    -------
    dict: spike_distance, kreuz_metric (= 1 - D_S), n_pairs, n_channels, subsampled
    """
    ch_ids = np.array(sorted(k for k, v in spike_times.items()
                             if np.asarray(v).size > 0))
    subsampled = False
    if max_channels is not None and ch_ids.size > max_channels:
        rng = np.random.default_rng(rng_seed)
        ch_ids = np.sort(rng.choice(ch_ids, size=max_channels, replace=False))
        subsampled = True

    if ch_ids.size < 2:
        return dict(spike_distance=float("nan"), kreuz_metric=float("nan"),
                    n_pairs=0, n_channels=int(ch_ids.size), subsampled=subsampled)

    vals = []
    for a in range(ch_ids.size - 1):
        ta = np.asarray(spike_times[ch_ids[a]], float)
        for b in range(a + 1, ch_ids.size):
            vals.append(spike_distance(ta,
                                       np.asarray(spike_times[ch_ids[b]], float),
                                       duration_s))
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return dict(spike_distance=float("nan"), kreuz_metric=float("nan"),
                    n_pairs=0, n_channels=int(ch_ids.size), subsampled=subsampled)

    ds = float(vals.mean())
    return dict(spike_distance=ds,
                kreuz_metric=1.0 - ds,
                n_pairs=int(vals.size),
                n_channels=int(ch_ids.size),
                subsampled=subsampled,
                backend=spike_distance_backend())


# ── Synchrony Index (Paiva / Schreiber) ────────────────────────────────────────
def synchrony_index_pair(t1: np.ndarray, t2: np.ndarray,
                         sigma_s: float,
                         n_sigma: float = 5.0) -> float:
    """
    Gaussian-kernel correlation r(S_i, S_j) in [0, 1] — Paiva et al. (2010) Eq 9,
    evaluated with the closed form of Eq 11:

        r = sum_mn exp[-(t_m^i - t_n^j)^2 / 2 sigma^2]
            / sqrt( sum_mn exp[-(t_m^i-t_n^i)^2/2s^2] * sum_mn exp[-(t_m^j-t_n^j)^2/2s^2] )

    1 = perfect synchrony. This is Axion's "Synchrony Index".

    The Gaussian decays fast, so only pairs within `n_sigma` sigma contribute
    measurably; the double sums are windowed with searchsorted, which turns the
    nominal O(N_i N_j) into roughly O(N_i k). At n_sigma = 5 the neglected terms
    are below exp(-12.5) ~ 4e-6 of the peak.

    Note (Paiva §2.3): the derived dissimilarity 1 - r is a semi-metric — symmetric
    and positive, but it does NOT satisfy the triangle inequality.
    """
    t1 = np.sort(np.asarray(t1, dtype=np.float64))
    t2 = np.sort(np.asarray(t2, dtype=np.float64))
    if t1.size == 0 or t2.size == 0 or sigma_s <= 0:
        return float("nan")

    two_s2 = 2.0 * sigma_s * sigma_s
    cut = n_sigma * sigma_s

    def _windowed_sum(a: np.ndarray, b: np.ndarray) -> float:
        """sum over pairs exp(-(a_m-b_n)^2 / 2 sigma^2), windowed to +/- cut."""
        lo = np.searchsorted(b, a - cut, side="left")
        hi = np.searchsorted(b, a + cut, side="right")
        total = 0.0
        for m in range(a.size):
            if hi[m] > lo[m]:
                d = a[m] - b[lo[m]:hi[m]]
                total += float(np.exp(-(d * d) / two_s2).sum())
        return total

    cross = _windowed_sum(t1, t2)
    a11 = _windowed_sum(t1, t1)
    a22 = _windowed_sum(t2, t2)
    denom = np.sqrt(a11 * a22)
    if denom <= 0:
        return float("nan")
    return float(cross / denom)


def pooled_synchrony_index(spike_times: Dict[int, np.ndarray],
                           sigma_s: float,
                           max_channels: Optional[int] = None,
                           rng_seed: int = 0) -> dict:
    """
    Synchrony Index averaged over all unique channel pairs.

    `sigma_s` is the timescale. Paiva stresses that this parameter sets what the
    measure reports: small relative to the mean ISI and it quantifies spike-timing
    coincidence; large and it approaches a firing-rate comparison. CASCADE defaults
    it to the synchrony window so the correlogram metrics and this one describe the
    same timescale.
    """
    ch_ids = np.array(sorted(k for k, v in spike_times.items()
                             if np.asarray(v).size > 0))
    subsampled = False
    if max_channels is not None and ch_ids.size > max_channels:
        rng = np.random.default_rng(rng_seed)
        ch_ids = np.sort(rng.choice(ch_ids, size=max_channels, replace=False))
        subsampled = True

    if ch_ids.size < 2:
        return dict(synchrony_index=float("nan"), n_pairs=0,
                    n_channels=int(ch_ids.size), subsampled=subsampled)

    vals = []
    for a in range(ch_ids.size - 1):
        ta = spike_times[ch_ids[a]]
        for b in range(a + 1, ch_ids.size):
            vals.append(synchrony_index_pair(ta, spike_times[ch_ids[b]], sigma_s))
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return dict(synchrony_index=float("nan"), n_pairs=0,
                    n_channels=int(ch_ids.size), subsampled=subsampled)

    return dict(synchrony_index=float(vals.mean()),
                n_pairs=int(vals.size),
                n_channels=int(ch_ids.size),
                subsampled=subsampled)


# ── assembly ───────────────────────────────────────────────────────────────────
def compute_synchrony(spike_times: Dict[int, np.ndarray],
                      duration_s: float,
                      window_s: float = 0.020,
                      bin_s: float = 0.001,
                      max_lag_s: float = 0.5,
                      sigma_s: Optional[float] = None,
                      metrics: str = "xcorr",
                      max_channels: Optional[int] = 64,
                      ) -> Tuple[dict, Optional[dict]]:
    """
    Compute the requested synchrony metrics.

    Parameters
    ----------
    metrics : "none" | "xcorr" | "kreuz" | "index" | "all"
              Mirrors Axion's "Additional Synchrony Metrics" dropdown. The
              correlogram family is cheap; kreuz and index are O(N^2) in channels.
    window_s : synchrony window (seconds). 0.005 emphasises millisecond-timescale
              synchrony, 0.100 captures slow timescales.
    sigma_s : timescale for the Synchrony Index; defaults to window_s.
    max_channels : cap on channels used for the O(N^2) measures. None = no cap.

    Returns
    -------
    (metrics_dict, correlogram_dict_or_None)
    """
    metrics = (metrics or "none").lower()
    if metrics not in {"none", "core", "xcorr", "kreuz", "index", "all"}:
        raise ValueError(f"unknown metrics selection: {metrics!r}")

    out: dict = dict(sync_window_s=float(window_s))
    if metrics == "none":
        return out, None

    # "core" is the default: the two metrics that both resolve on every
    # manufacturer and actually discriminate between recordings. Measured across
    # the five published sample recordings:
    #
    #   metric             spread(max/min)   CV    resolves on
    #   area_xcorr_norm          26x        0.78       5/5   <- core
    #   synchrony_index          11x        0.69       5/5   <- core
    #   area_xcorr (raw)         18x        0.75       5/5   (rate-confounded)
    #   fwhh_xcorr_norm         552x        1.47       2/5   (mostly artefact)
    #   kreuz_metric            1.2x        0.09       5/5   (near-constant)
    #
    # The others remain available explicitly. Kreuz is kept because it is
    # parameter-free, which makes it well suited to tracking change within one
    # culture over time even though it separates different cultures poorly.
    want_xcorr = metrics in {"core", "xcorr", "all"}
    want_kreuz = metrics in {"kreuz", "all"}
    want_index = metrics in {"core", "index", "all"}
    core_only = metrics == "core"

    cc = None
    if want_xcorr:
        cc = pooled_cross_correlogram(spike_times, duration_s,
                                      bin_s=bin_s, max_lag_s=max_lag_s)
        if cc is None:
            out.update(area_xcorr=float("nan"), area_xcorr_norm=float("nan"),
                       fwhh_xcorr=float("nan"), fwhh_xcorr_norm=float("nan"),
                       n_pairs_xcorr=0)
        else:
            f_raw = full_width_half_height(cc["lags_s"], cc["cc"])
            f_norm = full_width_half_height(cc["lags_s"], cc["cc_norm"])
            ok, why = fwhh_is_resolvable(cc["lags_s"], cc["cc_norm"],
                                         cc["bin_s"], f_norm)
            out.update(
                area_xcorr=area_under_correlogram(cc["lags_s"], cc["cc"], window_s),
                area_xcorr_norm=area_under_correlogram(cc["lags_s"], cc["cc_norm"],
                                                       window_s),
                fwhh_xcorr=f_raw,
                fwhh_xcorr_norm=f_norm,
                # The FWHH values above are reported exactly as computed. This flag
                # says whether the settings could actually resolve them, so a
                # resolution floor is not mistaken for a measurement.
                fwhh_resolution_limited=(not ok),
                fwhh_limit_reason=why,
                n_pairs_xcorr=cc["n_pairs"],
                xcorr_bin_s=cc["bin_s"],
                xcorr_max_lag_s=float(np.max(cc["lags_s"])) if cc["lags_s"].size else 0.0,
            )
            if core_only:
                # keep the headline table to what is trustworthy everywhere
                for k in ("area_xcorr", "fwhh_xcorr", "fwhh_xcorr_norm",
                          "fwhh_resolution_limited", "fwhh_limit_reason"):
                    out.pop(k, None)

    if want_kreuz:
        k = averaged_spike_distance(spike_times, duration_s,
                                    max_channels=max_channels)
        out.update(spike_distance=k["spike_distance"],
                   kreuz_metric=k["kreuz_metric"],
                   n_pairs_kreuz=k["n_pairs"],
                   kreuz_subsampled=k["subsampled"],
                   kreuz_backend=k.get("backend", "cascade"))

    if want_index:
        sg = float(sigma_s) if sigma_s else float(window_s)
        s = pooled_synchrony_index(spike_times, sg, max_channels=max_channels)
        out.update(synchrony_index=s["synchrony_index"],
                   sync_index_sigma_s=sg,
                   n_pairs_index=s["n_pairs"],
                   index_subsampled=s["subsampled"])

    return out, cc


def synchrony_to_dataframe(m: dict) -> pd.DataFrame:
    """Flatten the metrics dict into a two-column metric/value table."""
    rows = [(k, m[k]) for k in m]
    return pd.DataFrame(rows, columns=["metric", "value"])
