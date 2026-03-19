"""
core/avalanche.py
-----------------
Neuronal avalanche detection.

Method follows Habibollahi et al. (2022) / Beggs & Plenz (2003):
  1. Bin population spike counts with bin width Δt.
  2. Apply a percentile threshold to the population activity.
  3. An avalanche = run of consecutive bins ABOVE the threshold.
  4. Size S = total spikes in that run.
  5. Duration D = number of bins in that run.

Key parameters
--------------
dt   : bin width in seconds.
       Default: inferred as the mean inter-spike interval (ISI) of the
       population, clipped to [1 ms, 50 ms].  This is the standard
       data-driven choice (Beggs & Plenz 2003).

perc : percentile threshold applied to the distribution of bin counts.
       0.0  → threshold = 0, i.e. any non-empty bin contributes
              (appropriate for clean in-vitro MEA spike data).
       0.25 → threshold = 25th percentile (Habibollahi 2022 default,
              useful when there is a non-zero background firing rate).
"""

import numpy as np
import pandas as pd


# ── Bin-width inference ────────────────────────────────────────────────────────
def infer_dt(all_spike_times: np.ndarray,
             min_dt: float = 0.001,
             max_dt: float = 0.050) -> float:
    """
    Estimate avalanche bin width from the mean population ISI.

    Parameters
    ----------
    all_spike_times : sorted 1-D array of spike times in seconds
    min_dt / max_dt : hard clipping limits (seconds)

    Returns
    -------
    dt in seconds.
    """
    t = np.sort(all_spike_times)
    if t.size < 2:
        return 0.01
    return float(np.clip(np.mean(np.diff(t)), min_dt, max_dt))


# ── Population counts per bin ──────────────────────────────────────────────────
def counts_per_dt(all_spike_times: np.ndarray, dt: float):
    """
    Bin population spike times into fixed-width bins.

    Returns
    -------
    counts : 1-D int array  – number of spikes in each bin
    edges  : 1-D float array – bin edges (length = len(counts)+1)
    """
    if all_spike_times.size == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.float64)

    t0 = float(all_spike_times.min())
    t1 = float(all_spike_times.max())
    nb = max(1, int(np.ceil((t1 - t0) / dt)))
    edges  = np.linspace(t0, t0 + nb * dt, nb + 1)
    counts, _ = np.histogram(all_spike_times, bins=edges)
    return counts.astype(np.int64), edges


# ── Avalanche segment finder ───────────────────────────────────────────────────
def avalanche_segments(counts: np.ndarray, threshold: float = 0.0) -> list:
    """
    Find runs of consecutive bins where counts > threshold.

    Returns
    -------
    List of (start_idx, end_idx) pairs (end_idx is exclusive).
    """
    above = counts > threshold
    segs  = []
    i     = 0
    while i < above.size:
        if above[i]:
            j = i
            while j < above.size and above[j]:
                j += 1
            segs.append((i, j))
            i = j
        else:
            i += 1
    return segs


# ── Main avalanche computation ─────────────────────────────────────────────────
def compute_avalanches(all_spike_times: np.ndarray,
                       dt: float,
                       perc: float = 0.0) -> pd.DataFrame:
    """
    Detect neuronal avalanches from the population spike train.

    Parameters
    ----------
    all_spike_times : 1-D array of spike times in seconds (all channels)
    dt              : bin width in seconds
    perc            : percentile threshold (0.0 = no threshold)
                      Use 0.0 for clean sorted-spike MEA data.
                      Use 0.25 to match Habibollahi 2022 default.

    Returns
    -------
    DataFrame with columns:
        avalanche_size  (int)  – total spikes in the avalanche
        duration_bins   (int)  – number of Δt bins spanned
    """
    if all_spike_times.size == 0:
        return pd.DataFrame(columns=["avalanche_size", "duration_bins"])

    counts, _ = counts_per_dt(all_spike_times, dt)

    # percentile threshold on population activity
    if perc > 0.0:
        sorted_c  = np.sort(counts)
        idx       = int(round(len(sorted_c) * float(perc)))
        idx       = max(0, min(idx, len(sorted_c) - 1))
        threshold = float(sorted_c[idx])
    else:
        threshold = 0.0

    segs   = avalanche_segments(counts, threshold)
    sizes  = []
    durs   = []

    for start, end in segs:
        # size = total spikes above the threshold baseline
        raw = float(np.sum(counts[start:end]))
        dur = int(end - start)
        # subtract threshold contribution following Habibollahi / Poil 2012:
        # size = raw - threshold * (duration - 1)
        # Reference: AV_analysis_BurstT.py line 71
        #   fired = sum(activity) - threshold2*(J[i+1]-J[i]-2)
        #   where (J[i+1]-J[i]-2) = duration - 1
        # for perc=0 threshold=0 so behaviour is unchanged
        size = max(1, int(round(raw - threshold * (dur - 1))))
        sizes.append(size)
        durs.append(dur)

    return pd.DataFrame({"avalanche_size": sizes, "duration_bins": durs})
