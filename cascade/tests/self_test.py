"""
CASCADE built-in self-test
==========================
Generates a synthetic MEA recording and runs every core module.
No real data files required.

Usage:
    cascade --test
"""

import sys
import traceback
import numpy as np

from cascade.utils.common import MEARecording


# ── terminal colours ──────────────────────────────────────────────────────────
try:
    _tty = sys.stdout.isatty()
except AttributeError:
    _tty = False

_GRN  = "\033[32m" if _tty else ""
_RED  = "\033[31m" if _tty else ""
_YEL  = "\033[33m" if _tty else ""
_CYN  = "\033[36m" if _tty else ""
_RST  = "\033[0m"  if _tty else ""
_BOLD = "\033[1m"  if _tty else ""

_W = 55   # label column width


def _ok(label: str):
    pad = "." * max(2, _W - len(label))
    print(f"  {label} {pad} {_GRN}PASS{_RST}")


def _fail(label: str, err: Exception):
    pad = "." * max(2, _W - len(label))
    print(f"  {label} {pad} {_RED}FAIL{_RST}")
    for line in traceback.format_exception_only(type(err), err):
        print(f"    {_RED}{line.rstrip()}{_RST}")


def _skip(label: str, reason: str):
    pad = "." * max(2, _W - len(label))
    print(f"  {label} {pad} {_YEL}SKIP{_RST}  ({reason})")


# ── synthetic recording ────────────────────────────────────────────────────────
def _make_synthetic_recording() -> MEARecording:
    """
    Generate a 120-second, 20-channel synthetic recording.

    Design
    ------
    - Background Poisson firing at 8 Hz per channel
    - Synchronised bursts every 2 s: 80% of channels fire a tight
      cluster of spikes within a 10 ms window
    This gives clear avalanche structure and enough events for all
    criticality metrics to run without error.

    Amplitudes
    ----------
    Every spike also gets a synthetic peak-to-peak amplitude (µV), drawn
    around 50 µV. Channels 17-19 are seeded with one 500 µV spike each so
    the amplitude filter has a known-bad set to find.
    """
    rng      = np.random.default_rng(42)   # fixed seed — reproducible
    n_ch     = 20
    duration = 120.0
    sr       = 20_000.0

    spike_times: dict = {}
    for ch in range(n_ch):
        # background
        n_bg = rng.poisson(8 * duration)
        bg   = rng.uniform(0.5, duration - 0.5, n_bg)

        # synchronised bursts every 2 s
        centers = np.arange(1.0, duration, 2.0)
        burst   = []
        for c in centers:
            if rng.random() < 0.80:
                n_spk = rng.poisson(4) + 1
                t     = rng.normal(c, 0.010, n_spk)
                burst.extend(t.tolist())

        all_t = np.concatenate([bg, burst]) if burst else bg
        all_t = np.sort(all_t)
        all_t = all_t[(all_t >= 0) & (all_t < duration)]
        spike_times[ch] = all_t

    # per-spike peak-to-peak amplitudes (µV); ch 17-19 carry one artefact spike
    ARTEFACT_CH = (17, 18, 19)
    spike_amplitudes: dict = {}
    for ch, t in spike_times.items():
        amps = rng.normal(50.0, 8.0, t.size).astype(np.float32)
        if ch in ARTEFACT_CH and amps.size:
            amps[0] = 500.0
        spike_amplitudes[ch] = amps

    return MEARecording(
        sampling_rate = sr,
        n_channels    = n_ch,
        channel_ids   = list(range(n_ch)),
        spike_times   = spike_times,
        duration      = duration,
        metadata      = {"source": "CASCADE self-test (synthetic)"},
        spike_amplitudes = spike_amplitudes,
    )


def _all_spikes(rec: MEARecording) -> np.ndarray:
    """Concatenate and sort all spike times into a single 1-D array."""
    return np.sort(np.concatenate(list(rec.spike_times.values())))


# ── individual tests ───────────────────────────────────────────────────────────
def _test_recording(rec: MEARecording):
    assert rec.n_channels == 20
    assert rec.duration == 120.0
    assert len(rec.spike_times) == 20
    total = sum(len(v) for v in rec.spike_times.values())
    assert total > 500, f"Too few spikes: {total}"


def _test_spike_filtering(rec: MEARecording):
    from cascade.core.spikes import filter_active_electrodes
    active = filter_active_electrodes(rec.spike_times, ae_per_min=5, duration_s=rec.duration)
    assert len(active) > 0, "No active electrodes after filter"


def _test_amplitude_filter(rec: MEARecording):
    from cascade.core.spikes import filter_amplitude_electrodes
    # threshold between the 50 µV population and the 500 µV artefacts
    kept, n_removed, removed = filter_amplitude_electrodes(
        rec.spike_times, rec.spike_amplitudes, max_amplitude_uv=200.0)
    assert sorted(removed) == [17, 18, 19], f"Unexpected removals: {removed}"
    assert n_removed == 3, f"Expected 3 removals, got {n_removed}"
    assert len(kept) == len(rec.spike_times) - 3
    assert not any(ch in kept for ch in (17, 18, 19))
    # threshold of 0 must disable the filter entirely
    kept0, n0, _ = filter_amplitude_electrodes(
        rec.spike_times, rec.spike_amplitudes, max_amplitude_uv=0.0)
    assert n0 == 0 and len(kept0) == len(rec.spike_times), "0 did not disable filter"
    # a recording with no amplitude data must pass through untouched
    kept_none, n_none, _ = filter_amplitude_electrodes(
        rec.spike_times, None, max_amplitude_uv=200.0)
    assert n_none == 0 and len(kept_none) == len(rec.spike_times)


def _test_isi(rec: MEARecording):
    from cascade.core.spikes import compute_isi_metrics
    df = compute_isi_metrics(rec.spike_times, sampling_rate=rec.sampling_rate)
    assert df is not None and len(df) > 0


def _test_spike_binning(rec: MEARecording) -> object:
    from cascade.core.spikes import bin_spikes
    rate_df = bin_spikes(rec.spike_times, duration_s=rec.duration, bin_sec=1.0)
    assert len(rate_df) > 0, "bin_spikes returned empty DataFrame"
    return rate_df


def _test_bursts(rec: MEARecording):
    from cascade.core.bursts import compute_burst_metrics
    all_t = _all_spikes(rec)
    burst_df, summary = compute_burst_metrics(all_t, min_spikes=5, max_gap_s=0.1)
    assert summary is not None


def _test_active_area(rec: MEARecording):
    from cascade.core.bursts import compute_active_area
    result = compute_active_area(rec.spike_times, rec.duration, window_s=10.0)
    assert result is not None


def _test_avalanches(rec: MEARecording):
    from cascade.core.avalanche import infer_dt, compute_avalanches
    all_t = _all_spikes(rec)
    dt    = infer_dt(all_t)
    assert 0.001 <= dt <= 0.5, f"Unexpected dt: {dt:.4f} s"
    ava_df = compute_avalanches(all_t, dt=dt, perc=0.0)
    ava_df = ava_df[ava_df["avalanche_size"] >= 2]
    assert len(ava_df) > 10, f"Too few avalanches: {len(ava_df)}"
    return ava_df, dt, all_t


def _test_powerlaw_normalization():
    """
    The power-law MLE must recover a known exponent from data drawn out of a
    discrete truncated power law. This is a regression guard for the
    normalisation bug present in CASCADE <= 0.9.0 (and in Habibollahi 2022
    tplfit.py), where the normalising sum ran over only the observed unique
    values and biased the exponent downward at low avalanche counts.
    """
    from cascade.core.criticality import tplfit

    rng   = np.random.default_rng(7)
    TRUE  = 1.5
    xs    = np.arange(1, 301, dtype=np.float64)
    probs = xs ** (-TRUE)
    probs = probs / probs.sum()

    # n = 421 matches the smallest benchmark recording, where the old bias was
    # about -0.2; the corrected fit must land close to the true exponent.
    errs_fixed, errs_legacy = [], []
    for _ in range(5):
        data = rng.choice(xs, size=421, p=probs)
        errs_fixed.append(tplfit(data, limit=12)[0] - TRUE)
        errs_legacy.append(tplfit(data, limit=12, legacy_normalization=True)[0] - TRUE)

    mean_fixed  = float(np.mean(errs_fixed))
    mean_legacy = float(np.mean(errs_legacy))
    assert abs(mean_fixed) < 0.10, (
        f"Corrected fit is biased: mean error {mean_fixed:+.4f} (expected |err| < 0.10)")
    # the legacy path must still reproduce the OLD biased behaviour, so that
    # previously published numbers stay reproducible
    assert mean_legacy < -0.10, (
        f"legacy_normalization no longer reproduces the old bias "
        f"(mean error {mean_legacy:+.4f}, expected < -0.10)")

    # a fixed xmin (--pl_xmin_S) must give exactly the search's fit at that xmin
    for leg in (False, True):
        a_auto, x_auto, ks_auto, _ = tplfit(data, limit=12, legacy_normalization=leg)
        a_fix, x_fix, ks_fix, _ = tplfit(data, legacy_normalization=leg, xmin=x_auto)
        assert (x_fix, a_fix, ks_fix) == (x_auto, a_auto, ks_auto), (
            f"fixed xmin={x_auto} disagrees with the search (legacy={leg})")


def _test_synchrony(rec: MEARecording):
    """
    Synchrony metrics against cases whose answers are fixed by the definitions.

    Kreuz et al. (2013): the SPIKE-distance is zero only for identical spike trains
    and is time-scale invariant, so stretching both trains must not change it.
    Paiva et al. (2010) Eq 9: the Synchrony Index is bounded in [0,1] and equals 1
    for identical trains.
    """
    from cascade.core.synchrony import (spike_distance, synchrony_index_pair,
                                        pooled_cross_correlogram,
                                        area_under_correlogram,
                                        full_width_half_height,
                                        compute_synchrony)
    T = 10.0
    t = np.arange(0.5, 10.0, 0.5)

    # SPIKE-distance: zero only for identical trains
    assert abs(spike_distance(t, t, T)) < 1e-12, "D_S must be 0 for identical trains"
    d_small = spike_distance(t, t + 0.05, T)
    d_large = spike_distance(t, t + 0.25, T)
    assert 0.0 < d_small < d_large <= 1.0, \
        f"D_S must grow with mismatch and stay in [0,1]: {d_small}, {d_large}"

    # time-scale invariance (Kreuz 2013): stretching both trains changes nothing
    rng = np.random.default_rng(11)
    a = np.sort(rng.uniform(0, T, 25)); b = np.sort(rng.uniform(0, T, 25))
    d1 = spike_distance(a, b, T)
    d3 = spike_distance(a * 3.0, b * 3.0, T * 3.0)
    assert abs(d1 - d3) < 1e-9, f"D_S not time-scale invariant: {d1} vs {d3}"

    # Synchrony Index: exactly 1 for identical, decreasing with offset, in [0,1]
    sg = 0.020
    assert abs(synchrony_index_pair(t, t, sg) - 1.0) < 1e-12, \
        "Synchrony Index must be 1 for identical trains"
    r_near = synchrony_index_pair(t, t + 0.005, sg)
    r_far = synchrony_index_pair(t, t + 0.050, sg)
    assert 0.0 <= r_far < r_near <= 1.0, \
        f"Synchrony Index must decrease with offset and stay in [0,1]: {r_near}, {r_far}"

    # Pooled cross-correlogram: FWHH must track the injected jitter.
    # For Gaussian jitter of sd j per train the pooled peak has sd j*sqrt(2),
    # so FWHH = 2.355 * j * sqrt(2).
    ev = np.arange(0.5, 60.0, 0.5)
    for j in (0.002, 0.005):
        st = {i: np.sort(ev + rng.normal(0, j, ev.size)) for i in range(10)}
        cc = pooled_cross_correlogram(st, 60.0, bin_s=0.001, max_lag_s=0.3)
        w = full_width_half_height(cc["lags_s"], cc["cc_norm"])
        pred = 2.355 * j * np.sqrt(2)
        assert abs(w - pred) < 0.35 * pred, \
            f"FWHH {w:.4f} s is far from the predicted {pred:.4f} s for jitter {j}"

    # synchronous activity must give a larger normalised area than independent
    sync_st = {i: np.sort(ev + rng.normal(0, 0.002, ev.size)) for i in range(10)}
    indep_st = {i: np.sort(rng.uniform(0, 60.0, ev.size)) for i in range(10)}
    a_sync = area_under_correlogram(
        *[pooled_cross_correlogram(sync_st, 60.0, 0.001, 0.3)[k]
          for k in ("lags_s", "cc_norm")], 0.020)
    a_indep = area_under_correlogram(
        *[pooled_cross_correlogram(indep_st, 60.0, 0.001, 0.3)[k]
          for k in ("lags_s", "cc_norm")], 0.020)
    assert a_sync > a_indep, \
        f"synchronous area {a_sync} should exceed independent {a_indep}"

    # When PySpike (the reference implementation, by an author of the 2013 paper)
    # is installed, our SPIKE-distance must match it exactly, since we delegate.
    from cascade.core.synchrony import HAVE_PYSPIKE
    if HAVE_PYSPIKE:
        import pyspike as _spk
        for a_, b_ in [(t, t.copy()), (t, t + 0.05), (t, t + 0.25),
                       (np.sort(rng.uniform(0, T, 30)), np.sort(rng.uniform(0, T, 30)))]:
            mine = spike_distance(a_, b_, T)
            ref = _spk.spike_distance(_spk.SpikeTrain(a_, [0.0, T]),
                                      _spk.SpikeTrain(b_, [0.0, T]))
            assert abs(mine - ref) < 1e-12, \
                f"SPIKE-distance disagrees with PySpike: {mine} vs {ref}"

    # the assembly runs on the synthetic recording and 'none' is a clean no-op
    m, cc = compute_synchrony(rec.spike_times, rec.duration, metrics="all",
                              max_channels=8)
    for key in ("area_xcorr_norm", "kreuz_metric", "synchrony_index"):
        assert key in m, f"missing {key} from compute_synchrony(metrics='all')"
    m_none, cc_none = compute_synchrony(rec.spike_times, rec.duration,
                                        metrics="none")
    assert cc_none is None and "area_xcorr_norm" not in m_none


def _test_criticality(ava_df, dt: float, all_t: np.ndarray):
    from cascade.core.criticality import compute_exponents, compute_branching_ratio
    from cascade.core.avalanche import counts_per_dt

    sizes = ava_df["avalanche_size"].values
    durs  = ava_df["duration_bins"].values
    result = compute_exponents(sizes, durs)
    assert result is not None

    counts, _ = counts_per_dt(all_t, dt)
    if counts.size > 10:
        br = compute_branching_ratio(counts)
        assert br is not None


def _test_connectivity(rate_df):
    from cascade.core.connectivity import compute_connectivity
    C, ch_ids = compute_connectivity(rate_df)
    assert C is not None or ch_ids is not None  # None is a valid result for sparse data


def _test_entropy(rec: MEARecording):
    from cascade.core.entropy import compute_entropy
    result = compute_entropy(rec.spike_times, rec.duration, window_s=10.0)
    assert result is not None


# ── orchestrator ───────────────────────────────────────────────────────────────
def run_self_test() -> int:
    """
    Run all self-tests and print a formatted report.
    Returns 0 if all passed, 1 if any failed.
    """
    from cascade import __version__

    sep = "=" * (_W + 12)
    print()
    print(f"{_BOLD}{_CYN}CASCADE v{__version__} — Built-in Self-Test{_RST}")
    print(sep)
    print(f"  Synthetic recording: 20 ch · 120 s · Poisson + bursts · seed 42")
    print(sep)

    passed = 0
    failed = 0

    # helper to run one test step
    def run(label, fn, *a, **kw):
        nonlocal passed, failed
        try:
            result = fn(*a, **kw)
            _ok(label)
            passed += 1
            return result
        except Exception as e:
            _fail(label, e)
            failed += 1
            return None

    # 1 — build recording
    rec = run("Synthetic MEA recording", _make_synthetic_recording)
    if rec is None:
        print(f"\n  {_RED}Cannot continue without a valid recording.{_RST}\n")
        _summary(passed, failed)
        return 1

    # 2 — spike filtering
    run("Spike filtering (active-electrode filter)", _test_spike_filtering, rec)

    # 3 — amplitude filter
    run("Spike filtering (amplitude filter)", _test_amplitude_filter, rec)

    # 4 — ISI metrics
    run("ISI metrics", _test_isi, rec)

    # 5 — spike binning (returned for connectivity test)
    rate_df = run("Spike binning", _test_spike_binning, rec)

    # 6 — burst detection
    run("Network burst detection", _test_bursts, rec)

    # 7 — active area
    run("Active area", _test_active_area, rec)

    # 8 — avalanche detection (returned for criticality test)
    ava_result = run("Avalanche detection", _test_avalanches, rec)

    # 9 — synchrony
    run("Synchrony metrics (xcorr, Kreuz, Paiva)", _test_synchrony, rec)

    # 10 — power-law normalisation (regression guard)
    run("Power-law normalisation (known exponent)", _test_powerlaw_normalization)

    # 11 — criticality
    if ava_result is None:
        _skip("Criticality metrics (tau, alpha, DCC, BR)", "avalanche detection failed")
    else:
        ava_df, dt, all_t = ava_result
        run("Criticality metrics (tau, alpha, DCC, BR)", _test_criticality, ava_df, dt, all_t)

    # 12 — connectivity
    if rate_df is None:
        _skip("Functional connectivity", "spike binning failed")
    else:
        run("Functional connectivity", _test_connectivity, rate_df)

    # 13 — entropy
    run("Spatial entropy", _test_entropy, rec)

    print(sep)
    _summary(passed, failed)
    return 0 if failed == 0 else 1


def _summary(passed: int, failed: int):
    total = passed + failed
    if failed == 0:
        print(f"  {_GRN}{_BOLD}All {total} tests passed.{_RST} "
              f"CASCADE is installed correctly.\n")
    else:
        print(f"  {_RED}{_BOLD}{failed} of {total} tests failed.{_RST} "
              f"Check the errors above.\n")
