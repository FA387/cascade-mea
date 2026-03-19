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

    return MEARecording(
        sampling_rate = sr,
        n_channels    = n_ch,
        channel_ids   = list(range(n_ch)),
        spike_times   = spike_times,
        duration      = duration,
        metadata      = {"source": "CASCADE self-test (synthetic)"},
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

    # 3 — ISI metrics
    run("ISI metrics", _test_isi, rec)

    # 4 — spike binning (returned for connectivity test)
    rate_df = run("Spike binning", _test_spike_binning, rec)

    # 5 — burst detection
    run("Network burst detection", _test_bursts, rec)

    # 6 — active area
    run("Active area", _test_active_area, rec)

    # 7 — avalanche detection (returned for criticality test)
    ava_result = run("Avalanche detection", _test_avalanches, rec)

    # 8 — criticality
    if ava_result is None:
        _skip("Criticality metrics (tau, alpha, DCC, BR)", "avalanche detection failed")
    else:
        ava_df, dt, all_t = ava_result
        run("Criticality metrics (tau, alpha, DCC, BR)", _test_criticality, ava_df, dt, all_t)

    # 9 — connectivity
    if rate_df is None:
        _skip("Functional connectivity", "spike binning failed")
    else:
        run("Functional connectivity", _test_connectivity, rate_df)

    # 10 — entropy
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
