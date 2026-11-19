"""
core/spikes.py
--------------
Spike-level processing: filtering, binning, ISI metrics.
All functions take plain numpy arrays or dicts — no file I/O here.
"""

import numpy as np
import pandas as pd
from typing import Dict


# ── Stimulation blanking ───────────────────────────────────────────────────────
def blank_after_stims(spike_times: Dict[int, np.ndarray], stim_times,
                      window_s: float, amplitudes=None):
    """
    Drop every spike that falls within [stim, stim + window_s) of any
    stimulation, on all channels (stimulus artefacts and evoked responses
    reach the whole array). Returns (spike_times, amplitudes, n_removed).
    """
    if window_s <= 0 or stim_times is None or len(stim_times) == 0:
        return spike_times, amplitudes, 0
    st = np.sort(np.asarray(stim_times, dtype=np.float64))
    out_t, out_a, removed = {}, ({} if amplitudes is not None else None), 0
    for ch, t in spike_times.items():
        i = np.searchsorted(st, t, side="right") - 1
        hit = (i >= 0) & (t - st[np.maximum(i, 0)] < window_s)
        removed += int(hit.sum())
        out_t[ch] = t[~hit]
        if out_a is not None and ch in amplitudes:
            out_a[ch] = amplitudes[ch][~hit]
    return out_t, out_a, removed


# ── Active-electrode filter ────────────────────────────────────────────────────
def filter_active_electrodes(spike_times: Dict[int, np.ndarray],
                              ae_per_min: float,
                              duration_s: float) -> Dict[int, np.ndarray]:
    """
    Remove channels that fired fewer than *ae_per_min* spikes per minute.

    Parameters
    ----------
    spike_times : dict  {channel_id: array of spike times in seconds}
    ae_per_min  : float – minimum spikes/minute threshold (0 = keep all)
    duration_s  : float – recording duration used to compute the threshold

    Returns
    -------
    Filtered spike_times dict (same format, subset of channels).
    """
    if ae_per_min <= 0:
        return spike_times

    required = int(np.ceil(ae_per_min * duration_s / 60.0))
    return {ch: st for ch, st in spike_times.items() if len(st) >= required}


# ── Noisy-electrode filter ─────────────────────────────────────────────────────
def filter_noisy_electrodes(spike_times: Dict[int, np.ndarray],
                             duration_s: float,
                             max_rate_per_min: float = 0.0,
                             max_rpv_pct: float = 0.0) -> tuple:
    """
    Remove channels that are likely recording electrical noise rather than
    real neurons.  Two independent criteria are available; both are off by
    default (value 0).

    Parameters
    ----------
    spike_times      : dict  {channel_id: array of spike times in seconds}
    duration_s       : float – recording duration in seconds
    max_rate_per_min : float – remove channels whose firing rate exceeds this
                               value (spikes / minute).  0 = off.
                               Recommended: 1000 (≈ 17 Hz) for in-vitro MEA.
    max_rpv_pct      : float – remove channels where more than this percentage
                               of inter-spike intervals are shorter than 2 ms
                               (refractory period violation).  0 = off.
                               Recommended: 1.0–2.0 % for in-vitro MEA.

    Returns
    -------
    filtered   : dict  – spike_times with noisy channels removed
    n_removed  : int   – number of channels discarded
    removed_ch : list  – channel ids that were removed
    """
    if max_rate_per_min <= 0 and max_rpv_pct <= 0:
        return spike_times, 0, []

    removed = []
    for ch, times in spike_times.items():
        t = np.asarray(times, dtype=np.float64)

        # ── criterion 1: maximum firing rate ──────────────────────────────────
        if max_rate_per_min > 0 and duration_s > 0:
            rate = len(t) / duration_s * 60.0   # spikes / minute
            if rate > max_rate_per_min:
                removed.append(ch)
                continue

        # ── criterion 2: refractory period violation ──────────────────────────
        if max_rpv_pct > 0 and t.size >= 2:
            isi  = np.diff(np.sort(t))
            rpv  = float(np.mean(isi < 0.002) * 100.0)  # % ISIs < 2 ms
            if rpv > max_rpv_pct:
                removed.append(ch)
                continue

    filtered = {ch: st for ch, st in spike_times.items() if ch not in removed}
    return filtered, len(removed), removed


# ── Amplitude filter ───────────────────────────────────────────────────────────
def filter_amplitude_electrodes(spike_times: Dict[int, np.ndarray],
                                 spike_amplitudes: Dict[int, np.ndarray],
                                 max_amplitude_uv: float) -> tuple:
    """
    Remove channels whose maximum single-spike peak-to-peak amplitude exceeds
    *max_amplitude_uv* µV.  Channels with any spike that large are likely
    recording mechanical or electrical artefacts rather than real neurons.

    Parameters
    ----------
    spike_times      : dict  {channel_id: array of spike times in seconds}
    spike_amplitudes : dict  {channel_id: array of peak-to-peak amplitudes in µV}
                             One value per spike, same order as spike_times.
    max_amplitude_uv : float – threshold in µV. Channels where ANY spike exceeds
                               this value are removed.  0 = off.
                               Recommended starting value: 200 µV.

    Returns
    -------
    filtered   : dict  – spike_times with offending channels removed
    n_removed  : int   – number of channels discarded
    removed_ch : list  – channel ids that were removed
    """
    if max_amplitude_uv <= 0 or not spike_amplitudes:
        return spike_times, 0, []

    removed = [
        ch for ch, amps in spike_amplitudes.items()
        if ch in spike_times and np.max(amps) > max_amplitude_uv
    ]
    filtered = {ch: st for ch, st in spike_times.items() if ch not in removed}
    return filtered, len(removed), removed


# ── Spike binning (firing rates) ───────────────────────────────────────────────
def bin_spikes(spike_times: Dict[int, np.ndarray],
               duration_s: float,
               bin_sec: float = 1.0) -> pd.DataFrame:
    """
    Compute per-channel firing rate in consecutive time bins.

    Parameters
    ----------
    spike_times : dict  {channel_id: array of times in seconds}
    duration_s  : float – recording duration (sets the bin edges)
    bin_sec     : float – bin width in seconds (default 1 s)

    Returns
    -------
    DataFrame with columns: bin_idx (int), time_s (float, bin centre),
                             ch_idx (int), rate_hz (float)
    """
    if not spike_times:
        return pd.DataFrame(columns=["bin_idx", "time_s", "ch_idx", "rate_hz"])

    nb     = max(1, int(np.ceil(duration_s / bin_sec)))
    edges  = np.linspace(0.0, nb * bin_sec, nb + 1)
    centres = (edges[:-1] + edges[1:]) / 2.0

    rows = []
    for ch, times in spike_times.items():
        counts, _ = np.histogram(times, bins=edges)
        rate       = counts / bin_sec
        rows.extend([(i, float(centres[i]), ch, float(r))
                     for i, r in enumerate(rate)])

    return pd.DataFrame(rows, columns=["bin_idx", "time_s", "ch_idx", "rate_hz"])


# ── ISI metrics ────────────────────────────────────────────────────────────────
def compute_isi_metrics(spike_times: Dict[int, np.ndarray],
                        sampling_rate: float) -> pd.DataFrame:
    """
    Per-channel inter-spike-interval (ISI) statistics.

    Columns returned
    ----------------
    ch_idx              – channel identifier
    isi_mean_s          – mean ISI in seconds
    median_isi_ms       – median ISI in milliseconds
    isi_cv              – coefficient of variation = std(ISI) / mean(ISI)
    mean_rate_hz        – mean firing rate (spikes / recording span)
    refractory_viol_frac– fraction of ISIs < 2 ms (refractory period violation)
    """
    cols = ["ch_idx", "isi_mean_s", "median_isi_ms", "isi_cv",
            "mean_rate_hz", "refractory_viol_frac"]

    if not spike_times:
        return pd.DataFrame(columns=cols)

    rows = []
    for ch, times in spike_times.items():
        t = np.sort(times)
        if t.size < 2:
            rows.append([ch, np.nan, np.nan, np.nan, 0.0, np.nan])
            continue

        isi_s  = np.diff(t)
        isi_mu = float(isi_s.mean())
        isi_cv = float(isi_s.std() / isi_mu) if isi_mu > 0 else 0.0
        rv     = float(np.mean(isi_s < 0.002))
        rows.append([
            ch,
            isi_mu,
            float(np.median(isi_s) * 1000.0),
            isi_cv,
            float(t.size / (t[-1] - t[0])) if t[-1] > t[0] else 0.0,
            rv,
        ])

    return pd.DataFrame(rows, columns=cols)
