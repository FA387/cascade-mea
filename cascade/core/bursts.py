"""
core/bursts.py
--------------
Network burst detection and active-area computation.
Works on flat spike-time arrays — no file I/O.
"""

import numpy as np
import pandas as pd
from typing import Dict


# ── Network burst detection ────────────────────────────────────────────────────
def detect_network_bursts(all_spike_times: np.ndarray,
                          min_spikes: int = 20,
                          max_gap_s: float = 0.1) -> list:
    """
    Detect network bursts from the population spike train.

    A burst is a run of consecutive spikes where no gap between
    adjacent spikes exceeds *max_gap_s*, containing at least *min_spikes* spikes.

    Parameters
    ----------
    all_spike_times : 1-D array of spike times in seconds (all channels combined)
    min_spikes      : minimum number of spikes to qualify as a burst
    max_gap_s       : maximum allowed gap between consecutive spikes (seconds)

    Returns
    -------
    List of (start_s, end_s) tuples.
    """
    t = np.sort(all_spike_times)
    if t.size < min_spikes:
        return []

    bursts = []
    start  = t[0]
    count  = 1

    for i in range(1, t.size):
        gap = t[i] - t[i - 1]
        if gap <= max_gap_s:
            count += 1
        else:
            if count >= min_spikes:
                bursts.append((float(start), float(t[i - 1])))
            start = t[i]
            count = 1

    if count >= min_spikes:
        bursts.append((float(start), float(t[-1])))

    return bursts


def compute_burst_metrics(all_spike_times: np.ndarray,
                          min_spikes: int = 20,
                          max_gap_s: float = 0.1):
    """
    Detect bursts and compute summary statistics.

    Returns
    -------
    burst_df : DataFrame  – columns: burst_start_s, burst_end_s, duration_s
    summary  : dict       – n_bursts, mean_duration_s, mean_IBI_s,
                            frac_time_in_bursts
    """
    empty_df = pd.DataFrame(columns=["burst_start_s", "burst_end_s", "duration_s"])

    if all_spike_times.size == 0:
        return empty_df, dict(n_bursts=0, mean_duration_s=np.nan,
                              mean_IBI_s=np.nan, frac_time_in_bursts=0.0,
                              burst_frequency_hz=0.0)

    bursts   = detect_network_bursts(all_spike_times, min_spikes, max_gap_s)
    burst_df = pd.DataFrame(
        [[s, e, e - s] for s, e in bursts],
        columns=["burst_start_s", "burst_end_s", "duration_s"]
    )

    if burst_df.empty:
        return burst_df, dict(n_bursts=0, mean_duration_s=np.nan,
                              mean_IBI_s=np.nan, frac_time_in_bursts=0.0,
                              burst_frequency_hz=0.0)

    ibi = np.diff(burst_df["burst_start_s"].to_numpy())
    total_burst_time = burst_df["duration_s"].sum()
    recording_span   = all_spike_times[-1] - all_spike_times[0] if all_spike_times.size > 1 else 1.0

    summary = dict(
        n_bursts             = len(burst_df),
        mean_duration_s      = float(burst_df["duration_s"].mean()),
        mean_IBI_s           = float(ibi.mean()) if ibi.size > 0 else np.nan,
        frac_time_in_bursts= float(total_burst_time / recording_span)
                               if recording_span > 0 else np.nan,
        burst_frequency_hz   = float(len(burst_df) / recording_span)
                               if recording_span > 0 else 0.0,
    )
    return burst_df, summary


# ── Active-area curve ──────────────────────────────────────────────────────────
def compute_active_area(spike_times: Dict[int, np.ndarray],
                        duration_s: float,
                        window_s: float = 10.0) -> pd.DataFrame:
    """
    Fraction of channels active in consecutive time windows.

    Parameters
    ----------
    spike_times : dict  {channel_id: array of times in seconds}
    duration_s  : float – total recording duration
    window_s    : float – window size in seconds

    Returns
    -------
    DataFrame with columns: window_start_s, window_end_s, active_frac
    """
    if not spike_times:
        return pd.DataFrame(columns=["window_start_s", "window_end_s", "active_frac"])

    all_ch    = list(spike_times.keys())
    n_total   = len(all_ch)
    edges     = np.arange(0.0, duration_s + window_s, window_s)

    rows = []
    for i in range(len(edges) - 1):
        w0, w1 = edges[i], edges[i + 1]
        n_active = sum(
            1 for ch, times in spike_times.items()
            if np.any((times >= w0) & (times < w1))
        )
        rows.append([float(w0), float(w1),
                     float(n_active / n_total) if n_total > 0 else 0.0])

    return pd.DataFrame(rows, columns=["window_start_s", "window_end_s", "active_frac"])
