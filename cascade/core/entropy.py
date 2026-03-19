"""
core/entropy.py
---------------
Spatial entropy of population activity over time windows.

H(t) = -Σ p_i · log2(p_i)   where p_i = fraction of spikes on channel i
                               in the current window.

H_norm = H / log2(k)          where k = number of active channels.
H_norm ∈ [0, 1]:  1 = maximally distributed, 0 = all spikes on one channel.
"""

import numpy as np
import pandas as pd
from typing import Dict


def compute_entropy(spike_times: Dict[int, np.ndarray],
                    duration_s: float,
                    window_s: float = 10.0) -> pd.DataFrame:
    """
    Compute spatial entropy in consecutive time windows.

    Parameters
    ----------
    spike_times : dict {channel_id: array of spike times in seconds}
    duration_s  : total recording duration (sets window edges)
    window_s    : window width in seconds

    Returns
    -------
    DataFrame with columns:
        window_start_s   – window start (seconds)
        window_end_s     – window end (seconds)
        H_bits           – Shannon entropy in bits
        H_norm           – normalised entropy ∈ [0, 1]
    """
    cols = ["window_start_s", "window_end_s", "H_bits", "H_norm"]

    if not spike_times:
        return pd.DataFrame(columns=cols)

    all_ch  = list(spike_times.keys())
    n_total = len(all_ch)
    edges   = np.arange(0.0, duration_s + window_s, window_s)
    rows    = []

    for i in range(len(edges) - 1):
        w0, w1 = edges[i], edges[i + 1]
        counts = np.array([
            int(np.sum((st >= w0) & (st < w1)))
            for st in spike_times.values()
        ], dtype=np.float64)

        tot = counts.sum()
        if tot <= 0:
            rows.append([float(w0), float(w1), 0.0, 0.0])
            continue

        p  = counts / tot
        p  = p[p > 0]
        H  = float(-(p * np.log2(p)).sum())
        k  = int((counts > 0).sum())
        H_norm = float(H / np.log2(k)) if k > 1 else 0.0
        rows.append([float(w0), float(w1), H, H_norm])

    return pd.DataFrame(rows, columns=cols)
