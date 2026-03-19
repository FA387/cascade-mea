"""
core/connectivity.py
--------------------
Pearson correlation-based functional connectivity.
"""

import numpy as np
import pandas as pd


def compute_connectivity(rate_df: pd.DataFrame,
                         zscore: bool = True):
    """
    Compute the pairwise Pearson correlation matrix from per-channel firing rates.

    Parameters
    ----------
    rate_df : DataFrame with columns bin_idx, time_s, ch_idx, rate_hz
              (output of core.spikes.bin_spikes)
    zscore  : if True, z-score each channel's rate before correlating

    Returns
    -------
    C       : (N×N) float array – correlation matrix; None if insufficient data
    ch_ids  : 1-D int array – channel IDs (row/column order matches C)
    """
    if rate_df.empty:
        return None, None

    piv = (rate_df
           .pivot(index="ch_idx", columns="time_s", values="rate_hz")
           .fillna(0.0))
    X = piv.to_numpy(dtype=np.float64)

    if zscore:
        mu  = X.mean(axis=1, keepdims=True)
        sd  = X.std(axis=1, keepdims=True) + 1e-9
        X   = (X - mu) / sd

    C = np.corrcoef(X)
    return C, piv.index.to_numpy()


def summarize_connectivity(C: np.ndarray,
                           ch_ids: np.ndarray,
                           threshold: float = 0.3):
    """
    Threshold the correlation matrix to produce a binary graph and
    compute basic network statistics.

    Parameters
    ----------
    C         : (N×N) correlation matrix
    ch_ids    : 1-D array of channel IDs (same order as C rows/cols)
    threshold : |r| ≥ threshold → edge exists

    Returns
    -------
    deg_df  : DataFrame – {ch_idx, degree}
    summary : dict      – {edge_density, n_edges, n_nodes}
    """
    if C is None:
        return pd.DataFrame(), dict(edge_density=np.nan, n_edges=0, n_nodes=0)

    n      = C.shape[0]
    mask   = np.triu(np.ones_like(C, dtype=bool), k=1)
    n_edges = int((np.abs(C)[mask] >= threshold).sum())
    possible = n * (n - 1) // 2

    A = (np.abs(C) >= threshold).astype(int)
    np.fill_diagonal(A, 0)
    deg = A.sum(axis=1)
    deg_df = pd.DataFrame({"ch_idx": ch_ids, "degree": deg})

    summary = dict(
        edge_density = float(n_edges / possible) if possible > 0 else np.nan,
        n_edges      = n_edges,
        n_nodes      = n,
    )
    return deg_df, summary
