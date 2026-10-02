"""
core/plots.py
-------------
All figure generation.  Every function receives processed data (DataFrames,
arrays, dicts) and a file path — no raw data loading here.

Every function calls save_fig() from utils.common so PNG + optional SVG
are handled consistently.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # headless — no display required
import matplotlib.pyplot as plt
from typing import Dict

from cascade.utils.common import save_fig

# ── global plot style ──────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "figure.dpi"    : 300,
    "savefig.dpi"   : 300,
    "font.size"     : 10,
    "font.family"   : "DejaVu Sans",
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 1.0,
})


# ── Raster plot ────────────────────────────────────────────────────────────────
def plot_raster(spike_times: Dict[int, np.ndarray],
                out: str, title: str, save_svg: bool = False, save_jpg: bool = False,
                n_channels: int = None):
    fig, ax = plt.subplots(figsize=(8, 5))
    if not spike_times:
        ax.text(0.5, 0.5, "No spikes", ha="center", va="center")
    else:
        for ch, times in spike_times.items():
            ax.scatter(times, np.full(len(times), ch), s=1, c="k", linewidths=0)
        if n_channels is not None:
            ax.set_ylim(-1, n_channels)
        else:
            ax.set_ylim(-1, max(spike_times.keys()) + 1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Electrode")
    ax.set_title(title)
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Rate heatmap ───────────────────────────────────────────────────────────────
def plot_rate_heatmap(rate_df: pd.DataFrame,
                      out: str, title: str, save_svg: bool = False, save_jpg: bool = False,
                      n_channels: int = None):
    fig, ax = plt.subplots(figsize=(8, 5))
    if rate_df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
    else:
        piv = rate_df.pivot(index="ch_idx", columns="bin_idx",
                            values="rate_hz").fillna(0)
        n_bins = int(piv.columns.max()) + 1
        try:
            max_ch = int(piv.index.max())
            n_rows = max(n_channels, max_ch + 1) if n_channels is not None else max_ch + 1
            full = np.zeros((n_rows, n_bins))
            for ch_idx in piv.index:
                for bin_idx in piv.columns:
                    full[int(ch_idx), int(bin_idx)] = piv.loc[ch_idx, bin_idx]
            im = ax.imshow(full, aspect="auto", origin="lower")
            fig.colorbar(im, ax=ax, label="Hz")
        except (ValueError, TypeError):
            # string channel IDs (e.g. Axion 'B3_88') — use sequential rows
            ch_labels = list(piv.index)
            full = piv.to_numpy()
            im = ax.imshow(full, aspect="auto", origin="lower")
            fig.colorbar(im, ax=ax, label="Hz")
            step = max(1, len(ch_labels) // 10)
            ticks = list(range(0, len(ch_labels), step))
            ax.set_yticks(ticks)
            ax.set_yticklabels([ch_labels[i] for i in ticks], fontsize=7)
    ax.set_xlabel("Time bin")
    ax.set_ylabel("Electrode")
    ax.set_title(title)
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Population firing rate ─────────────────────────────────────────────────────
def plot_population_rate(rate_df: pd.DataFrame,
                         out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    fig, ax = plt.subplots(figsize=(8, 3))
    if not rate_df.empty:
        mean_rate = rate_df.groupby("time_s")["rate_hz"].mean()
        ax.plot(mean_rate.index, mean_rate.values)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean rate (Hz)")
    ax.set_title(f"Population firing rate — {label}")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Burst timeline ─────────────────────────────────────────────────────────────
def plot_burst_timeline(burst_df: pd.DataFrame,
                        duration_s: float,
                        out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    fig, ax = plt.subplots(figsize=(8, 1.2))
    if not burst_df.empty:
        for _, r in burst_df.iterrows():
            ax.axvspan(r["burst_start_s"], r["burst_end_s"],
                       color="tomato", alpha=0.7)
    ax.set_xlim(0, duration_s)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Time (s)")
    ax.set_yticks([])
    ax.set_title(f"Network bursts — {label}")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Active-area curve ──────────────────────────────────────────────────────────
def plot_active_area_curve(act_df: pd.DataFrame,
                           out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    fig, ax = plt.subplots(figsize=(8, 3))
    if not act_df.empty:
        ax.plot(act_df["window_start_s"],
                act_df["active_frac"] * 100, "-o", lw=1, ms=3)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Active electrodes (%)")
    ax.set_title(f"Active area — {label}")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Waveform templates (3Brain only) ──────────────────────────────────────────
def plot_waveforms(waveform_dict: dict,
                   out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    """
    waveform_dict: {ch_id: {"mean": array, "std": array, "t_ms": array}}
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    if not waveform_dict:
        ax.text(0.5, 0.5, "No waveform data", ha="center", va="center")
    else:
        for ch, d in waveform_dict.items():
            t  = d["t_ms"]
            mu = d["mean"]
            sd = d["std"]
            ax.plot(t, mu, label=f"ch{ch}")
            ax.fill_between(t, mu - sd, mu + sd, alpha=0.2)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude (ADC units)")
    ax.set_title(f"Mean spike waveforms — {label}")
    if len(waveform_dict) <= 12:
        ax.legend(fontsize=7, ncols=2)
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Correlation heatmap ────────────────────────────────────────────────────────
def plot_corr_heatmap(C: np.ndarray,
                      out: str, threshold: float, label: str,
                      save_svg: bool = False, save_jpg: bool = False,
                      channel_ids: list = None, n_channels: int = None):
    fig, ax = plt.subplots(figsize=(6, 5))
    if C is not None:
        display_C = C
        display_ids = list(channel_ids) if channel_ids is not None else None

        # expand to full electrode grid when channel IDs are integers
        if (channel_ids is not None and n_channels is not None):
            try:
                ids = [int(c) for c in channel_ids]
                grid = max(n_channels, max(ids) + 1)
                full = np.full((grid, grid), np.nan)
                for i, ci in enumerate(ids):
                    for j, cj in enumerate(ids):
                        full[ci, cj] = C[i, j]
                display_C = full
                display_ids = None  # ticks will be electrode indices directly
            except (ValueError, TypeError):
                pass  # string IDs — keep original Y×Y

        cmap = matplotlib.cm.coolwarm.copy()
        cmap.set_bad("lightgray")
        im = ax.imshow(display_C, vmin=-1, vmax=1, cmap=cmap, origin="lower")
        fig.colorbar(im, ax=ax, label="Pearson r")

        if display_ids is not None and len(display_ids) == display_C.shape[0]:
            n = len(display_ids)
            step = max(1, n // 10)
            ticks = list(range(0, n, step))
            ax.set_xticks(ticks)
            ax.set_xticklabels([display_ids[i] for i in ticks], rotation=45, ha="right")
            ax.set_yticks(ticks)
            ax.set_yticklabels([display_ids[i] for i in ticks])
    else:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
    ax.set_title(f"Connectivity — {label} (|r| ≥ {threshold})")
    ax.set_xlabel("Electrode")
    ax.set_ylabel("Electrode")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Degree histogram ───────────────────────────────────────────────────────────
def plot_degree_histogram(deg_df: pd.DataFrame,
                          out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    fig, ax = plt.subplots(figsize=(6, 3))
    if not deg_df.empty:
        ax.hist(deg_df["degree"], bins=30, color="steelblue", edgecolor="white")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Count")
    ax.set_title(f"Degree distribution — {label}")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Avalanche distributions ────────────────────────────────────────────────────
def plot_avalanche_distributions(av_df: pd.DataFrame,
                                 out_prefix: str, label: str,
                                 save_svg: bool = False, save_jpg: bool = False):
    """Creates two figures: <out_prefix>_size.png  and  <out_prefix>_duration.png"""
    for col, xlabel, title_word, suffix in [
        ("avalanche_size",  "Avalanche size (spikes)",   "size",     "_size"),
        ("duration_bins",   "Avalanche duration (bins)", "duration", "_duration"),
    ]:
        fig, ax = plt.subplots(figsize=(5, 4))
        if not av_df.empty:
            x = av_df[col].to_numpy()
            x = x[x > 0]
            if x.size > 0:
                lo  = max(1, int(x.min()))
                hi  = max(lo + 1, int(x.max()))
                bins = np.logspace(np.log10(lo), np.log10(hi), 20)
                hist, edges = np.histogram(x, bins=bins)
                centres = np.sqrt(edges[:-1] * edges[1:])
                ax.loglog(centres, np.maximum(hist, 1), "o", ms=4)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Count")
        ax.set_title(f"Avalanche {title_word} — {label}")
        fig.tight_layout()
        save_fig(fig, out_prefix + suffix + ".png", save_svg, save_jpg)
        plt.close(fig)


# ── Log-log linearisation (exponent visualisation) ────────────────────────────
def plot_loglog_linearization(data: np.ndarray,
                               out: str, title: str,
                               tau: float = None, xmin: float = None,
                               exclude_singletons: bool = False,
                               save_svg: bool = False, save_jpg: bool = False):
    """
    Plot empirical P(x) in log-log space with the fitted power-law line.

    Parameters
    ----------
    data : raw samples (e.g. all avalanche sizes)
    tau  : fitted exponent (for annotation)
    xmin : lower cut-off used in the fit (for annotation)
    """
    x = np.asarray(data, dtype=np.float64)
    x = x[np.isfinite(x) & (x > 0)]

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    if x.size == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        fig.tight_layout()
        save_fig(fig, out, save_svg, save_jpg)
        plt.close(fig)
        return

    vals, cnt = np.unique(x.astype(np.int64), return_counts=True)
    p_emp     = cnt.astype(np.float64) / cnt.sum()
    ax.scatter(np.log10(vals), np.log10(p_emp), s=12, label="empirical")

    if tau is not None and xmin is not None and np.isfinite(tau):
        fit_mask = vals >= xmin
        # --skip_rare_sizes / --skip_rare_durations: remove points where
        # the avalanche size/duration was observed only once (count=1).
        # These single-occurrence outliers can distort the visual fit line.
        if exclude_singletons:
            fit_mask = fit_mask & (cnt > 1)
        x_fit = vals[fit_mask]
        if x_fit.size >= 2:
            log_x = np.log10(x_fit.astype(np.float64))
            log_y = np.polyval(np.polyfit(log_x,
                                          np.log10(p_emp[fit_mask]), 1),
                               log_x)
            ax.plot(log_x, log_y, label=f"fit τ={tau:.3f}")

        ax.text(0.03, 0.97,
                f"xmin={xmin}\nMLE exponent={tau:.4f}",
                transform=ax.transAxes, va="top", fontsize=9)

    ax.set_xlabel("log₁₀(x)")
    ax.set_ylabel("log₁₀ P(x)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Shape collapse ─────────────────────────────────────────────────────────────
def plot_shape_collapse(details: dict,
                        out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    if not details or "mean_profiles" not in details \
            or not np.isfinite(details.get("beta_opt", np.nan)):
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No shape-collapse data", ha="center", va="center")
        fig.tight_layout()
        save_fig(fig, out, save_svg, save_jpg)
        plt.close(fig)
        return

    beta     = details["beta_opt"]
    gamma    = beta - 1.0
    x_common = details["x_common"]
    durations= details["durations"]
    profiles = details["mean_profiles"]

    fig, ax = plt.subplots(figsize=(6, 3.5))
    for D in durations:
        s  = profiles[D]
        F  = s / (float(D) ** gamma)
        x  = np.arange(D, dtype=np.float64) / float(D)
        Fi = np.interp(x_common, x, F)
        ax.plot(x_common, Fi, lw=1)

    ax.set_xlabel("Scaled time  t / D")
    ax.set_ylabel("Scaled shape  F(t/D)")
    ax.set_title(f"Shape collapse — {label}  (β = {beta:.3f})")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Entropy curve ──────────────────────────────────────────────────────────────
def plot_entropy_curve(ent_df: pd.DataFrame,
                       out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    fig, ax = plt.subplots(figsize=(8, 3))
    if not ent_df.empty:
        ax.plot(ent_df["window_start_s"], ent_df["H_norm"], "-o", lw=1, ms=3)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalised entropy")
    ax.set_title(f"Spatial entropy — {label}")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Epoch comparison ───────────────────────────────────────────────────────────
def plot_epoch_compare(stats: dict,
                       out: str, label: str, save_svg: bool = False, save_jpg: bool = False):
    """
    stats = {"A": {"pop_rate": float, "active_area": float},
             "B": {"pop_rate": float, "active_area": float}}
    """
    fig, ax = plt.subplots(figsize=(6, 3))
    labels = ["Pop. rate (Hz)", "Active area (%)"]
    A      = [stats["A"]["pop_rate"], stats["A"]["active_area"] * 100]
    B      = [stats["B"]["pop_rate"], stats["B"]["active_area"] * 100]
    x      = np.arange(len(labels))
    w      = 0.35
    ax.bar(x - w / 2, A, width=w, label="Epoch A")
    ax.bar(x + w / 2, B, width=w, label="Epoch B")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(f"Epoch comparison — {label}")
    ax.legend()
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)


# ── Synchrony ──────────────────────────────────────────────────────────────────
def plot_cross_correlogram(cc: dict,
                           out: str, label: str,
                           window_s: float = 0.020,
                           normalized: bool = True,
                           save_svg: bool = False, save_jpg: bool = False):
    """
    Well-wide pooled inter-electrode cross-correlogram.

    Follows the Axion Neural Metric Tool display: the correlogram with vertical
    grey lines marking the synchrony window, and the shaded area whose integral is
    the "area under cross-correlation" metric.

    cc : dict from core.synchrony.pooled_cross_correlogram
    """
    fig, ax = plt.subplots(figsize=(8, 3.2))
    if cc is not None:
        lags_ms = cc["lags_s"] * 1000.0
        y = cc["cc_norm"] if normalized else cc["cc"]
        ax.plot(lags_ms, y, "-", lw=1.0, color="#1f77b4")
        m = np.abs(cc["lags_s"]) <= window_s
        ax.fill_between(lags_ms[m], 0, y[m], color="#1f77b4", alpha=0.25)
        # synchrony-window guide lines, as Axion draws them
        for xv in (-window_s * 1000.0, window_s * 1000.0):
            ax.axvline(xv, color="grey", lw=0.9, ls="-")
        ax.set_xlim(lags_ms.min(), lags_ms.max())
    ax.set_xlabel("Phase lag (ms)")
    ax.set_ylabel("Normalised correlation" if normalized else "Correlation")
    kind = "normalised " if normalized else ""
    n_pairs = cc["n_pairs"] if cc else 0
    ax.set_title(f"Pooled {kind}cross-correlogram — {label} "
                 f"({n_pairs} pairs, window ±{window_s*1000:.0f} ms)")
    fig.tight_layout()
    save_fig(fig, out, save_svg, save_jpg)
    plt.close(fig)
