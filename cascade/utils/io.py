"""
utils/io.py
-----------
All file output: TSV tables, CSV tables, and the criticality text report.
The R script for LTP/LTD reads the TSV files produced here.
Keep column names stable — changing them breaks the R script.
"""

import io
import os
import re
import sys
import platform
import datetime
import numpy as np
import pandas as pd


def _strip_ansi(text: str) -> str:
    """Remove ANSI colour/style escape sequences from a string."""
    return re.compile(r'\x1b\[[0-9;]*[mABCDEFGHJKSTfhilmnprsu]').sub('', text)


# ── Tee — write to terminal and a buffer simultaneously ────────────────────────
class _Tee:
    """Duplicates writes to two streams: the real stdout and a capture buffer."""
    def __init__(self, original, buffer):
        self._original = original
        self._buffer   = buffer

    def write(self, data):
        self._original.write(data)
        self._buffer.write(data)

    def flush(self):
        self._original.flush()
        self._buffer.flush()

    def fileno(self):          # needed by some C-level code that checks the fd
        return self._original.fileno()

    def isatty(self):          # needed by code that checks for interactive terminal
        return self._original.isatty()


# ── Generic table savers ───────────────────────────────────────────────────────
def write_csv(df: pd.DataFrame, path: str):
    """Save a DataFrame as CSV (human-readable, used for intermediate data)."""
    df.to_csv(path, index=False)


def write_tsv(df: pd.DataFrame, path: str):
    """
    Save a DataFrame as TSV (tab-separated).
    The R LTP/LTD script reads these — do NOT change column names here.
    """
    df.to_csv(path, index=False, sep="\t")


# ── Criticality text report ────────────────────────────────────────────────────
def write_criticality_report(metrics: dict, out_path: str):
    """
    Write a human-readable criticality report to a .txt file.

    Expected keys in metrics
    ------------------------
    well        : str   – recording identifier
    dt          : float – avalanche bin width in seconds
    ava_size_thr: int   – minimum avalanche size threshold
    ava_dur_thr : int   – minimum avalanche duration threshold (bins)
    perc        : float – population activity percentile threshold (0 = off)
    tau         : float – size power-law exponent  P(S) ~ S^-tau
    alpha       : float – duration power-law exponent  P(T) ~ T^-alpha
    xmin_s      : int   – optimal xmin for size fit
    xmin_d      : int   – optimal xmin for duration fit
    beta_pred   : float – predicted scaling exponent = (alpha-1)/(tau-1)
    beta_fit    : float – empirical scaling exponent from <S>(T) ~ T^beta
    dcc         : float – deviation from criticality = |beta_pred - beta_fit|
    br          : float – branching ratio (Wilting & Priesemann 2018)
    br_intercept: float – intercept of the BR linear fit
    sc_beta     : float – shape-collapse optimal beta
    sce         : float – |sc_beta - beta_pred|
    sc_minerr   : float – minimum shape-collapse variance metric
    auto_fit    : dict  – optional; the automatic-xmin, untrimmed fit, present
                          when --pl_xmin_S/D or --pl_exclude changed the fit
    auto_fit_why: str   – "xmin" or "exclude": which of the two did
    xmax_s, xmax_d, n_fit_s, n_fit_d, n_total – upper cut-offs and avalanches
                          used, present only with --pl_exclude
    """
    lines = []
    w = metrics.get("well", "unknown")
    lines.append(f"CRITICALITY REPORT — {w}")
    lines.append("=" * 60)
    lines.append(f"Avalanche bin dt       : {metrics.get('dt', float('nan')):.6f} s"
                 f"  ({metrics.get('dt', float('nan'))*1000:.2f} ms)")
    lines.append(f"Avalanche threshold    : size >= {metrics.get('ava_size_thr', 2)}"
                 f"  |  duration >= {metrics.get('ava_dur_thr', 1)} bin(s)")
    lines.append(f"Population perc thresh : {metrics.get('perc', 0.0)*100:.1f}%"
                 "  (0 % = off, use raw spike counts)")
    lines.append("")
    lines.append("Power-law exponents  (truncated MLE, Habibollahi 2022 method)")
    lines.append(f"  tau   (size exponent)     : {metrics.get('tau',   float('nan')):.6f}"
                 f"   xmin_S = {metrics.get('xmin_s', 'n/a')}")
    lines.append(f"  alpha (duration exponent) : {metrics.get('alpha', float('nan')):.6f}"
                 f"   xmin_D = {metrics.get('xmin_d', 'n/a')}")
    if metrics.get("xmax_s") is not None:
        n_tot = metrics.get("n_total")
        lines.append(f"  ** UPPER TAIL TRIMMED (--pl_exclude, Habibollahi 2022 EXCLUDE)")
        lines.append(f"     size fit range     : [{metrics.get('xmin_s')}, {metrics['xmax_s']:g}]"
                     f"   {metrics.get('n_fit_s')}/{n_tot} avalanches used")
        lines.append(f"     duration fit range : [{metrics.get('xmin_d')}, {metrics['xmax_d']:g}]"
                     f"   {metrics.get('n_fit_d')}/{n_tot} avalanches used")
    auto = metrics.get("auto_fit")
    if auto:
        if metrics.get("auto_fit_why") == "exclude":
            lines.append("  ** Untrimmed fit (no --pl_exclude) for comparison:")
        else:
            lines.append("  ** xmin SET MANUALLY (--pl_xmin_S/--pl_xmin_D). Automatic fit for comparison:")
        lines.append(f"     tau = {auto.get('tau', float('nan')):.6f} (xmin_S = {auto.get('xmin_s')})"
                     f"   alpha = {auto.get('alpha', float('nan')):.6f} (xmin_D = {auto.get('xmin_d')})"
                     f"   DCC = {auto.get('dcc', float('nan')):.6f}")
    lines.append("")
    lines.append("Scaling relation + DCC  (Beggs & Plenz 2003 / Habibollahi 2022)")
    lines.append(f"  beta_pred = (alpha-1)/(tau-1) : {metrics.get('beta_pred', float('nan')):.6f}")
    lines.append(f"  beta_fit  from <S>(T) ~ T^beta : {metrics.get('beta_fit',  float('nan')):.6f}")
    lines.append(f"  DCC = |beta_pred - beta_fit|   : {metrics.get('dcc',       float('nan')):.6f}")
    lines.append("")
    lines.append("Branching Ratio  (Wilting & Priesemann 2018)")
    lines.append(f"  BR (slope)     : {metrics.get('br',           float('nan')):.6f}")
    lines.append(f"  h  (intercept) : {metrics.get('br_intercept', float('nan')):.6f}")
    lines.append("")
    lines.append("Shape Collapse  (Sethna 2001 / Friedman 2012)")
    lines.append(f"  beta_opt (best collapse)         : {metrics.get('sc_beta',   float('nan')):.6f}")
    lines.append(f"  SCE = |beta_opt - beta_pred|      : {metrics.get('sce',       float('nan')):.6f}")
    lines.append(f"  collapse metric var(F)/(max-min)²: {metrics.get('sc_minerr', float('nan')):.3e}")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ── Run log ────────────────────────────────────────────────────────────────────
def write_run_log(out_dir: str, cmd: str, status: str,
                  error: str = None, output: str = None):
    """
    Append a run entry to <out_dir>/cascade_run.log.

    Each entry records:
      - Timestamp (local time)
      - Full command line used
      - Result: SUCCESS or FAILED (with error message if applicable)
      - Full terminal output produced during the run

    Parameters
    ----------
    out_dir : str   – output folder (same as --out)
    cmd     : str   – full command line string
    status  : str   – 'SUCCESS' or 'FAILED'
    error   : str   – exception message, only written on failure
    output  : str   – captured terminal output
    """
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, "cascade_run.log")
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        "=" * 70,
        f"Timestamp  : {timestamp}",
        f"Command    : {cmd}",
        f"Python     : {sys.version.split()[0]}",
        f"Platform   : {platform.platform()}",
        f"Working dir: {os.getcwd()}",
        f"Result     : {status}",
    ]
    if error:
        lines.append(f"Error      : {error}")
    if output and output.strip():
        lines.append("Output     :")
        lines.append(_strip_ansi(output).rstrip())
    lines.append("")

    with open(log_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
