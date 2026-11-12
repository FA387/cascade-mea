#!/usr/bin/env python
"""
ltpltd.py
---------
LTP / LTD classification from CASCADE firing-rate outputs.

Compares mean per-electrode firing rate in one or more MEASURE recordings
against a BASELINE recording. Each electrode is classified as:

  LTP       — rate increased by >= threshold %  (default 20 %)
  LTD       — rate decreased by >= threshold %  (default 20 %)
  NoChange  — within threshold

Input files are the rates_*.tsv files produced by the main CASCADE pipeline.
No raw MEA data file is needed.

Usage examples
--------------
# Single comparison (before → after stimulation):
cascade-ltpltd --baseline baseline_out/rates_Well_A1.tsv \\
                 --measure  stim_out/rates_Well_A1.tsv \\
                 --output ltpltd_results/

# Multiple stimulation conditions in one run:
cascade-ltpltd --baseline baseline_out/rates_Well_A1.tsv \\
                 --measure  stim1_out/rates_Well_A1.tsv \\
                            stim2_out/rates_Well_A1.tsv \\
                            stim3_out/rates_Well_A1.tsv \\
                 --output ltpltd_results/

# Change the LTP/LTD threshold to 30 %:
cascade-ltpltd --baseline baseline_out/rates_Well_A1.tsv \\
                 --measure  stim_out/rates_Well_A1.tsv \\
                 --threshold 30 --output ltpltd_results/
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from cascade import __version__
from cascade.utils.citation import add_citation_argument
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── constants ──────────────────────────────────────────────────────────────────
EPS = 1e-12        # prevents division-by-zero when baseline rate ≈ 0
CALL_ORDER = ["LTP", "NoChange", "LTD"]
CALL_COLOURS = {"LTP": "#e74c3c", "NoChange": "#95a5a6", "LTD": "#2980b9"}


# ── core computation ───────────────────────────────────────────────────────────
def per_electrode_mean(rates_csv: str) -> pd.DataFrame:
    """
    Read a rates_*.tsv produced by CASCADE and return mean firing rate
    per electrode.

    Parameters
    ----------
    rates_csv : path to the rates CSV (columns: bin_idx, time_s, ch_idx, rate_hz)

    Returns
    -------
    DataFrame with columns: ch_idx, mean_rate_hz
    """
    # Read TSV format (CASCADE v0.9.2 outputs TSV, not CSV)
    df = pd.read_csv(rates_csv, sep='\t')
    required = {"ch_idx", "rate_hz"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"{rates_csv} is missing columns {required - set(df.columns)}. "
            "Expected a rates_*.tsv file produced by CASCADE."
        )
    return (
        df.groupby("ch_idx")["rate_hz"]
        .mean()
        .reset_index()
        .rename(columns={"rate_hz": "mean_rate_hz"})
    )


def classify(pct: pd.Series, threshold: float) -> pd.Series:
    """Classify % change values into LTP / LTD / NoChange."""
    return pd.cut(
        pct,
        bins=[-np.inf, -threshold, threshold, np.inf],
        labels=["LTD", "NoChange", "LTP"],
    ).astype(str).replace("nan", pd.NA)


def compare(baseline_csv: str, measure_csv: str,
            threshold: float) -> pd.DataFrame:
    """
    Compare mean firing rates between baseline and measure.

    Returns
    -------
    DataFrame: ch_idx, baseline_mean, measure_mean, pct_change, call
    """
    base = per_electrode_mean(baseline_csv).rename(
        columns={"mean_rate_hz": "baseline_mean"}
    )
    meas = per_electrode_mean(measure_csv).rename(
        columns={"mean_rate_hz": "measure_mean"}
    )

    merged = base.merge(meas, on="ch_idx", how="outer")

    # safe denominator
    denom = merged["baseline_mean"].copy()
    denom[denom.abs() < EPS] = np.nan

    merged["pct_change"] = (
        (merged["measure_mean"] - merged["baseline_mean"]) / denom * 100
    )
    merged["call"] = classify(merged["pct_change"], threshold)

    return merged[["ch_idx", "baseline_mean", "measure_mean",
                   "pct_change", "call"]]


# ── summary helpers ────────────────────────────────────────────────────────────
def summarise(result_df: pd.DataFrame, label: str, threshold: float):
    """Print a tidy terminal summary."""
    total = result_df["call"].notna().sum()
    print(f"\n{'='*50}")
    print(f"  {label}")
    print(f"  Electrodes analysed : {len(result_df)}")
    print(f"  Threshold           : ±{threshold} %")
    print(f"{'='*50}")
    for call in CALL_ORDER:
        n = (result_df["call"] == call).sum()
        pct = 100 * n / total if total > 0 else 0.0
        print(f"  {call:<10}  {n:>4}  ({pct:.1f} %)")
    print()


# ── figures ────────────────────────────────────────────────────────────────────
def plot_pct_change(result_df: pd.DataFrame,
                   label: str, threshold: float,
                   out_png: str, save_svg: bool = False, save_jpg: bool = False):
    """
    Bar chart: % change per electrode, coloured by LTP/LTD/NoChange.
    A horizontal dashed line marks ±threshold.
    """
    df = result_df.dropna(subset=["pct_change", "call"]).sort_values("ch_idx")

    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.15 + 2), 4))

    colours = df["call"].map(CALL_COLOURS).fillna("#95a5a6")
    ax.bar(range(len(df)), df["pct_change"].values,
           color=colours, width=1.0, linewidth=0)

    ax.axhline( threshold, color="gray", lw=1, ls="--")
    ax.axhline(-threshold, color="gray", lw=1, ls="--")
    ax.axhline(0,          color="black", lw=0.8)

    ax.set_xlabel("Electrode (sorted by ch_idx)")
    ax.set_ylabel("% change vs baseline")
    ax.set_title(f"LTP / LTD — {label}  (threshold ±{threshold} %)")

    # legend patches
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=CALL_COLOURS["LTP"],      label="LTP"),
        Patch(facecolor=CALL_COLOURS["NoChange"],  label="NoChange"),
        Patch(facecolor=CALL_COLOURS["LTD"],       label="LTD"),
    ]
    ax.legend(handles=legend_elements, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    if save_svg:
        fig.savefig(out_png.replace(".png", ".svg"))
    if save_jpg:
        fig.savefig(out_png.replace(".png", ".jpg"), dpi=300)
    plt.close(fig)


def plot_summary_pie(result_df: pd.DataFrame,
                     label: str, threshold: float,
                     out_png: str, save_svg: bool = False, save_jpg: bool = False):
    """Pie chart showing proportion of LTP / LTD / NoChange electrodes."""
    counts = {c: (result_df["call"] == c).sum() for c in CALL_ORDER}
    # drop zero slices
    counts = {k: v for k, v in counts.items() if v > 0}

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.pie(
        counts.values(),
        labels=counts.keys(),
        colors=[CALL_COLOURS[k] for k in counts],
        autopct="%1.1f%%",
        startangle=90,
    )
    ax.set_title(f"LTP / LTD summary — {label}\n(threshold ±{threshold} %)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    if save_svg:
        fig.savefig(out_png.replace(".png", ".svg"))
    if save_jpg:
        fig.savefig(out_png.replace(".png", ".jpg"), dpi=300)
    plt.close(fig)


# ── main ───────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        prog="cascade-ltpltd",
        description=(
            "CASCADE LTP/LTD — classify per-electrode firing-rate changes\n"
            "relative to a baseline recording.\n\n"
            "Input files must be rates_*.tsv files produced by the main\n"
            "CASCADE pipeline (cascade --input ... --output ...).\n\n"
            "Example:\n"
            "  cascade-ltpltd --baseline baseline_out/rates_Well_A1.tsv \\\n"
            "                   --measure  stim_out/rates_Well_A1.tsv \\\n"
            "                   --output ltpltd_results/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )

    # Flags are grouped by workflow step (same order as README "All flags").

    g = ap.add_argument_group("1 Setup — input / output")
    g.add_argument("--baseline",
                    required=True,
                    metavar="TSV",
                    help="Path to the rates_*.tsv from the BASELINE CASCADE run. "
                         "This is the reference against which all MEASURE files "
                         "are compared.")
    g.add_argument("--measure",
                    required=True,
                    nargs="+",
                    metavar="TSV",
                    help="Path(s) to one or more rates_*.tsv files to compare "
                         "against the baseline. Separate multiple paths with "
                         "spaces.")
    g.add_argument("--output",
                    default="ltpltd_out",
                    metavar="DIR",
                    help="Output folder for TSV results and figures. "
                         "Created automatically if it does not exist. "
                         "(default: ltpltd_out)")

    g = ap.add_argument_group("3 Analysis — LTP/LTD classification")
    g.add_argument("--threshold",
                    type=float,
                    default=20.0,
                    metavar="PCT",
                    help="LTP/LTD threshold as a percentage change from baseline. "
                         "Electrodes with change >= +threshold are called LTP; "
                         "<= -threshold are called LTD. (default: 20)")

    g = ap.add_argument_group("4 Output — figures and terminal")
    g.add_argument("--svg",
                    action="store_true",
                    help="Also save figures as SVG.")
    g.add_argument("--jpg",
                    action="store_true",
                    help="Also save figures as JPG (300 dpi).")
    g.add_argument("--verbose",
                    action="store_true",
                    help="Print per-electrode table to terminal.")

    g = ap.add_argument_group("5 Misc — information (no analysis)")
    g.add_argument("--help", action="help",
                    help="Show this help message and exit.")
    g.add_argument("--version", action="version", version=f"CASCADE v{__version__}",
                    help="Show the CASCADE version and exit.")
    add_citation_argument(g)

    args = ap.parse_args()

    # ── validate inputs ────────────────────────────────────────────────────────
    if not os.path.isfile(args.baseline):
        print(f"ERROR: baseline file not found: {args.baseline}", file=sys.stderr)
        sys.exit(1)

    missing = [p for p in args.measure if not os.path.isfile(p)]
    if missing:
        for p in missing:
            print(f"ERROR: measure file not found: {p}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    # ── run each measure file ──────────────────────────────────────────────────
    all_results = []

    for measure_csv in args.measure:
        label = os.path.splitext(os.path.basename(
            os.path.dirname(os.path.abspath(measure_csv))
        ))[0]
        if not label:
            label = os.path.splitext(os.path.basename(measure_csv))[0]

        print(f"\nComparing: {measure_csv}")
        print(f"      vs.: {args.baseline}")

        result = compare(args.baseline, measure_csv, args.threshold)
        result.insert(0, "label", label)

        summarise(result, label, args.threshold)

        if args.verbose:
            print(result.to_string(index=False))

        # save per-condition TSV
        safe_label = label.replace("/", "_").replace("\\", "_")
        csv_path = os.path.join(args.output, f"ltpltd_{safe_label}.tsv")
        result.to_csv(csv_path, index=False, sep='\t')
        print(f"  Saved: {csv_path}")

        # save figures
        bar_path = os.path.join(args.output, f"ltpltd_{safe_label}_bar.png")
        pie_path = os.path.join(args.output, f"ltpltd_{safe_label}_pie.png")
        plot_pct_change(result, label, args.threshold, bar_path, args.svg, args.jpg)
        plot_summary_pie(result, label, args.threshold, pie_path, args.svg, args.jpg)
        print(f"  Saved: {bar_path}")
        print(f"  Saved: {pie_path}")

        all_results.append(result)

    # ── combined summary TSV (all conditions in one file) ─────────────────────
    if len(all_results) > 1:
        combined = pd.concat(all_results, ignore_index=True)
        combined_path = os.path.join(args.output, "ltpltd_all_conditions.tsv")
        combined.to_csv(combined_path, index=False, sep='\t')
        print(f"\nCombined TSV saved: {combined_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
