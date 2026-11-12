#!/usr/bin/env python
"""
run_analysis.py
---------------
Unified MEA analysis pipeline.
Supports 3Brain (.bxr), Maxwell (.h5), MCS (.h5), Cortical Labs (.h5),
and Axion BioSystems (.csv) files.
Adding a new manufacturer = write one new loader in loaders/ only.

Usage examples
--------------
# 3Brain file:
python run_analysis.py --input recording.bxr --manufacturer 3brain --output results/

# Maxwell file:
python run_analysis.py --input data.raw.h5 --manufacturer maxwell --output results/

# CL1 file:
python run_analysis.py --input recording.h5 --manufacturer cl --output results/

# MCS file:
python run_analysis.py --input recording.h5 --manufacturer mcs --output results/

# Axion spike-list CSV:
python run_analysis.py --input recording_spike_list.csv --manufacturer axion --output results/

# With active-electrode filter and Habibollahi threshold:
python run_analysis.py --input recording.bxr --manufacturer 3brain --ae 5 --perc 0.25 --output results/

# With epoch comparison (pre vs post stimulation):
python run_analysis.py --input recording.bxr --manufacturer 3brain --epochA 0,60 --epochB 120,60 --output results/

Flags
-----
--input              path to data file (auto-searched if omitted)
--manufacturer      3brain | maxwell | cl | mcs | axion  (required)
--well              well/stream name (default: first available)
--start             start of analysis window in seconds (default: 0)
--dur               duration of analysis window in seconds (default: full recording)
--output               output folder (default: mea_out)
--bin               bin size for firing rates in seconds (default: 1.0)
--ae                active-electrode filter: keep channels with >= X spikes/min (default: 0 = off)
--max_rate          noisy-electrode filter: remove channels with > X spikes/min (default: 0 = off)
                    recommended: 1000 (≈ 17 Hz) for in-vitro MEA
--max_rpv           noisy-electrode filter: remove channels with > X% of ISIs < 2 ms (default: 0 = off)
                    recommended: 1.0–2.0% for in-vitro MEA
--ava_size          minimum avalanche size in spikes (default: 2)
--ava_dur           minimum avalanche duration in bins (default: 1 — no duration filter)
--perc              percentile threshold for avalanche detection (default: 0.0 = off)
                    0.0 = any spike counts; 0.25 = Habibollahi 2022 default
--aval_dt           avalanche bin width in seconds (default: inferred from mean ISI)
--burst_min_spikes  min spikes in a run to count as a network burst (default: 5)
--burst_max_gap     max gap between consecutive spikes within a burst in seconds (default: 0.1)
--burst_min         upper xmin limit for size power-law fit (default: 12)
--t_min             upper xmin limit for duration power-law fit (default: 4)
--conn_thr          |r| threshold for connectivity edges (default: 0.3)
--entropy_win       window size for entropy and active-area curves in seconds (default: 10)
--sc_Dmin           shape collapse: min avalanche duration in bins (default: 4)
--sc_Dmax           shape collapse: max avalanche duration in bins (default: 20)
--epochA            epoch A window: 'start,duration' in seconds  e.g. '0,60'
--epochB            epoch B window: 'start,duration' in seconds  e.g. '120,60'
--svg               also save all figures as SVG (default: PNG only)
--jpg               also save all figures as JPG (default: PNG only)
--verbose           print progress messages
"""

import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd

# ── imports from this package ──────────────────────────────────────────────────
from cascade.loaders.brainwave_loader import BrainwaveLoader
from cascade.loaders.maxwell_loader   import MaxwellLoader
from cascade.loaders.cl1_loader       import CL1Loader
from cascade.loaders.mcs_loader       import MCSLoader
from cascade.loaders.axion_loader     import AxionLoader
from cascade.core.spikes        import filter_active_electrodes, filter_noisy_electrodes, filter_amplitude_electrodes, bin_spikes, compute_isi_metrics
from cascade.core.bursts        import compute_burst_metrics, compute_active_area
from cascade.core.avalanche     import infer_dt, compute_avalanches, counts_per_dt
from cascade.core.criticality   import compute_exponents, compute_branching_ratio, shape_collapse
from cascade.core.connectivity  import compute_connectivity, summarize_connectivity
from cascade.core.entropy       import compute_entropy
from cascade.core.synchrony     import compute_synchrony, synchrony_to_dataframe
from cascade.core.plots         import (plot_raster, plot_rate_heatmap, plot_population_rate,
                                        plot_burst_timeline, plot_active_area_curve,
                                        plot_waveforms, plot_corr_heatmap, plot_degree_histogram,
                                        plot_avalanche_distributions, plot_loglog_linearization,
                                        plot_shape_collapse, plot_entropy_curve, plot_epoch_compare,
                                        plot_cross_correlogram)
from cascade.utils.common       import say, banner, step, CYAN, GRN, YEL, RED
from cascade.utils.io           import write_csv, write_tsv, write_criticality_report, write_run_log, _Tee


# ── manufacturer auto-detection ────────────────────────────────────────────────
def _peek_h5_manufacturer(path: str) -> str:
    """Peek inside an HDF5 file to distinguish CL1 / MCS / Maxwell."""
    try:
        import h5py
        with h5py.File(path, "r") as h5:
            sid = h5.attrs.get("system_id", b"")
            if isinstance(sid, (bytes, np.bytes_)):
                sid = sid.decode("ascii", errors="replace")
            if str(sid).upper().startswith("CL1"):
                return "cl"
            if "McsHdf5ProtocolType" in h5.attrs:
                return "mcs"
    except Exception:
        pass
    return "maxwell"


def detect_manufacturer(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    if ext == ".bxr":
        return "3brain"
    if ext in (".h5", ".hdf5"):
        return _peek_h5_manufacturer(path)
    return "unknown"


def get_loader(manufacturer: str):
    if manufacturer == "3brain":
        return BrainwaveLoader()
    if manufacturer == "maxwell":
        return MaxwellLoader()
    if manufacturer == "cl":
        return CL1Loader()
    if manufacturer == "mcs":
        return MCSLoader()
    if manufacturer == "axion":
        return AxionLoader()
    raise ValueError(f"Unknown manufacturer '{manufacturer}'. "
                     f"Supported: 3brain, maxwell, cl, mcs, axion")


# ── epoch comparison helper ────────────────────────────────────────────────────
def _parse_epoch(s: str):
    try:
        a, b = s.split(",")
        return float(a), float(b)
    except Exception:
        return None


def _epoch_stats(spike_times: dict, start: float, dur: float,
                 all_channels: list) -> dict:
    end = start + dur
    sub = {ch: t[(t >= start) & (t < end)]
           for ch, t in spike_times.items()}
    all_t = np.concatenate(list(sub.values())) if sub else np.array([])

    pop_rate = 0.0
    if all_t.size > 0 and dur > 0:
        pop_rate = float(all_t.size / dur / max(len(all_channels), 1))

    n_active  = sum(1 for t in sub.values() if t.size > 0)
    n_total   = max(len(all_channels), 1)
    return dict(pop_rate=pop_rate,
                active_area=float(n_active / n_total))


# ── main pipeline ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        prog="cascade",
        description=(
            "CASCADE — Criticality Avalanche Spike Cross-platform Analysis "
            "Detection Engine\n"
            "Unified MEA analysis pipeline supporting multiple manufacturers.\n\n"
            "Outputs: raster plot, rate heatmap, burst detection, active area,\n"
            "         connectivity, neuronal avalanches, criticality metrics\n"
            "         (DCC, branching ratio, shape collapse), entropy.\n\n"
            "Other CASCADE tools:\n"
            "  cascade-info      Quick recording inspector (duration, channels,\n"
            "                    active electrodes, time-window support)\n"
            "                    Usage: python scanner.py --help\n\n"
            "  cascade-ltpltd    LTP/LTD classification from firing-rate CSVs\n"
            "                    Usage: python ltpltd.py --help\n\n"
            "Example:\n"
            "  cascade --input recording.bxr --output results/ --ae 5 --ava_size 2 --verbose"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    ap.add_argument("--version", action="version", version="CASCADE v0.9.2")
    ap.add_argument("--test", action="store_true",
                    help="Run the built-in self-test with synthetic data. "
                         "No recording file required.")

    # ── input / output ─────────────────────────────────────────────────────────
    io_group = ap.add_argument_group("Input / Output")
    io_group.add_argument("--input",
                    default=None,
                    help="Path to the data file (.bxr or .h5). "
                         "If omitted, CASCADE searches the current folder.")
    io_group.add_argument("--manufacturer",
                    required=False, default=None,
                    choices=["3brain", "maxwell", "cl", "mcs", "axion"],
                    help="MEA manufacturer (required except with --test). "
                         "3brain (.bxr), maxwell (.h5), cl (.h5), mcs (.h5), axion (.csv)")
    io_group.add_argument("--well",
                    default=None,
                    help="Well or stream name to analyse (e.g. Well_A1 for 3Brain, "
                         "data0000 for Maxwell). Default: first available well.")
    io_group.add_argument("--output",
                    default="mea_out",
                    help="Output folder for all figures and data files. "
                         "Created automatically if it does not exist. (default: mea_out)")
    io_group.add_argument("--svg",
                    action="store_true",
                    help="Also save all figures as SVG.")
    io_group.add_argument("--jpg",
                    action="store_true",
                    help="Also save all figures as JPG (300 dpi).")
    io_group.add_argument("--verbose",
                    action="store_true",
                    help="Print detailed progress messages to the terminal.")

    # ── recording window ───────────────────────────────────────────────────────
    win_group = ap.add_argument_group("Recording Window")
    win_group.add_argument("--start",
                    type=float, default=0.0,
                    help="Start of the analysis window in seconds. (default: 0)")
    win_group.add_argument("--dur",
                    type=float, default=None,
                    help="Duration of the analysis window in seconds. "
                         "Default: full recording length.")

    # ── spike filtering ────────────────────────────────────────────────────────
    spike_group = ap.add_argument_group("Spike Filtering")
    spike_group.add_argument("--ae",
                    type=float, default=0.0,
                    help="Active electrode filter. Keep only channels that fired "
                         "at least X spikes per minute. Set to 0 to keep all "
                         "channels. (default: 0 = off)")
    spike_group.add_argument("--max_rate",
                    type=float, default=0.0,
                    help="Noisy-electrode filter: remove channels whose firing rate "
                         "exceeds this value (spikes per minute). Electrodes with "
                         "abnormally high rates are typically recording electrical "
                         "noise rather than real neurons. Set to 0 to disable. "
                         "Recommended starting value: 1000 (≈ 17 Hz). (default: 0 = off)")
    spike_group.add_argument("--max_rpv",
                    type=float, default=0.0,
                    help="Noisy-electrode filter: remove channels where more than "
                         "this percentage of inter-spike intervals are shorter than "
                         "2 ms (refractory period violation). Real neurons cannot "
                         "fire faster than ~500 Hz; a high violation rate indicates "
                         "electrical noise. Set to 0 to disable. "
                         "Recommended starting value: 1.0 (1%%). (default: 0 = off)")
    spike_group.add_argument("--max_amplitude",
                    type=float, default=0.0,
                    help="Noisy-electrode filter: remove channels where any single "
                         "spike has a peak-to-peak waveform amplitude exceeding this "
                         "value (µV). Electrodes with abnormally large spikes are "
                         "typically picking up mechanical artefacts or electrical "
                         "interference rather than real action potentials. "
                         "Set to 0 to disable. Recommended starting value: 200 (µV). "
                         "Supported: CL1, 3Brain, Maxwell (if stored), MCS, Axion "
                         "— all in µV. (default: 0 = off)")
    spike_group.add_argument("--bin",
                    type=float, default=1.0,
                    help="Bin size in seconds for computing per-channel firing "
                         "rates (used for heatmap and connectivity). (default: 1.0)")

    # ── avalanche detection ────────────────────────────────────────────────────
    ava_group = ap.add_argument_group("Avalanche Detection")
    ava_group.add_argument("--ava_size",
                    type=int, default=2,
                    help="Minimum avalanche size (total spikes). Avalanches with "
                         "fewer than X spikes are discarded. (default: 2)")
    ava_group.add_argument("--ava_dur",
                    type=int, default=1,
                    help="Minimum avalanche duration (bins). Avalanches spanning "
                         "fewer than X bins are discarded. Default is 1 (no "
                         "duration filter) because single-bin events are valid "
                         "avalanches — they represent synchronous bursts where "
                         "many channels fire within one Δt window. (default: 1)")
    ava_group.add_argument("--aval_dt",
                    type=float, default=None,
                    help="Avalanche bin width in seconds. Default: automatically "
                         "inferred from the mean inter-spike interval (ISI) of the "
                         "population, clipped to 1–50 ms.")
    ava_group.add_argument("--perc",
                    type=float, default=0.0,
                    help="Percentile threshold for avalanche detection. The bin "
                         "count must exceed this percentile of all bin counts to "
                         "be included in an avalanche. 0.0 = any spike counts "
                         "(recommended for clean in-vitro spike data). 0.25 = "
                         "Habibollahi 2022 default for noisy recordings. (default: 0.0)")

    # ── criticality metrics ────────────────────────────────────────────────────
    crit_group = ap.add_argument_group("Criticality Metrics")
    crit_group.add_argument("--burst_min",
                    type=int, default=12,
                    help="Upper limit for the xmin search in the avalanche SIZE "
                         "power-law fit. The algorithm tries xmin values from 1 "
                         "to this number and picks the best one. Recommended: "
                         "8–20. (default: 12)")
    crit_group.add_argument("--t_min",
                    type=int, default=4,
                    help="Upper limit for the xmin search in the avalanche "
                         "DURATION power-law fit. Recommended: 3–6. (default: 4)")
    crit_group.add_argument("--legacy_powerlaw",
                    action="store_true",
                    help="Normalise the power-law MLE over only the observed "
                         "avalanche values instead of every integer between the "
                         "cut-offs. This reproduces CASCADE <= 0.9.0 and "
                         "Habibollahi 2022 tplfit.py exactly, but the exponent is "
                         "biased (severely below ~1000 avalanches). Use ONLY to "
                         "reproduce previously published numbers, never for new "
                         "analyses. (default: off = correct normalisation)")
    crit_group.add_argument("--pl_xmin_S",
                    type=int, default=None,
                    help="Fix the lower cut-off (xmin) of the size power-law fit "
                         "instead of searching for it. tau, beta_pred, DCC, SCE "
                         "and the log-log plot all use it; the automatic fit is "
                         "still written to the report for comparison. "
                         "Default: auto-detected (recommended).")
    crit_group.add_argument("--pl_xmin_D",
                    type=int, default=None,
                    help="Fix the lower cut-off (xmin) of the duration power-law "
                         "fit instead of searching for it. Default: auto-detected.")
    crit_group.add_argument("--skip_rare_sizes",
                    action="store_true",
                    help="On the size log-log plot, exclude avalanche sizes that "
                         "occurred only once in the recording. These rare points "
                         "can distort the visual fit line. Does NOT affect the "
                         "MLE exponent value.")
    crit_group.add_argument("--skip_rare_durations",
                    action="store_true",
                    help="Same as --skip_rare_sizes but for the duration log-log "
                         "plot.")

    # ── shape collapse ─────────────────────────────────────────────────────────
    sc_group = ap.add_argument_group("Shape Collapse")
    sc_group.add_argument("--sc_Dmin",
                    type=int, default=4,
                    help="Minimum avalanche duration (in bins) to include in the "
                         "shape collapse analysis. (default: 4)")
    sc_group.add_argument("--sc_Dmax",
                    type=int, default=20,
                    help="Maximum avalanche duration (in bins) to include in the "
                         "shape collapse analysis. (default: 20)")

    # ── burst detection ────────────────────────────────────────────────────────
    burst_group = ap.add_argument_group("Burst Detection")
    burst_group.add_argument("--burst_min_spikes",
                    type=int, default=5,
                    help="Minimum number of spikes required for a run to be "
                         "counted as a network burst. Lower values detect more, "
                         "shorter bursts (MEAanalysis default: 5). "
                         "Higher values require denser events. (default: 5)")
    burst_group.add_argument("--burst_max_gap",
                    type=float, default=0.1,
                    help="Maximum allowed gap between consecutive spikes (seconds) "
                         "within a burst. Spikes further apart than this split the "
                         "burst into two. (default: 0.1)")

    # ── connectivity ───────────────────────────────────────────────────────────
    conn_group = ap.add_argument_group("Connectivity")
    conn_group.add_argument("--conn_thr",
                    type=float, default=0.3,
                    help="Pearson correlation threshold for drawing a functional "
                         "connectivity edge between two channels. Pairs with "
                         "|r| >= this value are connected. (default: 0.3)")

    # ── synchrony ──────────────────────────────────────────────────────────────
    sync_group = ap.add_argument_group("Synchrony")
    sync_group.add_argument("--sync_metrics",
                    choices=["none", "core", "xcorr", "kreuz", "index", "all"],
                    default="core",
                    help="Which synchrony metrics to compute. 'core' reports the two "
                         "that resolve on every manufacturer and discriminate best "
                         "between recordings: area under the normalised "
                         "cross-correlation, and the Synchrony Index. 'xcorr' adds "
                         "the raw area and both full-width-at-half-height values; "
                         "'kreuz' = 1 - Kreuz SPIKE-distance; 'index' = Synchrony "
                         "Index alone; 'all' = everything, mirroring the Axion "
                         "Neural Metric Tool. Note 'kreuz' and 'index' are O(N^2) in "
                         "channels and much slower than the correlogram. "
                         "(default: core)")
    sync_group.add_argument("--sync_window",
                    type=float, default=20.0,
                    help="Synchrony window in MILLISECONDS: the window around zero "
                         "lag used to compute the area under the cross-correlation. "
                         "A short window (5) emphasises synchrony on a millisecond "
                         "timescale; a long window (100) captures slow timescales. "
                         "(default: 20.0)")
    sync_group.add_argument("--sync_bin",
                    type=float, default=1.0,
                    help="Cross-correlogram bin width in MILLISECONDS. This is an "
                         "absolute duration and is deliberately NOT inferred from "
                         "the sampling rate, so that synchrony values are "
                         "comparable across manufacturers. (default: 1.0)")
    sync_group.add_argument("--sync_max_lag",
                    type=float, default=500.0,
                    help="Maximum lag of the cross-correlogram in MILLISECONDS. "
                         "Must be wide enough to contain the correlogram peak, or "
                         "the full width at half height cannot be measured. "
                         "(default: 500.0)")
    sync_group.add_argument("--sync_sigma",
                    type=float, default=None,
                    help="Gaussian timescale sigma in MILLISECONDS for the "
                         "Synchrony Index. Small relative to the mean ISI measures "
                         "spike-timing coincidence; large approaches a firing-rate "
                         "comparison. (default: same as --sync_window)")
    sync_group.add_argument("--sync_max_channels",
                    type=int, default=64,
                    help="Cap on the number of channels used for the O(N^2) "
                         "measures (--sync_metrics kreuz/index/all). Channels are "
                         "subsampled with a fixed seed, and the output records that "
                         "subsampling occurred. 0 = no cap. (default: 64)")

    # ── entropy and active area ────────────────────────────────────────────────
    ea_group = ap.add_argument_group("Entropy and Active Area")
    ea_group.add_argument("--entropy_win",
                    type=float, default=10.0,
                    help="Window size in seconds for computing spatial entropy "
                         "and active electrode fraction over time. (default: 10.0)")

    # ── epoch comparison ───────────────────────────────────────────────────────
    ep_group = ap.add_argument_group("Epoch Comparison (optional)")
    ep_group.add_argument("--epochA",
                    default=None,
                    help="First epoch window as 'start,duration' in seconds. "
                         "Example: --epochA 0,60 for the first 60 seconds.")
    ep_group.add_argument("--epochB",
                    default=None,
                    help="Second epoch window as 'start,duration' in seconds. "
                         "Example: --epochB 120,60 for 60 seconds starting at 2 min. "
                         "Both --epochA and --epochB must be set to compare.")

    args = ap.parse_args()

    # ── self-test ──────────────────────────────────────────────────────────────
    if args.test:
        from cascade.tests.self_test import run_self_test
        sys.exit(run_self_test())

    if args.manufacturer is None:
        ap.error("--manufacturer is required. Choose from: 3brain, maxwell, cl, mcs, axion")

    svg = args.svg
    jpg = args.jpg
    v   = args.verbose

    banner("CASCADE  •  Criticality Avalanche Spike Cross-platform Analysis Detection Engine  •  v0.9.2", enabled=v)

    # ── find file ──────────────────────────────────────────────────────────────
    path = args.input
    if not path:
        for pattern in ("*.bxr", "*.raw.h5", "*.h5"):
            hits = sorted(glob.glob(pattern))
            if hits:
                path = hits[0]
                break
    if not path:
        say("No data file found. Use --input <path>", "[FAIL]", RED)
        sys.exit(1)

    # ── load manufacturer ──────────────────────────────────────────────────────
    mfr    = args.manufacturer
    loader = get_loader(mfr)
    say(f"File: {path}", "[DIR]", CYAN, enabled=v)
    say(f"Manufacturer: {loader.manufacturer_name()}", "[CHECK]", CYAN, enabled=v)

    # ── load data ──────────────────────────────────────────────────────────────
    step("Loading data", enabled=v)
    rec = loader.load(path,
                      well       = args.well,
                      start_s    = args.start,
                      duration_s = args.dur)

    say(f"SR={rec.sampling_rate:.1f} Hz  |  "
        f"channels={rec.n_channels}  |  "
        f"duration={rec.duration:.1f} s  |  "
        f"spikes={rec.n_spikes():,}", "[OK]", GRN, enabled=v)

    os.makedirs(args.output, exist_ok=True)
    well_label = rec.metadata.get("well", "well")

    # ── active-electrode filter ────────────────────────────────────────────────
    step("Active-electrode filter", enabled=v)
    spike_times = filter_active_electrodes(rec.spike_times, args.ae, rec.duration)
    n_active    = len(spike_times)
    say(f"Channels kept (>= {args.ae} spikes/min): {n_active}", "[DATA]", CYAN, enabled=v)
    if n_active == 0:
        print("[WARN] No active channels found after filtering. "
              "The recording may have very low firing activity, or --ae is set too high. "
              "Try lowering --ae or removing it entirely.")
    elif n_active == 1:
        print("[WARN] Only 1 active channel found. Connectivity analysis will be "
              "skipped (requires at least 2 channels). Avalanche and criticality metrics "
              "may be unreliable with a single electrode. Consider lowering --ae.")

    # ── noisy-electrode filter ─────────────────────────────────────────────────
    if args.max_rate > 0 or args.max_rpv > 0:
        step("Noisy-electrode filter", enabled=v)
        spike_times, n_removed, removed_ch = filter_noisy_electrodes(
            spike_times, rec.duration,
            max_rate_per_min=args.max_rate,
            max_rpv_pct=args.max_rpv,
        )
        say(f"Noisy channels removed: {n_removed} "
            f"(max_rate={args.max_rate} spk/min, max_rpv={args.max_rpv}%)",
            "[DATA]", CYAN, enabled=v)
        if n_removed > 0 and v:
            say(f"Removed channel IDs: {removed_ch[:20]}"
                f"{'…' if len(removed_ch) > 20 else ''}",
                "[DATA]", CYAN, enabled=v)
        n_active = len(spike_times)

    # ── amplitude filter ───────────────────────────────────────────────────────
    if args.max_amplitude > 0:
        step("Amplitude filter", enabled=v)
        if rec.spike_amplitudes is None:
            say("--max_amplitude ignored: no amplitude data found in this file "
                "(Maxwell requires the optional amplitude field; "
                "3Brain requires SpikeForms).", "[WARN]", YEL, enabled=True)
        else:
            spike_times, n_removed, removed_ch = filter_amplitude_electrodes(
                spike_times, rec.spike_amplitudes,
                max_amplitude_uv=args.max_amplitude,
            )
            say(f"High-amplitude channels removed: {n_removed} "
                f"(max_amplitude={args.max_amplitude} µV)",
                "[DATA]", CYAN, enabled=v)
            if n_removed > 0 and v:
                say(f"Removed channel IDs: {removed_ch[:20]}"
                    f"{'…' if len(removed_ch) > 20 else ''}",
                    "[DATA]", CYAN, enabled=v)
            n_active = len(spike_times)

    all_spikes = rec.all_spike_times() if (args.ae <= 0 and args.max_rate <= 0
                                           and args.max_rpv <= 0 and args.max_amplitude <= 0) else \
        np.sort(np.concatenate(list(spike_times.values()))) \
        if spike_times else np.array([])

    # ── raster ────────────────────────────────────────────────────────────────
    step("Raster & rate heatmap", enabled=v)
    plot_raster(spike_times,
                f"{args.output}/raster_max_{well_label}.png",
                f"Spike raster — {well_label}", svg, jpg,
                n_channels=rec.n_channels)

    # active: remap electrode IDs to 0..N-1 so Y-axis max = number of active electrodes
    _sorted_chs = sorted(spike_times.keys())
    _ch_seq_map = {ch: i for i, ch in enumerate(_sorted_chs)}
    _st_active  = {_ch_seq_map[ch]: times for ch, times in spike_times.items()}
    plot_raster(_st_active,
                f"{args.output}/raster_active_{well_label}.png",
                f"Spike raster (active) — {well_label}", svg, jpg)

    # ── firing rates ──────────────────────────────────────────────────────────
    rate_df = bin_spikes(spike_times, rec.duration, args.bin)
    write_tsv(rate_df, f"{args.output}/rates_{well_label}.tsv")
    plot_rate_heatmap(rate_df,
                      f"{args.output}/rate_heatmap_max_{well_label}.png",
                      f"Rate heatmap — {well_label} (bin={args.bin}s)", svg, jpg,
                      n_channels=rec.n_channels)
    _rate_df_active = rate_df.copy()
    _rate_df_active["ch_idx"] = _rate_df_active["ch_idx"].map(_ch_seq_map)
    plot_rate_heatmap(_rate_df_active,
                      f"{args.output}/rate_heatmap_active_{well_label}.png",
                      f"Rate heatmap (active) — {well_label} (bin={args.bin}s)", svg, jpg)
    plot_population_rate(rate_df,
                         f"{args.output}/population_rate_{well_label}.png",
                         well_label, svg, jpg)

    # ── ISI metrics ────────────────────────────────────────────────────────────
    step("ISI metrics", enabled=v)
    isi_df = compute_isi_metrics(spike_times, rec.sampling_rate)
    write_tsv(isi_df, f"{args.output}/isi_metrics_{well_label}.tsv")

    # ── network bursts ─────────────────────────────────────────────────────────
    step("Network burst detection", enabled=v)
    burst_df, bsum = compute_burst_metrics(all_spikes,
                                           min_spikes=args.burst_min_spikes,
                                           max_gap_s=args.burst_max_gap)
    write_tsv(burst_df, f"{args.output}/bursts_{well_label}.tsv")
    plot_burst_timeline(burst_df, rec.duration,
                        f"{args.output}/bursts_timeline_{well_label}.png",
                        well_label, svg, jpg)
    say(f"Bursts: n={bsum['n_bursts']}  "
        f"mean_dur={bsum['mean_duration_s']:.2f}s", "[DATA]", CYAN, enabled=v)
    if bsum["n_bursts"] == 0:
        print("[WARN] No network bursts detected. This is normal for low-activity "
              "recordings or recordings with sparse, uncorrelated firing.")

    # ── active area ────────────────────────────────────────────────────────────
    act_df = compute_active_area(spike_times, rec.duration, args.entropy_win)
    write_tsv(act_df, f"{args.output}/active_area_{well_label}.tsv")
    plot_active_area_curve(act_df,
                           f"{args.output}/active_area_curve_{well_label}.png",
                           well_label, svg, jpg)

    # ── waveforms (3Brain only) ────────────────────────────────────────────────
    if mfr == "3brain" and rec.metadata.get("has_waveforms"):
        step("Spike waveforms", enabled=v)
        waves = loader.load_waveforms(path,
                                      well       = rec.metadata.get("well"),
                                      start_s    = args.start,
                                      duration_s = rec.duration)
        if waves:
            plot_waveforms(waves,
                           f"{args.output}/wave_templates_{well_label}.png",
                           well_label, svg, jpg)

    # ── connectivity ───────────────────────────────────────────────────────────
    step("Connectivity (Pearson correlation)", enabled=v)
    if rate_df.empty or rate_df["ch_idx"].nunique() < 2:
        say("Skipped -- fewer than 2 active channels", "[WARN]", YEL, enabled=v)
        C, ch_ids = None, []
    else:
        C, ch_ids = compute_connectivity(rate_df, zscore=True)
    if C is not None:
        np.savetxt(f"{args.output}/corr_{well_label}.tsv", C, delimiter="\t")
        plot_corr_heatmap(C,
                          f"{args.output}/corr_heatmap_max_{well_label}.png",
                          args.conn_thr, well_label, svg, jpg,
                          channel_ids=ch_ids,
                          n_channels=rec.n_channels)
        plot_corr_heatmap(C,
                          f"{args.output}/corr_heatmap_active_{well_label}.png",
                          args.conn_thr, well_label, svg, jpg,
                          channel_ids=ch_ids)
        deg_df, csum = summarize_connectivity(C, ch_ids, args.conn_thr)
        write_tsv(deg_df, f"{args.output}/degrees_{well_label}.tsv")
        plot_degree_histogram(deg_df,
                              f"{args.output}/degree_hist_{well_label}.png",
                              well_label, svg, jpg)
        say(f"Edges |r|>={args.conn_thr}: {csum['n_edges']}  "
            f"density={csum['edge_density']:.4f}", "[DATA]", CYAN, enabled=v)

    # ── synchrony ──────────────────────────────────────────────────────────────
    # Operates directly on spike times, so it is identical for every manufacturer.
    # All timing parameters are absolute (ms), never sampling-rate derived, so the
    # values are comparable across devices — see core/synchrony.py.
    if args.sync_metrics != "none":
        step(f"Synchrony ({args.sync_metrics})", enabled=v)
        sync_cap = args.sync_max_channels if args.sync_max_channels > 0 else None
        sync_m, sync_cc = compute_synchrony(
            spike_times,
            duration_s=rec.duration,
            window_s=args.sync_window / 1000.0,
            bin_s=args.sync_bin / 1000.0,
            max_lag_s=args.sync_max_lag / 1000.0,
            sigma_s=(args.sync_sigma / 1000.0) if args.sync_sigma else None,
            metrics=args.sync_metrics,
            max_channels=sync_cap,
        )
        write_tsv(synchrony_to_dataframe(sync_m),
                  f"{args.output}/synchrony_{well_label}.tsv")
        if sync_cc is not None:
            plot_cross_correlogram(sync_cc,
                                   f"{args.output}/crosscorrelogram_norm_{well_label}.png",
                                   well_label, window_s=args.sync_window / 1000.0,
                                   normalized=True, save_svg=svg, save_jpg=jpg)
            plot_cross_correlogram(sync_cc,
                                   f"{args.output}/crosscorrelogram_{well_label}.png",
                                   well_label, window_s=args.sync_window / 1000.0,
                                   normalized=False, save_svg=svg, save_jpg=jpg)
            _msg = (f"Pooled correlogram: {sync_cc['n_pairs']} pairs  "
                    f"area_norm={sync_m.get('area_xcorr_norm', float('nan')):.6g}")
            if "fwhh_xcorr_norm" in sync_m:
                _msg += f"  FWHH={sync_m['fwhh_xcorr_norm']*1000:.2f} ms"
            say(_msg, "[DATA]", CYAN, enabled=v)
            if sync_m.get("fwhh_resolution_limited"):
                say(f"FWHH is not resolved at these settings: "
                    f"{sync_m.get('fwhh_limit_reason','')}. The area metrics are "
                    f"unaffected.", "[WARN]", YEL, enabled=True)
        if "kreuz_metric" in sync_m:
            say(f"Kreuz metric (1-SPIKE-distance): {sync_m['kreuz_metric']:.4f}"
                f"{'  [channels subsampled]' if sync_m.get('kreuz_subsampled') else ''}",
                "[DATA]", CYAN, enabled=v)
        if "synchrony_index" in sync_m:
            say(f"Synchrony Index: {sync_m['synchrony_index']:.4f}"
                f"{'  [channels subsampled]' if sync_m.get('index_subsampled') else ''}",
                "[DATA]", CYAN, enabled=v)

    # ── avalanches ─────────────────────────────────────────────────────────────
    step("Neuronal avalanches", enabled=v)
    dt = args.aval_dt if args.aval_dt is not None else infer_dt(all_spikes)
    say(f"Δt = {dt*1000:.2f} ms  |  "
        f"perc threshold = {args.perc*100:.0f}%", "[FAST]", CYAN, enabled=v)

    av_df = compute_avalanches(all_spikes, dt, perc=args.perc)

    # apply size/duration threshold
    av_df = av_df[
        (av_df["avalanche_size"] >= args.ava_size) &
        (av_df["duration_bins"]  >= args.ava_dur)
    ].reset_index(drop=True)

    write_tsv(av_df, f"{args.output}/avalanches_{well_label}.tsv")
    plot_avalanche_distributions(av_df,
                                 f"{args.output}/avalanches_{well_label}",
                                 well_label, svg, jpg)
    say(f"Avalanches after threshold (size>={args.ava_size}, dur>={args.ava_dur}): {len(av_df)}", "[DATA]", CYAN, enabled=v)
    if av_df.empty:
        print("[WARN] No avalanches passed the size/duration threshold. "
              "Criticality metrics (tau, alpha, DCC, branching ratio) cannot be computed. "
              "Try lowering --ava_size (current: {}) or --ava_dur (current: {}) or --perc (current: {}).".format(
                  args.ava_size, args.ava_dur, args.perc))
    elif len(av_df) < 50:
        print(f"[WARN] Only {len(av_df)} avalanches detected. Power-law fitting "
              "requires many events for reliable results -- interpret criticality metrics "
              "with caution. Consider lowering --ava_size, --ava_dur, or --perc.")

    # ── criticality metrics ────────────────────────────────────────────────────
    step("Criticality metrics (exponents, DCC, BR, shape collapse)", enabled=v)
    sizes     = av_df["avalanche_size"].to_numpy() if not av_df.empty else np.array([])
    durations = av_df["duration_bins"].to_numpy()  if not av_df.empty else np.array([])

    if args.legacy_powerlaw:
        say("--legacy_powerlaw: using the biased pre-0.9.1 normalisation. "
            "Valid only for reproducing previously published numbers.",
            "[WARN]", YEL, enabled=True)
    exp_metrics = compute_exponents(sizes, durations,
                                    burst_min=args.burst_min,
                                    t_min=args.t_min,
                                    legacy_normalization=args.legacy_powerlaw)
    auto_metrics = None
    if args.pl_xmin_S is not None or args.pl_xmin_D is not None:
        # Keep the automatic fit so a hand-picked xmin is always reported next to it.
        auto_metrics = exp_metrics
        exp_metrics = compute_exponents(sizes, durations,
                                        burst_min=args.burst_min,
                                        t_min=args.t_min,
                                        legacy_normalization=args.legacy_powerlaw,
                                        xmin_s=args.pl_xmin_S,
                                        xmin_d=args.pl_xmin_D)
        say(f"xmin set manually (--pl_xmin_S={args.pl_xmin_S}, --pl_xmin_D={args.pl_xmin_D}): "
            f"tau/alpha/DCC/SCE use it. Automatic fit: tau={auto_metrics['tau']:.4f} "
            f"(xmin_S={auto_metrics['xmin_s']}), alpha={auto_metrics['alpha']:.4f} "
            f"(xmin_D={auto_metrics['xmin_d']}), DCC={auto_metrics['dcc']:.4f}",
            "[WARN]", YEL, enabled=True)

    # log–log plots with fitted lines
    # --skip_rare_sizes/durations excludes single-occurrence points from the OLS line
    xmin_s = exp_metrics.get("xmin_s")
    xmin_d = exp_metrics.get("xmin_d")

    plot_loglog_linearization(
        sizes,
        f"{args.output}/loglog_size_{well_label}.png",
        f"Log–log size — {well_label}",
        tau=exp_metrics.get("tau"), xmin=xmin_s,
        exclude_singletons=args.skip_rare_sizes,
        save_svg=svg, save_jpg=jpg)
    plot_loglog_linearization(
        durations,
        f"{args.output}/loglog_duration_{well_label}.png",
        f"Log–log duration — {well_label}",
        tau=exp_metrics.get("alpha"), xmin=xmin_d,
        exclude_singletons=args.skip_rare_durations,
        save_svg=svg, save_jpg=jpg)

    # branching ratio
    counts, _ = counts_per_dt(all_spikes, dt)
    br, br_h, br_pts = compute_branching_ratio(counts)
    write_tsv(br_pts, f"{args.output}/branching_points_{well_label}.tsv")

    # shape collapse
    Dmin_eff = max(args.sc_Dmin, args.ava_dur)
    sc_beta, sce, sc_minerr, sc_details = shape_collapse(
        counts, exp_metrics.get("beta_pred", np.nan),
        Dmin=Dmin_eff, Dmax=args.sc_Dmax, min_size=args.ava_size)
    plot_shape_collapse(sc_details,
                        f"{args.output}/shape_collapse_{well_label}.png",
                        well_label, svg, jpg)

    # write text report
    crit_metrics = dict(
        well         = well_label,
        dt           = dt,
        ava_size_thr = args.ava_size,
        ava_dur_thr  = args.ava_dur,
        perc         = args.perc,
        tau          = exp_metrics.get("tau", np.nan),
        xmin_s       = exp_metrics.get("xmin_s", np.nan),
        alpha        = exp_metrics.get("alpha", np.nan),
        xmin_d       = exp_metrics.get("xmin_d", np.nan),
        beta_pred    = exp_metrics.get("beta_pred", np.nan),
        beta_fit     = exp_metrics.get("beta_fit", np.nan),
        dcc          = exp_metrics.get("dcc", np.nan),
        br           = br,
        br_intercept = br_h,
        sc_beta      = sc_beta,
        sce          = sce,
        sc_minerr    = sc_minerr,
        auto_fit     = auto_metrics,
    )
    write_criticality_report(crit_metrics,
                             f"{args.output}/criticality_report_{well_label}.txt")
    say(f"DCC={exp_metrics.get('dcc', float('nan')):.4f}  "
        f"BR={br:.4f}  SCE={sce:.4f}", "[CRIT]", CYAN, enabled=v)

    # ── entropy ────────────────────────────────────────────────────────────────
    step("Spatial entropy", enabled=v)
    ent_df = compute_entropy(spike_times, rec.duration, args.entropy_win)
    write_tsv(ent_df, f"{args.output}/entropy_{well_label}.tsv")
    plot_entropy_curve(ent_df,
                       f"{args.output}/entropy_curve_{well_label}.png",
                       well_label, svg, jpg)

    # ── epoch comparison (optional) ────────────────────────────────────────────
    if args.epochA and args.epochB:
        step("Epoch comparison", enabled=v)
        epA = _parse_epoch(args.epochA)
        epB = _parse_epoch(args.epochB)
        if epA and epB:
            stats = {
                "A": _epoch_stats(spike_times, epA[0], epA[1], rec.channel_ids),
                "B": _epoch_stats(spike_times, epB[0], epB[1], rec.channel_ids),
            }
            plot_epoch_compare(stats,
                               f"{args.output}/epoch_compare_{well_label}.png",
                               well_label, svg, jpg)
            say(f"Epoch A PopRate={stats['A']['pop_rate']:.3f} Hz  "
                f"->  B PopRate={stats['B']['pop_rate']:.3f} Hz", "[DATA]", CYAN, enabled=v)
        else:
            say("Bad epoch format. Use 'start,dur' e.g. --epochA 0,60", "[WARN]", YEL)

    banner(f"ALL DONE  •  outputs in {args.output}", enabled=v)
    print(f"[OK] Done. Outputs saved in: {args.output}")


def _run_with_log():
    import io as _io

    _cmd = "cascade " + " ".join(sys.argv[1:])
    _buf = _io.StringIO()
    _tee = _Tee(sys.stdout, _buf)
    sys.stdout = _tee

    _out = "mea_out"
    for _i, _arg in enumerate(sys.argv):
        if _arg == "--output" and _i + 1 < len(sys.argv):
            _out = sys.argv[_i + 1]
            break

    _status = "SUCCESS"
    _error  = None
    try:
        main()
    except SystemExit as exc:
        # sys.exit(0) / sys.exit() (--test pass, --help) is a clean exit, not a failure
        if exc.code not in (0, None):
            _status = "FAILED"
        raise
    except Exception as exc:
        _status = "FAILED"
        _error  = str(exc)
        raise
    finally:
        sys.stdout = _tee._original
        write_run_log(_out, _cmd, _status, error=_error, output=_buf.getvalue())


if __name__ == "__main__":
    _run_with_log()
