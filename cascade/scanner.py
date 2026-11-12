#!/usr/bin/env python
"""
scanner.py
----------
CASCADE quick recording inspector.

Reports duration, total channels, and active electrodes for a MEA file.
Supports all five manufacturers (3Brain .bxr; Maxwell, Cortical Labs and MCS .h5;
Axion .csv) and plain spike CSV files.

The key feature is --start and --dur: you can inspect any time window inside
a recording. This is useful when a single file contains baseline + treatment:

  # Full recording
  cascade-info --input recording.bxr --manufacturer 3brain

  # Baseline window (0 – 60 s)
  cascade-info --input recording.bxr --manufacturer 3brain --start 0 --dur 60

  # Post-treatment window (61 – 120 s)
  cascade-info --input recording.bxr --manufacturer 3brain --start 61 --dur 59

  # Scan every recording of one manufacturer in a folder
  cascade-info --input /data/recordings/ --manufacturer cl --ae 5

Output is printed to the terminal and, with --output DIR, also saved as
DIR/info_summary.tsv (one row per file).
"""

import argparse
import math
import os
import sys
import numpy as np

# ── optional pandas (required for CSV input mode) ──────────────────────────────
try:
    import pandas as pd
    _PANDAS = True
except ImportError:
    _PANDAS = False


# ── CASCADE loaders ─────────────────────────────────────────────────────────────
# Resolve the package root so this script can be run directly or as a module.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from cascade import __version__
from cascade.utils.citation import add_citation_argument
from cascade.loaders.brainwave_loader import BrainwaveLoader
from cascade.loaders.maxwell_loader   import MaxwellLoader
from cascade.loaders.cl_loader        import CLLoader
from cascade.loaders.mcs_loader       import MCSLoader
from cascade.loaders.axion_loader     import AxionLoader
from cascade.loaders                  import read_error_message, EXTENSIONS


# ── helpers ─────────────────────────────────────────────────────────────────────
def _required_spikes(ae_per_min: float, duration_s: float) -> int:
    """Convert 'X spikes/min' into an absolute count for the window length."""
    if ae_per_min <= 0:
        return 1
    minutes = max(0.0, duration_s) / 60.0
    return max(1, int(math.ceil(ae_per_min * minutes)))


def _fmt_duration(s: float) -> str:
    """Return a human-readable duration string, e.g. '312.4 s  (5 min 12 s)'."""
    m, sec = divmod(int(s), 60)
    if m == 0:
        return f"{s:.2f} s"
    return f"{s:.2f} s  ({m} min {sec} s)"


# ── CSV mode ────────────────────────────────────────────────────────────────────
_TIME_ALIASES  = ["time_s", "t", "time", "timestamp", "sec", "seconds"]
_FRAME_ALIASES = ["frame", "frame_no", "frame_nos", "sample", "samples"]
_CH_ALIASES    = ["ch_idx", "ch", "channel", "channel_idx",
                  "electrode", "electrode_idx"]


def _pick_col(df, aliases):
    lower = {c.strip().lower(): c for c in df.columns}
    for a in aliases:
        if a in lower:
            return lower[a]
    return None


def _scan_csv(path: str, ae_per_min: float, sr: float,
              start_s: float, dur_s: float) -> dict:
    if not _PANDAS:
        raise RuntimeError(
            "pandas is required for CSV mode.  Install: pip install pandas"
        )
    df = pd.read_csv(path)
    if df.empty:
        return dict(manufacturer="CSV", well="-", sampling_rate=None,
                    full_duration=0.0, window_start=start_s, window_dur=0.0,
                    total_ch=0, active_ch=0, ae_threshold=ae_per_min)

    ch_col    = _pick_col(df, _CH_ALIASES)
    time_col  = _pick_col(df, _TIME_ALIASES)
    frame_col = _pick_col(df, _FRAME_ALIASES)

    if ch_col is None:
        raise ValueError("CSV must have a channel column (ch_idx, ch, channel, …).")
    if time_col is None and frame_col is None:
        raise ValueError("CSV must have a time_s column or a frame column (needs --sr).")

    ch = pd.to_numeric(df[ch_col], errors="coerce")
    ok = ch.notna()
    df, ch = df.loc[ok], ch[ok].astype(int).to_numpy()

    if time_col is not None:
        t = pd.to_numeric(df[time_col], errors="coerce").to_numpy()
    else:
        if sr is None or sr <= 0:
            raise ValueError("CSV uses frame column — provide --sr (sampling rate Hz).")
        t = pd.to_numeric(df[frame_col], errors="coerce").to_numpy() / float(sr)

    valid = np.isfinite(t)
    t, ch = t[valid], ch[valid]

    full_dur = float(t.max()) if t.size else 0.0
    total_ch = int(np.unique(ch).size)

    # apply time window
    end_s = (start_s + dur_s) if dur_s is not None else full_dur
    mask  = (t >= start_s) & (t < end_s)
    t_win, ch_win = t[mask], ch[mask]

    win_dur = float(t_win.max() - t_win.min()) if t_win.size > 1 else 0.0
    req     = _required_spikes(ae_per_min, win_dur)
    if ch_win.size:
        _, cnt = np.unique(ch_win, return_counts=True)
        active = int(np.sum(cnt >= req))
    else:
        active = 0

    return dict(manufacturer="CSV", well=os.path.basename(path),
                sampling_rate=sr,
                full_duration=full_dur,
                window_start=start_s, window_dur=win_dur,
                total_ch=total_ch, active_ch=active,
                ae_threshold=ae_per_min)


# ── loader factory ──────────────────────────────────────────────────────────────
def _get_loader(manufacturer: str):
    if manufacturer == "3brain":
        return BrainwaveLoader()
    if manufacturer == "maxwell":
        return MaxwellLoader()
    if manufacturer == "cl":
        return CLLoader()
    if manufacturer == "mcs":
        return MCSLoader()
    if manufacturer == "axion":
        return AxionLoader()
    raise ValueError(f"Unknown manufacturer '{manufacturer}'. "
                     f"Supported: 3brain, maxwell, cl, mcs, axion")


# ── HDF5 mode (via CASCADE loaders) ────────────────────────────────────────────
def _scan_hdf5(path: str, manufacturer: str, well: str, ae_per_min: float,
               start_s: float, dur_s: float) -> dict:
    """
    Use the appropriate CASCADE loader to get a MEARecording, then compute
    stats for the requested time window.
    """
    import h5py

    loader = _get_loader(manufacturer)

    if manufacturer == "3brain":
        with h5py.File(path, "r") as h5:
            all_wells = [k for k in h5.keys() if k.startswith("Well_")]
        chosen_well = well or (all_wells[0] if all_wells else None)
        manufacturer_label = "3Brain BrainWave"
    elif manufacturer == "cl":
        all_wells = ["default"]
        chosen_well = "default"
        manufacturer_label = "Cortical Labs"
    elif manufacturer == "mcs":
        with h5py.File(path, "r") as h5:
            all_wells = MCSLoader.list_wells(h5)
        chosen_well = well or (all_wells[0] if all_wells else None)
        manufacturer_label = "Multi Channel Systems"
    else:  # maxwell
        with h5py.File(path, "r") as h5:
            all_wells = sorted(
                k for k in h5.get("data_store", {}).keys()
                if k.startswith("data")
            )
        chosen_well = well or ("data0000" if "data0000" in all_wells
                               else (all_wells[0] if all_wells else None))
        manufacturer_label = "Maxwell Biosystems"

    if chosen_well is None:
        raise ValueError(f"No wells/streams found in {path}")

    # load full recording first to get total channels and true duration
    rec_full = loader.load(path, well=chosen_well,
                           start_s=0.0, duration_s=None)
    full_dur  = rec_full.duration
    total_ch  = rec_full.n_channels

    # load the requested window
    rec_win  = loader.load(path, well=chosen_well,
                           start_s=start_s, duration_s=dur_s)
    win_dur  = rec_win.duration

    # count active channels in the window
    req = _required_spikes(ae_per_min, win_dur)
    active = sum(
        1 for times in rec_win.spike_times.values()
        if len(times) >= req
    )

    return dict(
        manufacturer=manufacturer_label,
        well=chosen_well,
        all_wells=all_wells,
        sampling_rate=rec_full.sampling_rate,
        full_duration=full_dur,
        window_start=start_s,
        window_dur=win_dur,
        total_ch=total_ch,
        active_ch=active,
        ae_threshold=ae_per_min,
    )


# ── pretty print ────────────────────────────────────────────────────────────────
def _print_result(path: str, info: dict):
    fname  = os.path.basename(path)
    win_s  = info["window_start"]
    win_d  = info["window_dur"]
    ae     = info["ae_threshold"]

    print()
    print(f"  File        : {fname}")
    print(f"  Manufacturer: {info['manufacturer']}")
    print(f"  Well/stream : {info['well']}")
    if info.get("all_wells"):
        print(f"  All wells   : {', '.join(info['all_wells'])}")
    if info.get("sampling_rate"):
        print(f"  Sampling    : {info['sampling_rate']:.0f} Hz")
    print(f"  Full dur.   : {_fmt_duration(info['full_duration'])}")

    if win_s > 0 or (info.get("window_dur") != info.get("full_duration")):
        print(f"  Window      : {win_s:.1f} s → {win_s + win_d:.1f} s  "
              f"({_fmt_duration(win_d)})")

    print(f"  Total ch    : {info['total_ch']}")
    print(f"  Active ch   : {info['active_ch']}  (>= {ae:.0f} spikes/min "
          f"in window)")
    print()


# ── Axion scan (via AxionLoader) ─────────────────────────────────────────────────
def _scan_axion(path: str, well: str, ae_per_min: float,
                start_s: float, dur_s: float) -> dict:
    loader = AxionLoader()

    rec_full = loader.load(path, well=well, start_s=0.0, duration_s=None)
    full_dur  = rec_full.duration
    total_ch  = rec_full.n_channels

    rec_win  = loader.load(path, well=well, start_s=start_s, duration_s=dur_s)
    win_dur  = rec_win.duration

    req    = _required_spikes(ae_per_min, win_dur)
    active = sum(1 for t in rec_win.spike_times.values() if len(t) >= req)

    # infer available wells from electrode labels
    all_wells = sorted({ch.split("_")[0]
                        for ch in rec_full.channel_ids
                        if "_" in str(ch)})

    return dict(
        manufacturer    = loader.manufacturer_name(),
        well            = well or "all",
        all_wells       = all_wells,
        sampling_rate   = rec_full.sampling_rate,
        full_duration   = full_dur,
        window_start    = start_s,
        window_dur      = win_dur,
        total_ch        = total_ch,
        active_ch       = active,
        ae_threshold    = ae_per_min,
    )


# ── single-file scan ────────────────────────────────────────────────────────────
def scan_file(path: str, manufacturer: str, well: str, ae_per_min: float,
              start_s: float, dur_s: float, sr: float) -> dict:
    if manufacturer == "axion":
        return _scan_axion(path, well, ae_per_min, start_s, dur_s)
    if path.lower().endswith(".csv"):
        return _scan_csv(path, ae_per_min, sr, start_s, dur_s)
    return _scan_hdf5(path, manufacturer, well, ae_per_min, start_s, dur_s)


# ── summary table ───────────────────────────────────────────────────────────────
SUMMARY_NAME = "info_summary.tsv"


def _summary_row(path: str, info: dict, start_s: float, ae_per_min: float) -> dict:
    return {
        "file":           os.path.basename(path),
        "manufacturer":   info["manufacturer"],
        "well":           info["well"],
        "sampling_hz":    info.get("sampling_rate", ""),
        "full_dur_s":     round(info["full_duration"], 3),
        "window_start_s": start_s,
        "window_dur_s":   round(info["window_dur"], 3),
        "total_ch":       info["total_ch"],
        "active_ch":      info["active_ch"],
        "ae_per_min":     ae_per_min,
    }


def write_summary(rows: list, out_dir: str):
    """Write the summary rows to <out_dir>/info_summary.tsv (folder created if needed)."""
    if not out_dir or not rows:
        return
    if not _PANDAS:
        print("pandas not available — skipping TSV export.")
        return
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, SUMMARY_NAME)
    pd.DataFrame(rows).to_csv(out_path, index=False, sep='\t')
    print(f"Summary saved: {out_path}")


# ── folder scan ─────────────────────────────────────────────────────────────────
# Only the manufacturer's own file type (loaders.EXTENSIONS) is scanned, so a
# mixed folder skips the other devices' files instead of failing on each of them.
_EXTENSIONS = {k: set(v) for k, v in EXTENSIONS.items()}


def scan_folder(folder: str, manufacturer: str, ae_per_min: float,
                start_s: float, dur_s: float, sr: float) -> list:
    exts = _EXTENSIONS[manufacturer]
    files = sorted(
        os.path.join(folder, f)
        for f in os.listdir(folder)
        if os.path.splitext(f)[1].lower() in exts
    )
    if not files:
        print(f"No {'/'.join(sorted(exts))} files found in {folder}")
        return []

    print(f"\nScanning {len(files)} file(s) in: {folder}")
    print("=" * 60)

    rows = []
    for path in files:
        try:
            info = scan_file(path, manufacturer, None, ae_per_min, start_s, dur_s, sr)
            _print_result(path, info)
            rows.append(_summary_row(path, info, start_s, ae_per_min))
        except Exception as e:
            print(f"  SKIPPED: {os.path.basename(path)} — not readable as "
                  f"--manufacturer {manufacturer} ({e})\n")
    return rows


# ── CLI ──────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        prog="cascade-info",
        description=(
            "CASCADE recording inspector.\n\n"
            "Quickly reports duration, total channels, and active electrodes\n"
            "for any MEA file without running the full analysis pipeline.\n\n"
            "Use --start and --dur to inspect a specific time window — useful\n"
            "when one file contains baseline + treatment:\n\n"
            "  cascade-info --input recording.bxr --start 0  --dur 60   # baseline\n"
            "  cascade-info --input recording.bxr --start 61 --dur 59   # post-treatment\n\n"
            "Scan a whole folder (give a folder to --input):\n"
            "  cascade-info --input /data/recordings/ --manufacturer cl --ae 5 --output info_out"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )

    # Flags are grouped by workflow step (same order as README "All flags").

    g = ap.add_argument_group("1 Setup — input / output")
    g.add_argument("--manufacturer",
                    required=True, choices=["3brain", "maxwell", "cl", "mcs", "axion"],
                    help="MEA manufacturer. 3brain (.bxr), maxwell (.h5), cl (.h5), mcs (.h5), axion (.csv)")
    g.add_argument("--input",
                     default=None, metavar="PATH",
                     help="A recording file (.bxr, .h5 or .csv), or a folder: every "
                          "file of the --manufacturer's type in it is scanned "
                          "(3brain .bxr; maxwell/cl/mcs .h5; axion .csv).")
    # pre-0.9.5 spelling of a folder --input; hidden, kept so old scripts still run
    g.add_argument("--folder", default=None, help=argparse.SUPPRESS)
    g.add_argument("--well",
                     default=None, metavar="NAME",
                     help="Well or stream name (e.g. Well_A1 for 3Brain, "
                          "data0000 for Maxwell). Default: first available; "
                          "for Axion, all wells merged.")
    g.add_argument("--sr",
                     type=float, default=None, metavar="HZ",
                     help="Sampling rate in Hz. Only needed when reading a CSV "
                          "file that uses a frame column instead of time_s.")
    g.add_argument("--output",
                     default=None, metavar="DIR",
                     help=f"Output folder: saves the summary table as {SUMMARY_NAME} "
                          "(one row per file). Created automatically if it does not "
                          "exist. Default: print to the terminal only.")

    g = ap.add_argument_group("2 Preprocessing — time window and active channels")
    g.add_argument("--start",
                     type=float, default=0.0, metavar="SEC",
                     help="Start of the window to inspect in seconds. "
                          "(default: 0 = beginning of recording)")
    g.add_argument("--dur",
                     type=float, default=None, metavar="SEC",
                     help="Duration of the window in seconds. "
                          "Default: full recording length.")
    g.add_argument("--ae",
                     type=float, default=0.0, metavar="N",
                     help="Count channels with >= N spikes per minute as active. "
                          "Set to 0 to count all channels with at least 1 spike. "
                          "(default: 0)")

    g = ap.add_argument_group("5 Misc — information (no analysis)")
    g.add_argument("--help", action="help",
                    help="Show this help message and exit.")
    g.add_argument("--version", action="version", version=f"CASCADE v{__version__}",
                    help="Show the CASCADE version and exit.")
    add_citation_argument(g)

    args = ap.parse_args()

    if args.folder:
        if args.input:
            ap.error("--folder is the old spelling of --input; give only --input.")
        args.input = args.folder
    if args.input is None:
        ap.error("--input is required (a recording file or a folder).")

    if os.path.isdir(args.input):
        rows = scan_folder(args.input, args.manufacturer, args.ae, args.start,
                           args.dur, args.sr)
    elif os.path.isfile(args.input):
        try:
            info = scan_file(args.input, args.manufacturer, args.well, args.ae,
                             args.start, args.dur, args.sr)
        except Exception as e:
            print("ERROR: " + read_error_message(args.input, args.manufacturer, e),
                  file=sys.stderr)
            sys.exit(1)
        print("=" * 60)
        _print_result(args.input, info)
        print("=" * 60)
        rows = [_summary_row(args.input, info, args.start, args.ae)]
    else:
        print(f"ERROR: not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    write_summary(rows, args.output)


if __name__ == "__main__":
    main()
