"""
loaders/axion_loader.py
-----------------------
Loader for Axion BioSystems Maestro MEA recordings.

Axion AxIS exports spike data as a CSV with an interleaved layout:
  - Columns 0-1: metadata key-value pairs (left side)
  - Columns 2-4: Time (s), Electrode, Amplitude (mV) (right side)
  Both sides share the same rows; metadata rows and spike rows co-exist.

Supported export files:
  *_spike_list.csv  – all detected spikes (preferred)
  *(000).csv        – Axion main export (same format, subset of spikes)

Electrode labels follow the pattern WELL_ELECTRODE (e.g. 'B1_36').
The `well` parameter filters to a single well (e.g. 'B1').
If `well` is None, spikes from all wells are merged into one recording.
"""

import os
import re
import numpy as np
import pandas as pd

from .base_loader import BaseLoader
from cascade.utils.common import MEARecording


_SR_RE = re.compile(r"([\d.]+)\s*kHz", re.IGNORECASE)


def _parse_sr(val: str) -> float:
    """Parse '12.5 kHz' → 12500.0.  Falls back to the Maestro default."""
    m = _SR_RE.search(str(val))
    if m:
        return float(m.group(1)) * 1000.0
    return 12500.0   # Maestro hardware default


class AxionLoader(BaseLoader):

    @staticmethod
    def manufacturer_name() -> str:
        return "Axion BioSystems"

    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        """
        Load spikes from an Axion AxIS spike-list CSV.

        Parameters
        ----------
        path       : str   – path to *_spike_list.csv or *(000).csv
        well       : str   – well to filter, e.g. 'B1'; None = all wells merged
        start_s    : float – start of analysis window in seconds (default 0)
        duration_s : float – duration of window in seconds; None = full recording

        Returns
        -------
        MEARecording with spike_times keyed by electrode label (str).
        """
        # ── read file ─────────────────────────────────────────────────────────────
        # Axion CSV files end with a "Well Information" footer that has a
        # different number of columns than the 5-column spike data section.
        # We read the raw text to find the split point, then parse only the
        # spike data section with pandas.
        with open(path, encoding="utf-8-sig") as fh:
            raw_lines = fh.readlines()

        # find where the Well Information footer begins (if present)
        split_idx = len(raw_lines)
        for i, line in enumerate(raw_lines):
            if line.strip().startswith("Well Information"):
                split_idx = i
                break

        import io
        spike_section = "".join(raw_lines[:split_idx])
        df = pd.read_csv(io.StringIO(spike_section), header=0,
                         dtype=str, low_memory=False)

        # parse the Well Information footer for metadata (optional)
        well_info: dict = {}
        if split_idx < len(raw_lines):
            import csv as csv_mod
            footer_lines = raw_lines[split_idx + 1:]   # skip "Well Information" header
            if footer_lines:
                reader = csv_mod.reader(footer_lines)
                for row in reader:
                    if len(row) >= 2 and row[0].strip():
                        well_info[row[0].strip()] = row[1:]

        if df.shape[1] < 5:
            raise ValueError(
                f"{path} has only {df.shape[1]} columns; "
                "expected the Axion 5-column interleaved format."
            )

        # column names after reading the header row
        time_col = df.columns[2]   # 'Time (s)'
        elec_col = df.columns[3]   # 'Electrode'
        amp_col  = df.columns[4]   # 'Amplitude (mV)'

        # ── extract metadata from left two columns ────────────────────────────────
        meta_raw: dict = {}
        for _, row in df.iterrows():
            k = str(row.iloc[0]).strip()
            v = str(row.iloc[1]).strip()
            if k and k not in ("nan", "") and v and v not in ("nan", ""):
                meta_raw[k] = v

        sr = 12500.0
        for k, v in meta_raw.items():
            if "sampling" in k.lower() or "frequency" in k.lower():
                sr = _parse_sr(v)
                break

        recording_name = meta_raw.get(
            "Recording Name", os.path.splitext(os.path.basename(path))[0]
        )
        plate_type = next(
            (v for k, v in meta_raw.items() if "plate type" in k.lower()), ""
        )

        # ── extract spike rows from right three columns ───────────────────────────
        t_num     = pd.to_numeric(df[time_col], errors="coerce")
        spike_idx = t_num.notna()

        t_arr    = t_num[spike_idx].to_numpy(dtype=np.float64)
        elec_arr = df.loc[spike_idx, elec_col].to_numpy(dtype=str)
        # amplitudes: abs value, mV → µV
        amp_arr  = (np.abs(pd.to_numeric(df[amp_col], errors="coerce")[spike_idx]
                           .to_numpy(dtype=np.float64)) * 1000.0)

        if t_arr.size == 0:
            raise ValueError(f"No spike rows found in {path}.")

        # ── well filter ───────────────────────────────────────────────────────────
        if well is not None:
            # electrode labels: 'B1_36' → prefix 'B1_'
            keep = np.array(
                [e == well or e.startswith(well + "_") for e in elec_arr]
            )
            if not keep.any():
                available = sorted({e.split("_")[0] for e in elec_arr})
                raise ValueError(
                    f"Well '{well}' not found. Available wells: {available}"
                )
            t_arr    = t_arr[keep]
            elec_arr = elec_arr[keep]
            amp_arr  = amp_arr[keep]

        # ── determine full recording duration before windowing ────────────────────
        full_dur = float(t_arr.max()) if t_arr.size > 0 else 0.0
        if duration_s is None:
            duration_s = full_dur - start_s

        # ── time window filter ────────────────────────────────────────────────────
        end_s    = start_s + duration_s
        win_mask = (t_arr >= start_s) & (t_arr < end_s)
        t_win    = t_arr[win_mask] - start_s   # times relative to window start
        elec_win = elec_arr[win_mask]
        amp_win  = amp_arr[win_mask]

        # ── build spike dict ──────────────────────────────────────────────────────
        # channel_ids = every electrode seen in the full (possibly well-filtered)
        # data, even if silent in the requested window.
        all_electrodes = sorted(np.unique(elec_arr).tolist())

        spike_dict: dict = {}
        amp_dict:   dict = {}
        for elec in np.unique(elec_win):
            m          = elec_win == elec
            sort_order = np.argsort(t_win[m])
            spike_dict[elec] = t_win[m][sort_order].astype(np.float64)
            amp_dict[elec]   = amp_win[m][sort_order].astype(np.float32)

        metadata = {
            "manufacturer"   : self.manufacturer_name(),
            "file"           : path,
            "recording_name" : recording_name,
            "plate_type"     : plate_type,
            "well"           : well or "all",
            "start_s"        : start_s,
            "raw_metadata"   : meta_raw,
            "well_info"      : well_info,
        }

        return MEARecording(
            sampling_rate    = sr,
            n_channels       = len(all_electrodes),
            channel_ids      = all_electrodes,
            spike_times      = spike_dict,
            duration         = float(duration_s),
            metadata         = metadata,
            spike_amplitudes = amp_dict if amp_dict else None,
        )
