"""
loaders/brainwave_loader.py
---------------------------
Loader for 3Brain BrainWave .bxr files.

HDF5 structure expected
-----------------------
/ (root attrs)
    SamplingRate              – float, Hz
    ExperimentDateTimeUtc     – str
    PlateModel                – str

/Well_A1/
    SpikeTimes   – int array, frame indices
    SpikeChIdxs  – int array, channel index per spike
    SpikeUnits   – int array (optional), unit ID per spike
    SpikeForms   – float array (optional), waveform data
    StoredChIdxs – int array, active channel indices
    TOC          – (N×2) array, [[start_frame, end_frame], …]
"""

import numpy as np
import h5py
from .base_loader import BaseLoader
from cascade.utils.common import MEARecording, frames_to_seconds


class BrainwaveLoader(BaseLoader):

    @staticmethod
    def manufacturer_name() -> str:
        return "3Brain BrainWave"

    # ── well listing ──────────────────────────────────────────────────────────
    @staticmethod
    def list_wells(h5: h5py.File) -> list:
        return sorted([k for k in h5.keys() if k.startswith("Well_")])

    # ── duration from TOC ─────────────────────────────────────────────────────
    @staticmethod
    def _recording_duration(h5: h5py.File, sr: float) -> float:
        if "TOC" in h5:
            toc = np.array(h5["TOC"])
            if toc.size >= 2:
                return float(toc[-1, 1] - toc[0, 0]) / sr
        return None

    # ── main entry point ──────────────────────────────────────────────────────
    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        """
        Load spikes from a 3Brain .bxr file.

        Returns
        -------
        MEARecording with spike_times in seconds, keyed by channel index.
        """
        with h5py.File(path, "r") as h5:
            # sampling rate
            sr = float(h5.attrs["SamplingRate"])

            # well selection
            wells = self.list_wells(h5)
            if not wells:
                raise ValueError(f"No wells found in {path}")
            well = well or wells[0]
            if well not in h5:
                raise ValueError(f"Well '{well}' not found. Available: {wells}")

            # recording duration
            if duration_s is None:
                duration_s = self._recording_duration(h5, sr) or 60.0

            # raw arrays
            spike_frames = h5[f"{well}/SpikeTimes"][:]
            spike_ch     = h5[f"{well}/SpikeChIdxs"][:]

            # stored (active) channels
            stored_ch = (h5[f"{well}/StoredChIdxs"][:]
                         if f"{well}/StoredChIdxs" in h5
                         else np.unique(spike_ch))

            # apply time window
            s_frame = int(round(start_s * sr))
            e_frame = int(round((start_s + duration_s) * sr))
            mask    = (spike_frames >= s_frame) & (spike_frames < e_frame)

            spike_frames = spike_frames[mask]
            spike_ch     = spike_ch[mask]

            # convert to seconds relative to window start
            spike_times_s = frames_to_seconds(spike_frames - s_frame, sr)

            # build per-channel dict
            spike_dict = {}
            for ch in np.unique(spike_ch):
                spike_dict[int(ch)] = spike_times_s[spike_ch == ch].astype(np.float64)

            # metadata
            metadata = {
                "manufacturer" : self.manufacturer_name(),
                "file"         : path,
                "well"         : well,
                "start_s"      : start_s,
                "has_waveforms": f"{well}/SpikeForms" in h5,
            }
            for attr_key in ("ExperimentDateTimeUtc", "PlateModel", "SamplingRate"):
                if attr_key in h5.attrs:
                    metadata[attr_key] = h5.attrs[attr_key]

        return MEARecording(
            sampling_rate = sr,
            n_channels    = int(max(stored_ch)) + 1 if len(stored_ch) > 0 else 0,
            channel_ids   = [int(c) for c in sorted(stored_ch)],
            spike_times   = spike_dict,
            duration      = float(duration_s),
            metadata      = metadata,
        )

    # ── optional: load mean waveforms ─────────────────────────────────────────
    def load_waveforms(self, path: str, well: str = None,
                       start_s: float = 0.0, duration_s: float = 60.0,
                       max_spikes: int = 200) -> dict:
        """
        Load mean spike waveforms from the SpikeForms dataset.

        Returns
        -------
        dict  {channel_id: {"mean": np.ndarray, "std": np.ndarray, "t_ms": np.ndarray}}
        or empty dict if no waveform data is present.
        """
        result = {}
        with h5py.File(path, "r") as h5:
            sr    = float(h5.attrs["SamplingRate"])
            wells = self.list_wells(h5)
            well  = well or wells[0]

            if f"{well}/SpikeForms" not in h5:
                return result

            spike_frames = h5[f"{well}/SpikeTimes"][:]
            spike_ch     = h5[f"{well}/SpikeChIdxs"][:]
            wl           = int(h5[f"{well}/SpikeForms"].attrs.get("Wavelength", 0))
            if wl <= 0:
                return result

            s_frame = int(round(start_s * sr))
            e_frame = int(round((start_s + duration_s) * sr))
            mask    = (spike_frames >= s_frame) & (spike_frames < e_frame)
            idxs    = np.where(mask)[0]
            if idxs.size == 0:
                return result

            forms_all = h5[f"{well}/SpikeForms"][:]
            if forms_all.shape[0] < (idxs[-1] + 1) * wl:
                return result

            use     = idxs[:max_spikes]
            waves   = np.vstack([forms_all[i * wl:(i + 1) * wl] for i in use])
            ch_use  = spike_ch[use]
            t_ms    = np.arange(wl) / sr * 1000.0

            for ch in np.unique(ch_use):
                W = waves[ch_use == ch]
                result[int(ch)] = {
                    "mean" : W.mean(axis=0),
                    "std"  : W.std(axis=0, ddof=1) if W.shape[0] > 1 else np.zeros(wl),
                    "t_ms" : t_ms,
                }
        return result
