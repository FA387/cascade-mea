"""
loaders/cl1_loader.py
---------------------
Loader for Cortical Labs CL1 .h5 files.

Both the standard recording and the game recording are supported.
The game recording contains additional /data_stream/ groups (ball position,
paddle, motor regions, etc.) which are silently ignored — only /spikes is used.

HDF5 structure (PyTables v2 format)
------------------------------------
/ (root attrs)
    sampling_frequency    – int64, Hz  (identical to frames_per_second)
    channel_count         – int64, total channels on array (typically 64)
    duration_frames       – int64, recording length in frames
    duration_seconds      – float64, recording length in seconds
    system_id             – bytes, e.g. b'CL1-2544-058'
    uV_per_sample_unit    – float64, µV scale factor for waveform snippets

/spikes                   – PyTables TABLE, structured dtype:
    timestamp             – int64, frame index relative to recording start
    channel               – uint8,  0-based channel index (0 … channel_count-1)
    samples               – float32 (75,), spike waveform snippet in µV
                            (stored for reference; not used by CASCADE)

/stims                    – PyTables TABLE (stimulation timestamps; same layout
                            minus the samples field). Not used by CASCADE.

Game recordings additionally have:
/data_stream/gamestate/   – ball/paddle positions per frame  (ignored)
/data_stream/motor_regions/ – electrode→game-region mapping  (ignored)
"""

import numpy as np
import h5py

from .base_loader import BaseLoader
from cascade.utils.common import MEARecording, frames_to_seconds


class CL1Loader(BaseLoader):

    @staticmethod
    def manufacturer_name() -> str:
        return "Cortical Labs CL1"

    # ── file identification ────────────────────────────────────────────────────
    @staticmethod
    def is_cl1_file(path: str) -> bool:
        """Return True if the HDF5 file looks like a CL1 recording."""
        try:
            with h5py.File(path, "r") as h5:
                sid = h5.attrs.get("system_id", b"")
                if isinstance(sid, (bytes, np.bytes_)):
                    sid = sid.decode("ascii", errors="replace")
                return str(sid).upper().startswith("CL1")
        except Exception:
            return False

    # ── well listing (CL1 is single-well; exposed for API compatibility) ──────
    @staticmethod
    def list_wells(h5: h5py.File) -> list:
        """CL1 has no well concept.  Returns ['default'] for display purposes."""
        return ["default"]

    # ── main entry point ──────────────────────────────────────────────────────
    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        """
        Load spikes from a Cortical Labs CL1 .h5 file.

        The `well` parameter is accepted for API compatibility but is unused —
        CL1 is a single-well system.

        Parameters
        ----------
        path       : str   – path to the .h5 file
        well       : str   – ignored
        start_s    : float – start of analysis window in seconds (default 0)
        duration_s : float – duration of analysis window; None = full recording

        Returns
        -------
        MEARecording with spike_times in seconds, keyed by 0-based channel index.
        """
        with h5py.File(path, "r") as h5:

            # ── metadata ──────────────────────────────────────────────────────
            sr   = float(h5.attrs["sampling_frequency"])
            n_ch = int(h5.attrs["channel_count"])

            if duration_s is None:
                dur_attr = h5.attrs.get("duration_seconds", 0.0)
                duration_s = float(dur_attr) if float(dur_attr) > 0 else (
                    float(h5.attrs.get("duration_frames", 0)) / sr
                )

            # ── spikes table ──────────────────────────────────────────────────
            if "spikes" not in h5:
                raise KeyError(f"No /spikes table found in {path}")

            spikes = h5["spikes"][:]   # structured array: timestamp, channel, samples

            # ── metadata dict (read while file is open) ───────────────────────
            sid = h5.attrs.get("system_id", b"")
            if isinstance(sid, (bytes, np.bytes_)):
                sid = sid.decode("ascii", errors="replace")
            metadata = {
                "manufacturer" : self.manufacturer_name(),
                "file"         : path,
                "start_s"      : start_s,
                "system_id"    : str(sid),
                "channel_count": n_ch,
            }

        # ── empty recording guard ─────────────────────────────────────────────
        if spikes.size == 0:
            return MEARecording(
                sampling_rate = sr,
                n_channels    = n_ch,
                channel_ids   = list(range(n_ch)),
                spike_times   = {},
                duration      = float(duration_s),
                metadata      = metadata,
            )

        timestamps = spikes["timestamp"].astype(np.int64)
        channels   = spikes["channel"].astype(np.int32)

        # ── apply time window ─────────────────────────────────────────────────
        # timestamps are already in frames relative to recording start
        s_frame = int(round(start_s * sr))
        e_frame = int(round((start_s + duration_s) * sr))
        mask    = (timestamps >= s_frame) & (timestamps < e_frame)

        timestamps = timestamps[mask]
        channels   = channels[mask]

        # convert to seconds relative to window start
        spike_times_s = frames_to_seconds(timestamps - s_frame, sr)

        # ── build per-channel dict ────────────────────────────────────────────
        spike_dict = {}
        for ch in np.unique(channels):
            spike_dict[int(ch)] = spike_times_s[channels == ch].astype(np.float64)

        return MEARecording(
            sampling_rate = sr,
            n_channels    = n_ch,
            channel_ids   = list(range(n_ch)),
            spike_times   = spike_dict,
            duration      = float(duration_s),
            metadata      = metadata,
        )
