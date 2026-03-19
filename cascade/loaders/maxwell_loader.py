"""
loaders/maxwell_loader.py
--------------------------
Loader for Maxwell Biosystems .h5 files (data.raw.h5).

HDF5 structure expected
-----------------------
/data_store/<stream>/
    spikes                    – structured array with fields:
                                  frameno   (int64)
                                  channel   (int32)
                                  amplitude (float32, optional)
    settings/
        sampling              – float scalar, Hz
        gain                  – float scalar
        lsb                   – float scalar
        mapping               – structured array (channel → electrode mapping)

/assay/inputs/
    record_time               – float scalar, seconds
    electrodes                – int, number of active electrodes
    spike_only                – bool
"""

import numpy as np
import h5py
from .base_loader import BaseLoader
from cascade.utils.common import MEARecording, frames_to_seconds


class MaxwellLoader(BaseLoader):

    @staticmethod
    def manufacturer_name() -> str:
        return "Maxwell Biosystems"

    # ── stream listing ────────────────────────────────────────────────────────
    @staticmethod
    def list_wells(h5: h5py.File) -> list:
        """Maxwell calls them 'streams'; we expose them as 'wells' to the pipeline."""
        if "data_store" not in h5:
            return []
        return sorted([k for k in h5["data_store"].keys() if k.startswith("data")])

    # ── sampling rate ─────────────────────────────────────────────────────────
    @staticmethod
    def _get_sampling_rate(h5: h5py.File, stream: str) -> float:
        path = f"data_store/{stream}/settings/sampling"
        if path not in h5:
            raise KeyError(f"Sampling rate not found at {path}")
        return float(np.array(h5[path]).flat[0])

    # ── recording duration ────────────────────────────────────────────────────
    @staticmethod
    def _get_duration(h5: h5py.File) -> float | None:
        path = "assay/inputs/record_time"
        if path not in h5:
            return None
        val = np.array(h5[path])
        return float(val.flat[0])

    # ── time-zero frame ───────────────────────────────────────────────────────
    @staticmethod
    def _get_start_frame(h5: h5py.File, stream: str) -> int:
        """
        Determine the first frame index for the stream so spike times
        are expressed relative to the recording start (not absolute frame 0).
        Tries routed frame_nos first, then falls back to first spike.
        """
        fn_path = f"data_store/{stream}/groups/routed/frame_nos"
        if fn_path in h5:
            fn = np.array(h5[fn_path])
            if fn.size > 0:
                return int(np.min(fn))

        sp_path = f"data_store/{stream}/spikes"
        if sp_path in h5:
            sp = h5[sp_path][:]
            if sp.size > 0 and "frameno" in sp.dtype.fields:
                return int(np.min(sp["frameno"]))

        return 0

    # ── main entry point ──────────────────────────────────────────────────────
    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        """
        Load spikes from a Maxwell .h5 file.

        Returns
        -------
        MEARecording with spike_times in seconds, keyed by channel index.
        """
        with h5py.File(path, "r") as h5:
            # stream selection
            streams = self.list_wells(h5)
            if not streams:
                raise ValueError(f"No data streams found in {path} "
                                 f"(expected /data_store/data0000, …)")
            stream = well or streams[0]
            if f"data_store/{stream}" not in h5:
                raise ValueError(f"Stream '{stream}' not found. Available: {streams}")

            # sampling rate
            sr = self._get_sampling_rate(h5, stream)

            # recording duration
            if duration_s is None:
                duration_s = self._get_duration(h5) or 60.0

            # spikes dataset
            sp_path = f"data_store/{stream}/spikes"
            if sp_path not in h5:
                raise KeyError(f"Spikes dataset not found at {sp_path}")

            sp = h5[sp_path][:]

            if sp.size == 0:
                # empty recording → return empty MEARecording
                return MEARecording(
                    sampling_rate = sr,
                    n_channels    = 0,
                    channel_ids   = [],
                    spike_times   = {},
                    duration      = float(duration_s),
                    metadata      = {"manufacturer": self.manufacturer_name(),
                                     "file": path, "well": stream},
                )

            start_frame0  = self._get_start_frame(h5, stream)
            frames_abs    = sp["frameno"].astype(np.int64)
            channels      = sp["channel"].astype(np.int32)

            # apply time window
            s_frame = start_frame0 + int(round(start_s * sr))
            e_frame = start_frame0 + int(round((start_s + duration_s) * sr))
            mask    = (frames_abs >= s_frame) & (frames_abs < e_frame)

            frames_rel    = frames_abs[mask] - start_frame0
            channels_sel  = channels[mask]

            spike_times_s = frames_to_seconds(frames_rel, sr)

            # build per-channel dict
            spike_dict = {}
            for ch in np.unique(channels_sel):
                spike_dict[int(ch)] = spike_times_s[channels_sel == ch].astype(np.float64)

            # active channel list from mapping table (if present)
            map_path = f"data_store/{stream}/settings/mapping"
            if map_path in h5:
                mapping = h5[map_path][:]
                channel_ids = [int(r["channel"]) for r in mapping
                               if "channel" in mapping.dtype.names]
            else:
                channel_ids = [int(c) for c in np.unique(channels_sel)]

            # metadata
            metadata = {
                "manufacturer": self.manufacturer_name(),
                "file"        : path,
                "well"        : stream,
                "start_s"     : start_s,
            }
            for key in ("gain", "lsb"):
                p = f"data_store/{stream}/settings/{key}"
                if p in h5:
                    metadata[key] = float(np.array(h5[p]).flat[0])

        return MEARecording(
            sampling_rate = sr,
            n_channels    = max(channel_ids) + 1 if channel_ids else 0,
            channel_ids   = sorted(channel_ids),
            spike_times   = spike_dict,
            duration      = float(duration_s),
            metadata      = metadata,
        )
