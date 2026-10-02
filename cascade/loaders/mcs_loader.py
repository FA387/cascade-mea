"""
loaders/mcs_loader.py
---------------------
Loader for Multi Channel Systems (MCS) HDF5 files produced by
Multi Channel Experimenter / Multi Channel DataManager.

HDF5 structure (MCS HDF5 Protocol v3)
---------------------------------------
/ (root attrs)
    McsHdf5ProtocolType    – b'RawData'
    McsHdf5ProtocolVersion – int32

/Data/Recording_<R>/AnalogStream/Stream_0/
    ChannelData               – int32 (n_channels × n_samples), raw voltage traces
    ChannelDataTimeStamps     – int64 (1 × 3): [FirstTimeStamp_µs, FirstIndex, LastIndex]
    InfoChannel               – structured array, one row per channel:
        ChannelID             – int32, 0-based channel index
        Label                 – bytes, electrode name (e.g. b'47')
        Tick                  – int64, µs per sample  (e.g. 50 → 20 000 Hz)
        Unit                  – bytes, physical unit of ChannelData (normally b'V')
        ConversionFactor      – int64, ADC → Unit scale factor
        Exponent              – int32, power-of-ten applied with ConversionFactor
                                value_in_Unit = raw_adc * ConversionFactor * 10**Exponent

/Data/Recording_<R>/SegmentStream/Stream_0/
    InfoSegment               – structured array (one row per segment/channel):
        SegmentID             – int32
        Label                 – bytes, electrode label
        SourceChannelIDs      – bytes, ChannelID of the source analog channel
    SegmentData_<N>           – int32 (n_waveform_samples × n_spikes), waveform cutouts
        @SourceChannelID      – int32, ChannelID in AnalogStream this segment came from
    SegmentData_ts_<N>        – int64 (1 × n_spikes), spike timestamps in µs
                                relative to recording start (t0 from ChannelDataTimeStamps)
"""

import re
import numpy as np
import h5py

from .base_loader import BaseLoader
from cascade.utils.common import MEARecording, frames_to_seconds


class MCSLoader(BaseLoader):

    @staticmethod
    def manufacturer_name() -> str:
        return "Multi Channel Systems"

    # ── file identification ────────────────────────────────────────────────────
    @staticmethod
    def is_mcs_file(path: str) -> bool:
        """Return True if the HDF5 file looks like an MCS recording."""
        try:
            with h5py.File(path, "r") as h5:
                return "McsHdf5ProtocolType" in h5.attrs
        except Exception:
            return False

    # ── recording listing (analogous to 'wells' in other loaders) ─────────────
    @staticmethod
    def list_wells(h5: h5py.File) -> list:
        """Return available recording names, e.g. ['Recording_0', 'Recording_1']."""
        if "Data" not in h5:
            return []
        return sorted(k for k in h5["Data"].keys() if k.startswith("Recording_"))

    # ── internal helpers ──────────────────────────────────────────────────────
    @staticmethod
    def _analog_path(recording: str) -> str:
        return f"Data/{recording}/AnalogStream/Stream_0"

    @staticmethod
    def _segment_path(recording: str) -> str:
        return f"Data/{recording}/SegmentStream/Stream_0"

    # ── main entry point ──────────────────────────────────────────────────────
    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        """
        Load spikes from an MCS Multi Channel Experimenter .h5 file.

        The `well` parameter selects the recording block (e.g. 'Recording_0').
        If None, the first available recording is used.

        Parameters
        ----------
        path       : str   – path to the .h5 file
        well       : str   – recording name (default: first available)
        start_s    : float – start of analysis window in seconds (default 0)
        duration_s : float – duration of analysis window; None = full recording

        Returns
        -------
        MEARecording with spike_times in seconds, keyed by ChannelID.
        """
        with h5py.File(path, "r") as h5:

            # ── recording selection ───────────────────────────────────────────
            recordings = self.list_wells(h5)
            if not recordings:
                raise ValueError(f"No Recording_N groups found in {path}")
            recording = well or recordings[0]
            if f"Data/{recording}" not in h5:
                raise ValueError(
                    f"Recording '{recording}' not found. Available: {recordings}"
                )

            a_path = self._analog_path(recording)
            s_path = self._segment_path(recording)

            # ── timing & channel info (from AnalogStream) ─────────────────────
            info_ch = h5[f"{a_path}/InfoChannel"][:]
            tick_us = int(info_ch[0]["Tick"])          # µs per sample
            sr      = 1e6 / tick_us                    # Hz

            ts_header = h5[f"{a_path}/ChannelDataTimeStamps"][:]
            t0_us     = int(ts_header[0, 0])           # recording start in µs (usually 0)

            n_samples  = h5[f"{a_path}/ChannelData"].shape[1]
            full_dur_s = n_samples * tick_us / 1e6

            if duration_s is None:
                duration_s = full_dur_s

            n_ch       = len(info_ch)
            channel_ids = [int(r["ChannelID"]) for r in info_ch]

            # ── build label map: ChannelID → electrode label ──────────────────
            label_map = {
                int(r["ChannelID"]): r["Label"].decode("ascii", errors="replace").strip()
                for r in info_ch
            }

            # ── ADC → µV scale per channel ────────────────────────────────────
            # MCS stores waveform cutouts as raw ADC integers. The physical value
            # is  raw * ConversionFactor * 10**Exponent  in `Unit` (always 'V' for
            # analog streams), so µV needs a further 1e6. Peak-to-peak is a
            # difference, so the ADZero offset cancels and is not needed.
            # If a file lacks these fields we return no amplitudes at all rather
            # than raw ADC counts, so --max_amplitude cannot silently compare a
            # µV threshold against ADC units.
            fields    = info_ch.dtype.names or ()
            amps_ok   = {"ConversionFactor", "Exponent"}.issubset(fields)
            scale_map = {}
            if amps_ok:
                for r in info_ch:
                    unit = (r["Unit"].decode("ascii", "replace").strip()
                            if "Unit" in fields else "V")
                    if unit != "V":
                        amps_ok = False
                        break
                    scale_map[int(r["ChannelID"])] = (
                        float(r["ConversionFactor"]) * (10.0 ** int(r["Exponent"])) * 1e6
                    )

            # ── spike timestamps + waveforms (from SegmentStream) ────────────
            spike_dict: dict = {}
            amp_dict:   dict = {}

            if s_path in h5:
                seg = h5[s_path]
                # find all SegmentData_ts_<N> datasets
                ts_keys = sorted(
                    k for k in seg.keys() if k.startswith("SegmentData_ts_")
                )

                # time window in µs relative to t0
                s_us = int(round(start_s * 1e6))
                e_us = int(round((start_s + duration_s) * 1e6))

                for ts_key in ts_keys:
                    # matching waveform key (just drop "_ts")
                    seg_key = ts_key.replace("SegmentData_ts_", "SegmentData_")
                    if seg_key not in seg:
                        continue

                    src_ch = int(seg[seg_key].attrs["SourceChannelID"])

                    # timestamps: shape (1, n_spikes) in µs absolute
                    ts_us = seg[ts_key][0].astype(np.int64)   # squeeze row dim

                    # waveform cutouts: shape (n_waveform_samples, n_spikes)
                    waveforms = seg[seg_key][:]
                    ptp = (waveforms.max(axis=0)
                           - waveforms.min(axis=0)).astype(np.float64)
                    if amps_ok:
                        # 0.0 for an unmapped channel: fail open (never filtered)
                        # rather than risk discarding a channel on a bogus scale.
                        ptp *= scale_map.get(src_ch, 0.0)   # ADC → µV

                    # make relative to t0, apply time window
                    ts_rel  = ts_us - t0_us
                    mask    = (ts_rel >= s_us) & (ts_rel < e_us)
                    ts_win  = ts_rel[mask]
                    ptp_win = ptp[mask]

                    if ts_win.size == 0:
                        continue

                    # convert to seconds relative to window start
                    spike_times_s = (ts_win - s_us) / 1e6

                    if src_ch in spike_dict:
                        # sort the merged times AND carry the amplitudes through
                        # the same permutation, so spike i still pairs with
                        # amplitude i after two segments are combined.
                        merged_t = np.concatenate([spike_dict[src_ch], spike_times_s])
                        merged_a = np.concatenate([amp_dict[src_ch], ptp_win])
                        order    = np.argsort(merged_t, kind="stable")
                        spike_dict[src_ch] = merged_t[order].astype(np.float64)
                        amp_dict[src_ch]   = merged_a[order].astype(np.float32)
                    else:
                        order = np.argsort(spike_times_s, kind="stable")
                        spike_dict[src_ch] = spike_times_s[order].astype(np.float64)
                        amp_dict[src_ch]   = ptp_win[order].astype(np.float32)

            metadata = {
                "manufacturer" : self.manufacturer_name(),
                "file"         : path,
                "recording"    : recording,
                "start_s"      : start_s,
                "label_map"    : label_map,
                "amplitude_unit": "uV" if amps_ok else None,
            }

        return MEARecording(
            sampling_rate    = sr,
            n_channels       = n_ch,
            channel_ids      = channel_ids,
            spike_times      = spike_dict,
            duration         = float(duration_s),
            metadata         = metadata,
            spike_amplitudes = amp_dict if (amp_dict and amps_ok) else None,
        )
