"""
loaders/finalspark_loader.py
----------------------------
Loader for FinalSpark Neuroplatform whole-life exports (fsNNN_package.hdf5 /
fsNNN_raw.hdf5).

Export layout (pandas HDFStore "table" format, read here with h5py)
-------------------------------------------------------------------
fsNNN_package.hdf5
    /fsNNN_wholelife_events/table    one row per detected spike
        time_of_event   int64, ns since 1970 (UTC)
        electrode       int64, 0-127 = 32 * MEA + 8 * site + n (one MEA per export)
        values_block_0  float64 (1,), max_voltage_uv (signed single peak)
    /fsNNN_wholelife_stimulations/table   time_of_stim (ns) -> metadata["stim_times"]
    incubator tables                      (not used)
fsNNN_raw.hdf5
    /fsNNN_wholelife_raw             3 ms, 30 kHz snippets around each event
                                     (1 ms before, 2 ms after the crossing)

--input is the package file. The raw file only holds snippets of the same
events, so it adds no spikes and is not read.

Time 0 is the first event, as in FinalSpark's demo_utils. Wells: the 4x8 MEA
holds 4 sites of 8 electrodes, one neurosphere group per site; --well 0-3 keeps electrodes
base+8*w ... base+8*w+7 (base = 32 * MEA, as in np-utils), no --well keeps all 32.
Channel ids are FinalSpark's own electrode numbers.
"""

import os
import numpy as np
import h5py

from .base_loader import BaseLoader
from cascade.utils.common import MEARecording

SAMPLING_RATE = 30000.0   # Intan headstage, per the export README
SITE_SIZE = 8


def _events_key(h5) -> str | None:
    return next((k for k in h5.keys() if k.endswith("_wholelife_events")), None)


def is_finalspark_file(path: str) -> bool:
    """Format check only (reads the table names): true for both files of an export."""
    try:
        with h5py.File(path, "r") as h5:
            return any(k.endswith(("_wholelife_events", "_wholelife_raw")) for k in h5.keys())
    except Exception:
        return False


def events_file(path: str) -> str:
    """Check that *path* is the package file (the one with the spike table)."""
    if os.path.basename(path).endswith("_raw.hdf5"):
        raise ValueError(
            f"{os.path.basename(path)} only holds 3 ms snippets; the spikes are in "
            f"the *_package.hdf5 file of the same export. Pass that file to --input.")
    with h5py.File(path, "r") as h5:
        if _events_key(h5):
            return path
    raise KeyError(f"no *_wholelife_events table in {os.path.basename(path)}")


class FinalSparkLoader(BaseLoader):

    @staticmethod
    def manufacturer_name() -> str:
        return "FinalSpark"

    @staticmethod
    def list_wells(path: str = None) -> list:
        return ["all", "0", "1", "2", "3"]

    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        ev_path = events_file(path)
        with h5py.File(ev_path, "r") as h5:
            key = _events_key(h5)
            table = h5[key]["table"]
            # values_block_N is named in its *_kind attribute (pickled column list)
            amp_field = next(
                n for n in table.dtype.names
                if b"max_voltage_uv" in bytes(table.attrs.get(f"{n}_kind", b"")))
            t_ns = table["time_of_event"]
            elec = table["electrode"].astype(np.int32)
            amp = np.abs(table[amp_field][:, 0]).astype(np.float32)
            stim_key = key.replace("_wholelife_events", "_wholelife_stimulations")
            stim_ns = (np.sort(h5[stim_key]["table"]["time_of_stim"])
                       if stim_key in h5 else np.zeros(0, np.int64))

        fs_name = key[:-len("_wholelife_events")]
        base = 32 * (int(elec.min()) // 32) if elec.size else 0   # first electrode of this MEA
        t0 = int(t_ns.min()) if t_ns.size else 0
        t = (t_ns - t0) / 1e9
        span = float(t.max()) if t.size else 0.0

        if well in (None, "", "all"):
            label, channels = fs_name, list(range(base, base + 4 * SITE_SIZE))
        else:
            if str(well) not in ("0", "1", "2", "3"):
                raise ValueError(f"--well for FinalSpark is 0, 1, 2 or 3 (site of 8 electrodes), got '{well}'")
            w = int(well)
            label, channels = f"{fs_name}_site{w}", list(range(base + w * SITE_SIZE, base + (w + 1) * SITE_SIZE))

        if duration_s is None:
            duration_s = max(span - start_s, 0.0) + 1 / SAMPLING_RATE  # keep the last event
        keep = (t >= start_s) & (t < start_s + duration_s) & np.isin(elec, channels)
        t, elec, amp = t[keep] - start_s, elec[keep], amp[keep]

        # stimulation times (any electrode) on the same clock; 1 s before the window
        # is kept so a stim just before --start still blanks the first spikes
        st = (stim_ns - t0) / 1e9 - start_s
        stim_times = st[(st >= -1.0) & (st < duration_s)]

        order = np.lexsort((t, elec))
        t, elec, amp = t[order], elec[order], amp[order]
        chans, first = np.unique(elec, return_index=True)
        spike_times = dict(zip(chans.tolist(), np.split(t, first[1:])))
        amps = dict(zip(chans.tolist(), np.split(amp, first[1:])))

        return MEARecording(
            sampling_rate    = SAMPLING_RATE,
            n_channels       = len(channels),
            channel_ids      = channels,
            spike_times      = spike_times,
            duration         = float(duration_s),
            metadata         = {
                "manufacturer": self.manufacturer_name(),
                "file": path,
                "events_file": ev_path,
                "well": label,
                "start_s": start_s,
                "t0_utc": np.datetime64(t0, "ns").astype(str) + "Z",
                "stim_times": stim_times,   # s, window clock; used by --stim_blank
            },
            spike_amplitudes = amps,
        )
