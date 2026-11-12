"""CASCADE manufacturer loaders - v0.9.5"""
from .brainwave_loader import BrainwaveLoader
from .maxwell_loader   import MaxwellLoader
from .cl_loader        import CLLoader
from .mcs_loader       import MCSLoader
from .axion_loader     import AxionLoader


import os as _os

# File types each manufacturer writes; .h5 is shared by three devices.
EXTENSIONS = {
    "3brain":  (".bxr",),
    "maxwell": (".h5", ".hdf5"),
    "cl":      (".h5", ".hdf5"),
    "mcs":     (".h5", ".hdf5"),
    "axion":   (".csv",),
}


def detect_manufacturer(path: str) -> str:
    """Best guess of the --manufacturer for a file, or "unknown".

    .h5 files are told apart by their root attributes: Cortical Labs writes a
    system_id starting with "CL1", MCS writes McsHdf5ProtocolType; any other
    .h5 is assumed to be Maxwell.
    """
    ext = _os.path.splitext(path)[1].lower()
    if ext == ".bxr":
        return "3brain"
    if ext == ".csv":
        return "axion"
    if ext not in (".h5", ".hdf5"):
        return "unknown"
    try:
        import h5py
        with h5py.File(path, "r") as h5:
            sid = h5.attrs.get("system_id", b"")
            if isinstance(sid, bytes):
                sid = sid.decode("ascii", errors="replace")
            if str(sid).upper().startswith("CL1"):
                return "cl"
            if "McsHdf5ProtocolType" in h5.attrs:
                return "mcs"
    except Exception:
        return "unknown"
    return "maxwell"


def read_error_message(path: str, manufacturer: str, exc: Exception) -> str:
    """Plain-language message for a file the chosen loader could not read."""
    msg = (f"Could not read {_os.path.basename(path)} as --manufacturer "
           f"{manufacturer} ({type(exc).__name__}: {exc}).")
    guess = detect_manufacturer(path)
    if guess not in ("unknown", manufacturer):
        msg += (f"\nThis file looks like a '{guess}' recording - try "
                f"--manufacturer {guess}.")
    elif _os.path.splitext(path)[1].lower() not in EXTENSIONS.get(manufacturer, ()):
        msg += (f"\n--manufacturer {manufacturer} expects "
                f"{' or '.join(EXTENSIONS[manufacturer])} files.")
    return msg
