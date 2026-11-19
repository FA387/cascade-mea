"""CASCADE manufacturer loaders - v0.9.6"""
from .brainwave_loader import BrainwaveLoader
from .maxwell_loader   import MaxwellLoader
from .cl_loader        import CLLoader
from .mcs_loader       import MCSLoader
from .axion_loader     import AxionLoader
from .finalspark_loader import FinalSparkLoader


import os as _os

# --manufacturer name -> loader class. The single registration point for a new manufacturer.
LOADERS = {
    "3brain":     BrainwaveLoader,
    "maxwell":    MaxwellLoader,
    "cl":         CLLoader,
    "mcs":        MCSLoader,
    "axion":      AxionLoader,
    "finalspark": FinalSparkLoader,
}


def get_loader(manufacturer: str):
    if manufacturer not in LOADERS:
        raise ValueError(f"Unknown manufacturer '{manufacturer}'. "
                         f"Supported: {', '.join(LOADERS)}")
    return LOADERS[manufacturer]()


# File types each manufacturer writes; .h5 is shared by three devices.
EXTENSIONS = {
    "3brain":  (".bxr",),
    "maxwell": (".h5", ".hdf5"),
    "cl":      (".h5", ".hdf5"),
    "mcs":     (".h5", ".hdf5"),
    "axion":   (".csv",),
    "finalspark": (".hdf5", ".h5"),
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
    import h5py
    try:
        if not h5py.is_hdf5(path):
            return "unknown"
    except OSError:  # missing or unreadable file
        return "unknown"
    from .finalspark_loader import is_finalspark_file
    if is_finalspark_file(path):
        return "finalspark"
    if CLLoader.is_cl_file(path):
        return "cl"
    if MCSLoader.is_mcs_file(path):
        return "mcs"
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
