"""
loaders/base_loader.py
----------------------
Abstract base class that every manufacturer loader must inherit from.

To add a new manufacturer:
1. Create loaders/<name>_loader.py
2. Subclass BaseLoader
3. Implement load()  — that's it.
   All analysis code in core/ will work automatically.
"""

from abc import ABC, abstractmethod
from cascade.utils.common import MEARecording


class BaseLoader(ABC):
    """
    Abstract loader.  Every manufacturer-specific loader inherits this class
    and implements the load() method to return a MEARecording.
    """

    @abstractmethod
    def load(self, path: str, well: str = None,
             start_s: float = 0.0, duration_s: float = None) -> MEARecording:
        """
        Read raw data from *path* and return a normalised MEARecording.

        Parameters
        ----------
        path       : str   – full path to the data file
        well       : str   – well/stream identifier (manufacturer-specific)
                             If None, the first available well is used.
        start_s    : float – start of the analysis window in seconds (default 0)
        duration_s : float – duration of the analysis window in seconds
                             If None, use the full recording length.

        Returns
        -------
        MEARecording  – normalised, manufacturer-independent recording object
        """
        ...

    @staticmethod
    def manufacturer_name() -> str:
        """Human-readable name shown in reports (override in subclass)."""
        return "Unknown"
