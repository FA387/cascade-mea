"""
utils/common.py
---------------
Shared data structures, logging helpers, and figure saving.
All other modules import from here.
"""

import sys
import time
import dataclasses
from typing import Dict, List
import numpy as np
import matplotlib.pyplot as plt

# ── Colour codes for terminal output ──────────────────────────────────────────
RESET = "\033[0m";  BOLD  = "\033[1m"
CYAN  = "\033[36m"; MAG   = "\033[35m"
YEL   = "\033[33m"; GRN   = "\033[32m"
RED   = "\033[31m"; BLU   = "\033[34m"


# ── Normalised recording object (every loader must produce this) ───────────────
@dataclasses.dataclass
class MEARecording:
    """
    Manufacturer-independent representation of one MEA recording.

    Fields
    ------
    sampling_rate  : float  – samples per second (Hz)
    n_channels     : int    – total number of channels on the array
    channel_ids    : list   – list of channel identifiers (ints or strings)
    spike_times    : dict   – {channel_id: np.ndarray of spike times in SECONDS}
    duration       : float  – recording duration in seconds
    metadata       : dict   – manufacturer-specific extras (well name, gain, …)
    """
    sampling_rate : float
    n_channels    : int
    channel_ids   : List
    spike_times   : Dict[int, np.ndarray]   # channel_id → times in seconds
    duration      : float
    metadata      : dict

    def all_spike_times(self) -> np.ndarray:
        """Return every spike time as one flat, sorted array (seconds)."""
        if not self.spike_times:
            return np.array([], dtype=np.float64)
        return np.sort(np.concatenate(list(self.spike_times.values())))

    def n_spikes(self) -> int:
        return int(self.all_spike_times().size)

    def active_channels(self) -> list:
        """Return only channels that fired at least once."""
        return [ch for ch, st in self.spike_times.items() if len(st) > 0]


# ── Figure saving ──────────────────────────────────────────────────────────────
def save_fig(fig: plt.Figure, out_png: str, save_svg: bool = False,
             save_jpg: bool = False, dpi: int = 300):
    """
    Save figure as PNG (always, 300 dpi) and optionally as SVG and/or JPG.

    Parameters
    ----------
    fig      : matplotlib Figure
    out_png  : full path ending in .png
    save_svg : if True, also save .svg with the same base name
    save_jpg : if True, also save .jpg at the same dpi
    dpi      : resolution for raster formats
    """
    fig.savefig(out_png, dpi=dpi)
    if save_svg and out_png.lower().endswith(".png"):
        fig.savefig(out_png[:-4] + ".svg", format="svg")
    if save_jpg and out_png.lower().endswith(".png"):
        fig.savefig(out_png[:-4] + ".jpg", dpi=dpi)


# ── Terminal logging ───────────────────────────────────────────────────────────
def say(msg: str, icon: str = "[INFO]", color: str = CYAN,
        end: str = "\n", enabled: bool = True):
    if not enabled:
        return
    ts = time.strftime("%H:%M:%S")
    sys.stdout.write(f"{color}{icon} [{ts}] {msg}{RESET}{end}")
    sys.stdout.flush()


def banner(text: str, enabled: bool = True):
    if not enabled:
        return
    line = "═" * (len(text) + 2)
    print(f"{MAG}{BOLD}╔{line}╗\n║ {text} ║\n╚{line}╝{RESET}")


def step(title: str, enabled: bool = True):
    if enabled:
        print(f"{BLU}{BOLD}—— {title} ——{RESET}")


# ── Unit conversion helper ─────────────────────────────────────────────────────
def frames_to_seconds(frames: np.ndarray, sampling_rate: float) -> np.ndarray:
    """Convert raw frame indices to seconds."""
    return np.asarray(frames, dtype=np.float64) / float(sampling_rate)
