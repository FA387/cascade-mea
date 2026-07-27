# CASCADE v0.8.0

**Criticality Avalanche Spike Cross-platform Analysis Detection Engine**

<img width="750" height="750" alt="CASCADE_Logo" src="https://github.com/user-attachments/assets/7a029bf6-93b7-4185-8859-888672b4bb65" />

CASCADE is a modular Python pipeline for analysing multi-electrode array (MEA) recordings from multiple manufacturers. It reads manufacturer-specific spike files, converts them into a common internal format, and runs a standardised sequence of analyses: firing-rate characterisation, network burst detection, neuronal avalanche segmentation, criticality metrics, functional connectivity, and spatial entropy. All figures are saved as PNG (300 dpi) and all numerical results as TSV files. The pipeline is designed to produce reproducible results suitable for peer-reviewed publication.

---

## Tools

| Command | Description |
|---|---|
| `cascade-info` | Inspect recording: duration, channels, active electrodes, statistics |
| `cascade` | Full analysis pipeline: raster, rates, bursts, avalanches, criticality, connectivity, entropy |
| `cascade-ltpltd` | LTP/LTD classification by comparing firing rates between recordings |

---

## Supported manufacturers

| Manufacturer | File format | Notes |
|---|---|---|
| 3Brain BrainWave | `.bxr` | Reads SpikeTimes, SpikeChIdxs, SpikeUnits, SpikeForms |
| Maxwell Biosystems | `.h5` | Reads data_store/stream/spikes (frameno, channel, amplitude) |
| Multi Channel Systems | `.h5` | MCS HDF5 v3; reads SegmentStream spike timestamps |
| Cortical Labs | `.h5` | PyTables HDF5; reads /spikes table (timestamp, channel) |
| Axion BioSystems | `.csv` | Reads spike timestamps from AxIS spike-list CSV exports |

---

## Installation

```bash
git clone https://github.com/FA387/cascade-mea.git
cd cascade-mea
conda env create -f environment.yml
conda activate cascade-mea
pip install .
```

### Verify installation

```bash
cascade --test
```

**Expected output:**

```
================================
  CASCADE v0.8.0 — Self-test
================================
  1. Synthetic recording generation... PASS
  2. Spike filtering (active-electrode)... PASS
  3. ISI metrics computation... PASS
  4. Spike binning and rates... PASS
  5. Network burst detection... PASS
  6. Active area analysis... PASS
  7. Avalanche detection... PASS
  8. Criticality metrics (tau, alpha, DCC, BR, SCE)... PASS
  9. Functional connectivity... PASS
 10. Spatial entropy... PASS

All 10 tests passed ✓
```

If any test fails, check your environment (missing dependencies, incorrect Python version, etc.).

---

## Quick start

### Manufacturer spike loaders

CASCADE supports five MEA manufacturers. All tools require `--manufacturer` to identify the file format and load spikes correctly. Choose the value that matches your recording:

| Manufacturer | Flag | File format |
|---|---|---|
| 3Brain BrainWave | `3brain` | `.bxr` |
| Maxwell Biosystems | `maxwell` | `.h5` |
| Cortical Labs | `cl` | `.h5` |
| Multi Channel Systems | `mcs` | `.h5` |
| Axion BioSystems | `axion` | `.csv` |

### cascade-info

Inspect a recording before running the full analysis. Prints summary statistics: duration, channel count, active electrodes, and per-well metadata.

```bash
# Full recording
cascade-info --manufacturer [manufacturer] --input [input_recording]

# Specific time window (0–60 seconds)
cascade-info --manufacturer [manufacturer] --input [input_recording] --start 0 --dur 60

# Active-electrode filter (channels with ≥5 spikes/min)
cascade-info --manufacturer [manufacturer] --input [input_recording] --ae 5

# Batch audit: scan folder, save summary (duration, channels, active count per file)
cascade-info --folder [input_directory] --manufacturer [manufacturer] --ae 5 --output [output_summary.tsv]
```

### cascade

Full analysis pipeline: generates raster plots, firing rates, burst detection, avalanche analysis, criticality metrics (tau, alpha, DCC, branching ratio, shape collapse error), functional connectivity, and spatial entropy.

```bash
# Basic: all defaults
cascade --manufacturer [manufacturer] --input [input_recording] --output [output_dir]

# Habibollahi 2022 method (percentile threshold + active-electrode filter)
cascade --manufacturer [manufacturer] --input [input_recording] --ae 5 --perc 0.25 --output [output_dir]

# Analyze specific time window (0–120 seconds)
cascade --manufacturer [manufacturer] --input [input_recording] --start 0 --dur 120 --output [output_dir]

# Single well from multi-well Axion recording
cascade --manufacturer axion --input [input_recording] --well [well_id] --output [output_dir]

# Custom burst detection + shape collapse parameters
cascade --manufacturer [manufacturer] --input [input_recording] --burst_min_spikes 20 --sc_Dmin 4 --sc_Dmax 20 --output [output_dir]

# Compare two time windows (baseline vs. post-treatment)
cascade --manufacturer [manufacturer] --input [input_recording] --epochA 0,60 --epochB 120,60 --output [output_dir]

# Save figures as SVG in addition to PNG
cascade --manufacturer [manufacturer] --input [input_recording] --svg --output [output_dir]

# Test mode (synthetic data, no real file needed)
cascade --test
```

### cascade-ltpltd

LTP/LTD (long-term potentiation/depression) classification by comparing firing rates between two recordings. Outputs per-electrode calls (LTP/LTD/NoChange) and summary figures.

```bash
# Single condition
cascade-ltpltd --baseline [baseline_rates.tsv] \
                --measure [stim_rates.tsv] \
                --output [output_dir]

# Custom threshold (default 20%)
cascade-ltpltd --baseline [baseline_rates.tsv] --measure [stim_rates.tsv] --threshold 30 --output [output_dir]

# Multiple conditions
cascade-ltpltd --baseline [baseline_rates.tsv] \
                --measure [stim1_rates.tsv] [stim2_rates.tsv] [stim3_rates.tsv] \
                --output [output_dir]
```

---

## Output files

### cascade output

```
mea_out/
├── raster_max_<well>.png              # spike raster — all electrode slots
├── raster_active_<well>.png           # spike raster — active electrodes only
├── rate_heatmap_max_<well>.png        # firing rate heatmap — all electrode slots
├── rate_heatmap_active_<well>.png     # firing rate heatmap — active electrodes only
├── population_rate_<well>.png         # population mean firing rate over time
├── bursts_timeline_<well>.png         # network burst events over time
├── avalanches_<well>_size.png         # avalanche size distribution
├── avalanches_<well>_duration.png     # avalanche duration distribution
├── loglog_size_<well>.png             # log-log size distribution with power-law fit
├── loglog_duration_<well>.png         # log-log duration distribution with power-law fit
├── shape_collapse_<well>.png          # avalanche shape collapse
├── corr_heatmap_max_<well>.png        # correlation matrix — all electrode slots
├── corr_heatmap_active_<well>.png     # correlation matrix — active electrodes only
├── degree_hist_<well>.png             # node degree distribution
├── entropy_curve_<well>.png           # spatial entropy over time
├── active_area_curve_<well>.png       # active electrode fraction over time
├── epoch_compare_<well>.png           # (optional) epoch A vs B comparison
│
├── rates_<well>.tsv                   # per-channel firing rates (binned)
├── bursts_<well>.tsv                  # network burst event table
├── avalanches_<well>.tsv              # avalanche catalogue
├── corr_<well>.tsv                    # raw Pearson correlation matrix
├── isi_metrics_<well>.tsv             # per-channel ISI statistics
├── active_area_<well>.tsv             # active electrode fraction per window
├── entropy_<well>.tsv                 # spatial entropy per window
├── degrees_<well>.tsv                 # node degree per channel
├── branching_points_<well>.tsv        # N(t) vs N(t+1) pairs for branching ratio fit
├── criticality_report_<well>.txt      # plain-text criticality metrics summary
└── cascade_run.log                    # command, parameters, and run status
```

#### Figures

**`raster_max_<well>.png`**
Spike raster plot with the Y-axis spanning all electrode slots on the array, including silent ones. Each dot is a single spike event. Useful for assessing spatial coverage and identifying which regions of the array are active.

**`raster_active_<well>.png`**
Spike raster restricted to electrodes that passed the active-electrode filter, with the Y-axis remapped to 0…N. Provides a cleaner view of firing patterns when the array has many silent slots.

**`rate_heatmap_max_<well>.png`**
Firing rate heatmap with all electrode slots on the Y-axis and time on the X-axis. Colour encodes spikes per second in each bin. Shows how activity is distributed spatially across the full array over time.

**`rate_heatmap_active_<well>.png`**
Same heatmap restricted to active electrodes, remapped to contiguous row indices. Produces a denser, more readable plot when the proportion of active electrodes is low.

**`population_rate_<well>.png`**
Mean firing rate averaged across all active channels, plotted over time at the `--bin` resolution. Provides a single-trace overview of global activity dynamics and makes bursting epochs and silent periods immediately visible.

**`bursts_timeline_<well>.png`**
Timeline showing all detected network burst events as shaded intervals across the recording. Start time, duration, and inter-burst interval are visible at a glance. An empty plot indicates no bursts were detected above the minimum spike threshold.

**`avalanches_<well>_size.png` and `avalanches_<well>_duration.png`**
Histograms of avalanche sizes (total spikes per avalanche) and durations (number of Δt bins per avalanche). Show the range and frequency of avalanche events after size and duration thresholding. Reviewed before fitting power laws to confirm there are sufficient events.

**`loglog_size_<well>.png` and `loglog_duration_<well>.png`**
Log–log plots of the avalanche size and duration distributions with the fitted power-law line overlaid. The fitted exponents τ (size) and α (duration) and their respective xmin values are annotated. A straight line on a log–log scale is the primary visual signature of scale-free, critical dynamics.

**`shape_collapse_<well>.png`**
Avalanche profiles of different durations are rescaled by their predicted scaling exponent and overlaid. A tight collapse onto a single universal curve is a strong indicator of criticality beyond what the exponents alone can confirm. The shape collapse error (SCE) quantifies how well the profiles align.

**`corr_heatmap_max_<well>.png`**
Pearson correlation matrix across all electrode slots, with pairs above `--conn_thr` considered functionally connected. The full-grid layout makes it possible to see where connected and silent electrodes fall relative to each other on the array.

**`corr_heatmap_active_<well>.png`**
Correlation matrix restricted to active electrodes. Provides a cleaner view of the functional connectivity structure without silent electrode slots diluting the colour scale.

**`degree_hist_<well>.png`**
Histogram of node degree — the number of functional connections per electrode above the correlation threshold. Characterises whether connectivity is sparse or dense and whether hub electrodes with many connections are present.

**`entropy_curve_<well>.png`**
Normalised spatial entropy of the spike distribution computed in sliding windows of `--entropy_win` seconds. High entropy indicates spikes are distributed evenly across active electrodes; low entropy indicates spatially localised activity. Tracks changes in spatial organisation over the course of the recording.

**`active_area_curve_<well>.png`**
Fraction of active electrodes that fired in each time window. Reveals how broadly or focally the network recruits its active population over time and complements the entropy curve in characterising spatial dynamics.

**`epoch_compare_<well>.png`**
Bar chart comparing population firing rate and active area fraction between epoch A and epoch B. Only produced when both `--epochA` and `--epochB` are provided. Provides a concise visual summary of the effect of an intervention or stimulation.

#### Data files

**`rates_<well>.tsv`**
Per-channel firing rate in Hz, binned at the `--bin` interval. Each row is one channel–time-bin combination. This file is the primary input for `cascade-ltpltd` and for any downstream R-based analysis.

**`bursts_<well>.tsv`**
Table of detected network burst events with start time, end time, duration, and spike count per burst. An empty table indicates that no bursts exceeded the minimum spike threshold set by `--burst_min_spikes`.

**`avalanches_<well>.tsv`**
Catalogue of all avalanche events that passed the size and duration thresholds. Each row records the avalanche start time, total spike count (size), and duration in Δt bins. This is the dataset used for all criticality computations.

**`corr_<well>.tsv`**
Raw Pearson correlation matrix saved as a tab-separated square array. Rows and columns correspond to active channels in the order listed in `degrees_<well>.tsv`. Included for custom connectivity analyses beyond what CASCADE computes.

**`isi_metrics_<well>.tsv`**
Per-channel inter-spike interval (ISI) statistics: mean ISI, median ISI, coefficient of variation, mean firing rate, and refractory period violation fraction. Useful for characterising single-channel spiking regularity and for identifying channels with unusually high violation rates before applying noise filters.

**`active_area_<well>.tsv`**
Numerical counterpart of `active_area_curve_<well>.png`. Each row is one time window with start time, end time, and the fraction of active electrodes that fired. Used for quantitative comparisons of spatial recruitment across recordings or conditions.

**`entropy_<well>.tsv`**
Spatial entropy values per time window. Columns record window start, window end, and normalised entropy. Numerical counterpart of `entropy_curve_<well>.png`, included for statistical comparisons across recordings.

**`degrees_<well>.tsv`**
Node degree for each active channel: the number of other channels with which it shares a correlation above `--conn_thr`. Includes channel identifiers for cross-referencing with the raster and heatmap outputs.

**`branching_points_<well>.tsv`**
Population spike count pairs N(t) and N(t+1) for each Δt bin, used as input for the branching ratio linear regression. Included so the BR fit can be reproduced, visualised, or checked independently of the main pipeline.

**`criticality_report_<well>.txt`**
Plain-text summary of all criticality metrics: τ, α, β_pred, β_fit, DCC, branching ratio (BR), BR intercept, shape collapse error (SCE), and the xmin values used in the power-law fits. This is the primary result file for reporting and manual inspection.

**`cascade_run.log`**
Records the exact command, all parameters, timestamp, Python version, platform, and success or failure status for each run. Supports reproducibility by preserving the full analysis context alongside the outputs.

---

### cascade-info output

Results printed to terminal. With `--output [file.tsv]` in folder mode, saves summary table.

### cascade-ltpltd output

```
ltpltd_out/
├── ltpltd_[label].tsv               # per-electrode LTP/LTD call
├── ltpltd_[label]_bar.png           # bar chart (% change)
├── ltpltd_[label]_pie.png           # pie chart (LTP/LTD/NoChange proportions)
└── ltpltd_all_conditions.tsv        # (optional) all conditions combined
```

---

## All flags

### cascade-info

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `maxwell`, `cl1`, `mcs`, or `axion` |
| `--input` | — | Path to `.bxr`, `.h5`, or `.csv` file |
| `--folder` | — | Folder to scan (all `.bxr` and `.h5` files) |
| `--well` | first | Well or stream name |
| `--start` | `0` s | Window start (seconds) |
| `--dur` | full | Window duration (seconds) |
| `--ae` | off | Active-electrode filter: ≥N spikes/min |
| `--output` | — | Save summary to this TSV path (folder mode only) |
| `--sr` | — | Sampling rate Hz (for plain CSV input only) |

### cascade

**Input/Output**

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `maxwell`, `cl1`, `mcs`, or `axion` |
| `--input` | auto-search | Path to `.bxr`, `.h5`, or `.csv` file |
| `--well` | first | Well or stream name |
| `--output` | `mea_out` | Output folder |
| `--svg` / `--jpg` | off | Save figures as SVG/JPG in addition to PNG |
| `--verbose` | off | Print progress to terminal |

**Recording window**

| Flag | Default | Description |
|---|---|---|
| `--start` | `0` | Analysis window start (seconds) |
| `--dur` | full | Analysis window duration (seconds) |

**Spike filtering**

| Flag | Default | Description |
|---|---|---|
| `--ae` | off | Active-electrode: keep channels with ≥N spikes/min |
| `--max_rate` | off | Noise filter: remove if rate >N spikes/min (try `1000`) |
| `--max_rpv` | off | Noise filter: remove if >N% of ISIs <2 ms (try `1.0`) |
| `--bin` | `1.0` | Bin size for rate computation (seconds) |

**Avalanche detection**

| Flag | Default | Description |
|---|---|---|
| `--ava_size` | `2` | Minimum avalanche size (spikes). Single spikes are not avalanches. |
| `--ava_dur` | `1` | Minimum avalanche duration (bins). Default=1 allows synchronous bursts. |
| `--aval_dt` | inferred | Avalanche bin width (seconds); auto = mean ISI clipped to 1–50 ms |
| `--perc` | off | Population activity percentile threshold. Use `0.25` for Habibollahi 2022 method |

**Criticality metrics**

| Flag | Default | Description |
|---|---|---|
| `--burst_min` | `12` | Upper xmin search limit (size exponent fit) |
| `--t_min` | `4` | Upper xmin search limit (duration exponent fit) |
| `--pl_xmin_S` | auto | Manual xmin override (size fit) |
| `--pl_xmin_D` | auto | Manual xmin override (duration fit) |
| `--skip_rare_sizes` | off | Exclude singletons from log-log fit (does not affect MLE) |
| `--skip_rare_durations` | off | Exclude singletons from log-log fit (does not affect MLE) |

**Shape collapse**

| Flag | Default | Description |
|---|---|---|
| `--sc_Dmin` | `4` | Minimum duration (bins) in shape collapse |
| `--sc_Dmax` | `20` | Maximum duration (bins) in shape collapse |

**Burst detection**

| Flag | Default | Description |
|---|---|---|
| `--burst_min_spikes` | `5` | Minimum spikes per burst event |
| `--burst_max_gap` | `0.1` | Maximum inter-spike gap within burst (seconds) |

**Connectivity**

| Flag | Default | Description |
|---|---|---|
| `--conn_thr` | `0.3` | Pearson \|r\| threshold for functional connection |

**Entropy**

| Flag | Default | Description |
|---|---|---|
| `--entropy_win` | `10` | Window size for entropy and active-area curves (seconds) |

**Epoch comparison**

| Flag | Default | Description |
|---|---|---|
| `--epochA` | off | Epoch A as `start,duration` (e.g. `0,60` = 0–60 s). Both A and B required. |
| `--epochB` | off | Epoch B as `start,duration` (e.g. `120,60` = 120–180 s). Both A and B required. |

### cascade-ltpltd

| Flag | Default | Description |
|---|---|---|
| `--baseline` | required | Path to `rates_*.tsv` from baseline CASCADE run |
| `--measure` | required | Path(s) to one or more `rates_*.tsv` files to compare |
| `--threshold` | `20` | LTP/LTD threshold (% change from baseline) |
| `--output` | `ltpltd_out` | Output folder |
| `--svg` / `--jpg` | off | Save figures as SVG/JPG |
| `--verbose` | off | Print per-electrode table to terminal |

---

## Key metrics

| Metric | Definition | Reference |
|---|---|---|
| **τ** (tau) | Avalanche size exponent | Power-law fit; tplfit KS-minimisation MLE |
| **α** (alpha) | Avalanche duration exponent | Power-law fit; tplfit KS-minimisation MLE |
| **β_pred** | Predicted scaling exponent | (α − 1) / (τ − 1) |
| **β_fit** | Fitted scaling exponent | Slope of log(mean size) vs. log(duration), [t_min, t_min+60] |
| **DCC** | Deviation from Criticality Coefficient | \|β_pred − β_fit\| |
| **BR** | Branching ratio | E[N(t+1)\|N(t)] slope (linear regression) |
| **SCE** | Shape Collapse Error | \|β_opt − β_pred\| |

---

## Noise filtering

MEA recordings sometimes contain electrodes picking up electrical noise rather than real spikes. CASCADE provides two independent filters:

```bash
# Remove high-rate channels (>1000 spikes/min ≈ 17 Hz; typical noise signature)
cascade --manufacturer [manufacturer] --input [input_recording] --max_rate 1000 --output [output_dir]

# Remove channels with >1% refractory period violations (<2 ms inter-spike intervals; biological impossibility)
cascade --manufacturer [manufacturer] --input [input_recording] --max_rpv 1.0 --output [output_dir]

# Both filters together (recommended for noisy recordings)
cascade --manufacturer [manufacturer] --input [input_recording] --ae 5 --max_rate 1000 --max_rpv 1.0 --output [output_dir]
```

Use `--verbose` to see which channels were removed and why.

---

## Project structure

```
cascade/
├── run_analysis.py          cascade command
├── scanner.py               cascade-info command
├── ltpltd.py                cascade-ltpltd command
├── environment.yml
│
├── loaders/                 manufacturer-specific I/O
│   ├── base_loader.py
│   ├── brainwave_loader.py  (3Brain .bxr)
│   ├── maxwell_loader.py    (Maxwell .h5)
│   ├── cl1_loader.py        (Cortical Labs .h5)
│   ├── mcs_loader.py        (Multi Channel Systems .h5)
│   └── axion_loader.py      (Axion .csv)
│
├── core/                    analysis algorithms
│   ├── spikes.py            filtering, binning, ISI metrics
│   ├── bursts.py            burst detection
│   ├── avalanche.py         avalanche segmentation
│   ├── criticality.py       power-law fitting, metrics
│   ├── connectivity.py      functional connectivity
│   ├── entropy.py           spatial entropy
│   └── plots.py             figure generation
│
├── utils/
│   ├── common.py            MEARecording dataclass, helpers
│   ├── io.py                CSV/TSV writers
│   └── inspect_file.py      HDF5 structure inspection
│
└── tests/                   unit tests (in progress)
```

Loaders convert manufacturer files → `MEARecording` object. Core modules only operate on that object. To add a manufacturer: write one new loader file and register it in `run_analysis.py`.

---

## Adding a new manufacturer

1. Create `loaders/your_manufacturer_loader.py`
2. Subclass `BaseLoader` from `loaders/base_loader.py`
3. Implement `load(path, well=None, start_s=0.0, duration_s=None) → MEARecording`
4. Register in `run_analysis.py` in the `get_loader()` function

Minimal template:

```python
from .base_loader import BaseLoader
from ..utils.common import MEARecording

class YourLoader(BaseLoader):
    @staticmethod
    def manufacturer_name() -> str:
        return "yourmanufacturer"

    def load(self, path, well=None, start_s=0.0, duration_s=None) -> MEARecording:
        # Read the file and return MEARecording(
        #   sampling_rate, n_channels, channel_ids,
        #   spike_times={ch_id: np.array(times_in_seconds)},
        #   duration, metadata={...}
        # )
        ...
```

Run `python -m utils.inspect_file dump [yourfile]` to explore file structure.

---

## References

- Habibollahi, F., Kagan, B. J., Burkitt, A. N., & French, C. (2023). Critical dynamics arise during structured information presentation within embodied in vitro neuronal networks. *Nature Communications*, **14**(1), 5287. https://doi.org/10.1038/s41467-023-41020-3

- Beggs, J. M., & Plenz, D. (2003). Neuronal Avalanches in Neocortical Circuits. *The Journal of Neuroscience*, **23**(35), 11167–11177. https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003

- Wilting, J., & Priesemann, V. (2018). Inferring collective dynamical states from widely unobserved systems. *Nature Communications*, **9**(1), 2325. https://doi.org/10.1038/s41467-018-04725-4

- Friedman, N., Ito, S., Brinkman, B. A. W., Shimono, M., DeVille, R. E. L., Dahmen, K. A., Beggs, J. M., & Butler, T. C. (2012). Universal Critical Dynamics in High Resolution Neuronal Avalanche Data. *Physical Review Letters*, **108**(20), 208102. https://doi.org/10.1103/PhysRevLett.108.208102

- Poil, S.-S., Hardstone, R., Mansvelder, H. D., & Linkenkaer-Hansen, K. (2012). Critical-State Dynamics of Avalanches and Oscillations Jointly Emerge from Balanced Excitation/Inhibition in Neuronal Networks. *Journal of Neuroscience*, **32**(29), 9817–9823. https://doi.org/10.1523/JNEUROSCI.5990-11.2012

---

## Funding

This work was supported by:

- **Korea Institute of Toxicology (KIT)** (Project No. KK-2607)
- **Korea Ministry of Food and Drug Safety** (RS-2024-00331678)

---

## Contributing

Contributions welcome, especially new manufacturer loaders. Open an issue or PR on GitHub.

---

## Licence

Apache 2.0
