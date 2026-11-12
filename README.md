# CASCADE v0.9.5

**Criticality Avalanche Spike Cross-platform Analysis Detection Engine**

<div align="center">
<img width="250" height="250" alt="CASCADE_Logo" src="https://github.com/user-attachments/assets/7a029bf6-93b7-4185-8859-888672b4bb65" />
</div>

CASCADE is a modular Python pipeline for analysing multi-electrode array (MEA) recordings from multiple manufacturers. It reads manufacturer-specific spike files, converts them into a common internal format, and runs a standardised sequence of analyses: firing-rate characterisation, network burst detection, neuronal avalanche segmentation, criticality metrics, functional connectivity, and spatial entropy. All figures are saved as PNG (300 dpi) and all numerical results as TSV files. The pipeline is designed to produce reproducible results suitable for peer-reviewed publication.

Published in *Bioinformatics* (2026), [doi:10.1093/bioinformatics/btag570](https://doi.org/10.1093/bioinformatics/btag570) — see [Citing CASCADE](#citing-cascade). The article describes v0.8.0; see [Differences from the published article](#differences-from-the-published-article) for what has changed since.

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
CASCADE v0.9.5 — Built-in Self-Test
===================================================================
  Synthetic recording: 20 ch · 120 s · Poisson + bursts · seed 42
===================================================================
  Synthetic MEA recording ................................ PASS
  Spike filtering (active-electrode filter) .............. PASS
  Spike filtering (amplitude filter) ..................... PASS
  ISI metrics ............................................ PASS
  Spike binning .......................................... PASS
  Network burst detection ................................ PASS
  Active area ............................................ PASS
  Avalanche detection .................................... PASS
  Synchrony metrics (xcorr, Kreuz, Paiva) ................ PASS
  Power-law normalisation (known exponent) ............... PASS
  Criticality metrics (tau, alpha, DCC, BR) .............. PASS
  Branching ratio (MR estimator, subsampled) ............. PASS
  Functional connectivity ................................ PASS
  Spatial entropy ........................................ PASS
===================================================================
  All 14 tests passed. CASCADE is installed correctly.
```

Check the installed version with `cascade --version`; `cascade --citation` prints how to cite it.

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

# 25th-percentile avalanche threshold + active-electrode filter
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

# Reproduce the published v0.8.0 numbers (see "Reproducing the published v0.8.0 results")
cascade --manufacturer [manufacturer] --input [input_recording] --simple_PL --simple_BR --output [output_dir]
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
├── wave_templates_<well>.png          # (3Brain with stored waveforms) mean spike waveform per channel
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
├── synchrony_<well>.tsv               # synchrony metrics (Axion metric set)
├── crosscorrelogram_<well>.png        # pooled cross-correlogram
├── crosscorrelogram_norm_<well>.png   # auto-correlation-normalised correlogram
├── isi_metrics_<well>.tsv             # per-channel ISI statistics
├── active_area_<well>.tsv             # active electrode fraction per window
├── entropy_<well>.tsv                 # spatial entropy per window
├── degrees_<well>.tsv                 # node degree per channel
├── branching_mr_<well>.tsv            # MR estimator: r_k per lag with both fits
├── branching_mr_<well>.png            # r_k vs lag with the fitted exponentials
├── branching_points_<well>.tsv        # mean N(t+1) per N(t), for the --simple_BR fit
├── criticality_report_<well>.txt      # plain-text criticality metrics summary
└── cascade_run.log                    # command, CASCADE version, and run status
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

**`wave_templates_<well>.png`**
Mean spike waveform per channel. Produced only for 3Brain `.bxr` files that store spike waveforms (`SpikeForms`); useful for checking that detected spikes look like action potentials rather than artefacts.

**`epoch_compare_<well>.png`**
Bar chart comparing population firing rate and active area fraction between epoch A and epoch B. Only produced when both `--epochA` and `--epochB` are provided. Provides a concise visual summary of the effect of an intervention or stimulation.

#### Data files

**`rates_<well>.tsv`**
Per-channel firing rate in Hz, binned at the `--bin` interval. Each row is one channel–time-bin combination. This file is the primary input for `cascade-ltpltd` and for any downstream R-based analysis.

**`bursts_<well>.tsv`**
Table of detected network burst events with start time, end time, duration, and spike count per burst. An empty table indicates that no bursts exceeded the minimum spike threshold set by `--burst_min_spikes`.

**`avalanches_<well>.tsv`**
Catalogue of all avalanche events that passed the size and duration thresholds. Each row records the avalanche start time, total spike count (size), and duration in Δt bins. This is the dataset used for all criticality computations.

**`synchrony_<well>.tsv`**
Two-column metric/value table of the synchrony metrics, matching the Axion Neural Metric Tool set: `area_xcorr`, `area_xcorr_norm`, `fwhh_xcorr`, `fwhh_xcorr_norm` and, when requested, `kreuz_metric` and `synchrony_index`. Also records the parameters used (`sync_window_s`, `xcorr_bin_s`, `sync_index_sigma_s`), the pair counts, and `kreuz_subsampled` / `index_subsampled` flags so a subsampled result is never mistaken for a full one. Durations are in seconds. See the Synchrony metrics section for interpretation.

**`crosscorrelogram_<well>.png`** and **`crosscorrelogram_norm_<well>.png`**
Well-wide pooled inter-electrode cross-correlogram, standard and auto-correlation-normalised. Grey vertical lines mark the synchrony window and the shaded region is the area reported as `area_xcorr` / `area_xcorr_norm`. A tall narrow peak means high synchrony; a flat correlogram means none, in which case the FWHH value describes noise rather than network structure.

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

**`branching_mr_<well>.tsv`** / **`branching_mr_<well>.png`**
The multi-lag (MR) branching ratio of Wilting & Priesemann (2018), the default since v0.9.4: the regression slope r_k of a(t+k) on a(t) for every lag k, with the exponential fit r_k = b·m^k (m = BR) and an exponential+offset fit used as a drift check. Inspect the plot: r_k should decay exponentially to 0.

**`branching_points_<well>.tsv`**
Population spike count pairs N(t) and N(t+1) for each Δt bin, used as input for the simple single-lag regression reported with `--simple_BR`. Included so the BR fit can be reproduced, visualised, or checked independently of the main pipeline.

**`criticality_report_<well>.txt`**
Plain-text summary of all criticality metrics: τ, α and their xmin values, β_pred, β_fit, DCC, the branching ratio (BR, with the MR fit's b, autocorrelation time τ, lag range and its three quality checks, plus the simple single-lag value for comparison), and the shape collapse error (SCE). Any `--simple_PL`, `--simple_BR`, `--pl_xmin_S/D` or `--pl_exclude` choice is stated in the report. This is the primary result file for reporting and manual inspection.

**`cascade_run.log`**
Records the exact command (and therefore every parameter), the CASCADE version, timestamp, Python version, platform, working directory, success or failure status, and the terminal output for each run. Report the CASCADE version with your results: criticality metrics differ between versions (see `CHANGELOG.md`).

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
| `--manufacturer` | required | `3brain`, `maxwell`, `cl`, `mcs`, or `axion` |
| `--input` | — | Path to `.bxr`, `.h5`, or `.csv` file |
| `--folder` | — | Folder to scan (all `.bxr` and `.h5` files) |
| `--well` | first | Well or stream name |
| `--start` | `0` s | Window start (seconds) |
| `--dur` | full | Window duration (seconds) |
| `--ae` | off | Active-electrode filter: ≥N spikes/min |
| `--output` | — | Save summary to this TSV path (folder mode only) |
| `--sr` | — | Sampling rate Hz (for plain CSV input only) |
| `--version` / `--citation` | — | Print the version / how to cite CASCADE, and exit |

### cascade

**General**

| Flag | Default | Description |
|---|---|---|
| `--test` | — | Run the 14 built-in self-tests on synthetic data (no recording needed) |
| `--version` | — | Print the CASCADE version and exit |
| `--citation` | — | Print how to cite CASCADE (text and BibTeX) and exit |

**Input/Output**

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `maxwell`, `cl`, `mcs`, or `axion` |
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
| `--max_amplitude` | off | Noise filter: remove if any spike exceeds N µV peak-to-peak (try `200`). Supported by CL1, 3Brain, Maxwell, MCS, and Axion (all in µV). |
| `--bin` | `1.0` | Bin size for rate computation (seconds) |

**Avalanche detection**

| Flag | Default | Description |
|---|---|---|
| `--ava_size` | `2` | Minimum avalanche size (spikes). Single spikes are not avalanches. |
| `--ava_dur` | `1` | Minimum avalanche duration (bins). Default=1 allows synchronous bursts. |
| `--aval_dt` | inferred | Avalanche bin width (seconds); auto = mean ISI clipped to 1–50 ms |
| `--perc` | off | Population activity percentile threshold. `0.25` is the default of the Habibollahi et al. (2023) reference code (their paper reports 40% of the median) |

**Criticality metrics**

| Flag | Default | Description |
|---|---|---|
| `--burst_min` | `12` | Upper xmin search limit (size exponent fit) |
| `--t_min` | `4` | Upper xmin search limit (duration exponent fit) |
| `--simple_PL` | off | Use the power-law normalisation of v0.8.0–v0.9.0 (biased below ~1,000 avalanches). Use **only** to regenerate previously published numbers; add `--simple_BR` as well to reproduce v0.8.0 exactly |
| `--pl_xmin_S` | auto | Fix the size-fit xmin instead of searching. τ, DCC, SCE and the plot use it; the automatic fit is still written to the report |
| `--pl_xmin_D` | auto | Fix the duration-fit xmin instead of searching (α, β_fit, DCC) |
| `--pl_exclude` | off | Also trim the upper tail before fitting, with the `EXCLUDE` step of the Habibollahi et al. (2023) Python reference code (the MATLAB original differs). For comparison with that pipeline only: add `--simple_PL` to match its code exactly. The untrimmed fit is still reported — see `CHANGELOG.md` |
| `--skip_rare_sizes` | off | Exclude singletons from log-log fit (does not affect MLE) |
| `--skip_rare_durations` | off | Exclude singletons from log-log fit (does not affect MLE) |

**Shape collapse**

| Flag | Default | Description |
|---|---|---|
| `--sc_Dmin` | `4` | Minimum duration (bins) in shape collapse |
| `--sc_Dmax` | `20` | Maximum duration (bins) in shape collapse |

**Branching ratio**

| Flag | Default | Description |
|---|---|---|
| `--simple_BR` | off | Report the simple single-lag branching ratio of v0.8.0–v0.9.3 instead of the MR estimator. Biased toward 0 under subsampling; use **only** to reproduce previously published numbers; add `--simple_PL` as well to reproduce v0.8.0 exactly |
| `--br_tmax` | `2000` | Longest lag (ms) for the MR estimator. Should span several autocorrelation times; the report warns when it does not |

**Burst detection**

| Flag | Default | Description |
|---|---|---|
| `--burst_min_spikes` | `5` | Minimum spikes per burst event |
| `--burst_max_gap` | `0.1` | Maximum inter-spike gap within burst (seconds) |

**Connectivity**

| Flag | Default | Description |
|---|---|---|
| `--conn_thr` | `0.3` | Pearson \|r\| threshold for functional connection |

**Synchrony**

| Flag | Default | Description |
|---|---|---|
| `--sync_metrics` | `core` | `none` / `core` / `xcorr` / `kreuz` / `index` / `all`. `core` = normalised correlogram area + Synchrony Index. The correlogram is cheap (FFT); `kreuz` and `index` are O(N²) in channels |
| `--sync_window` | `20.0` | Synchrony window in **ms** — the window around zero lag used for the area metrics. `5` for millisecond timescale, `100` for slow |
| `--sync_bin` | `1.0` | Cross-correlogram bin width in **ms**. Absolute, never inferred from sampling rate |
| `--sync_max_lag` | `500.0` | Maximum correlogram lag in **ms**. Must be wide enough to contain the peak or FWHH cannot be measured |
| `--sync_sigma` | = window | Gaussian σ in **ms** for the Synchrony Index |
| `--sync_max_channels` | `64` | Channel cap for the O(N²) measures; `0` = no cap |

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
| `--version` / `--citation` | — | Print the version / how to cite CASCADE, and exit |

All three tools name their output folder or file `--output`; `cascade` and `cascade-info` read one recording with `--input`, while `cascade-ltpltd` takes two kinds of input, `--baseline` and `--measure`.

---

## Key metrics

| Metric | Definition | Reference |
|---|---|---|
| **τ** (tau) | Avalanche size exponent | Power-law fit; tplfit KS-minimisation MLE |
| **α** (alpha) | Avalanche duration exponent | Power-law fit; tplfit KS-minimisation MLE |
| **β_pred** | Predicted scaling exponent | (α − 1) / (τ − 1) |
| **β_fit** | Fitted scaling exponent | Slope of log(mean size) vs. log(duration) for xmin_D < D < xmin_D + 60 (xmin_D = fitted duration cut-off) |
| **DCC** | Deviation from Criticality Coefficient | \|β_pred − β_fit\| |
| **BR** | Branching ratio m | Multi-lag MR estimator: r_k = b·m^k (Wilting & Priesemann 2018). `--simple_BR`: single-lag E[N(t+1)\|N(t)] slope, biased toward 0 under subsampling |
| **SCE** | Shape Collapse Error | \|β_opt − β_pred\| |

---

## Noise filtering

MEA recordings sometimes contain electrodes picking up electrical noise rather than real spikes. CASCADE provides three independent filters:

```bash
# Remove high-rate channels (>1000 spikes/min ≈ 17 Hz; typical noise signature)
cascade --manufacturer [manufacturer] --input [input_recording] --max_rate 1000 --output [output_dir]

# Remove channels with >1% refractory period violations (<2 ms inter-spike intervals; biological impossibility)
cascade --manufacturer [manufacturer] --input [input_recording] --max_rpv 1.0 --output [output_dir]

# Remove channels where any spike exceeds 200 µV peak-to-peak (artefact rejection)
cascade --manufacturer [manufacturer] --input [input_recording] --max_amplitude 200 --output [output_dir]

# All three filters together (recommended for noisy recordings)
cascade --manufacturer [manufacturer] --input [input_recording] --ae 5 --max_rate 1000 --max_rpv 1.0 --max_amplitude 200 --output [output_dir]
```

Use `--verbose` to see which channels were removed and why.

Amplitudes are reported in µV for every manufacturer. Note that Maxwell stores a single-peak magnitude rather than a true peak-to-peak value, so a threshold is not strictly equivalent across devices — calibrate per device with `cascade-info` before settling on one.

---

## Synchrony metrics

CASCADE computes the synchrony metric set of the Axion Neural Metric Tool, applied
uniformly to **every** supported manufacturer. These are spike-timing measures and
are distinct from the Pearson connectivity above, which correlates firing rates in
1 s bins and therefore cannot see millisecond structure.

### The two default metrics

`--sync_metrics core` (the default) reports just two:

| Metric | Output key | Source |
|---|---|---|
| Area under normalised cross-correlation | `area_xcorr_norm` | Axion NMT §4.4 |
| Synchrony Index | `synchrony_index` | Paiva et al. 2010 |

These are the two that **resolve on every manufacturer** and **actually separate
recordings**. Measured across the five published sample recordings:

| metric | spread (max/min) | CV | resolves on | verdict |
|---|---:|---:|---|---|
| `area_xcorr_norm` | 26× | 0.78 | 5/5 | **core** |
| `synchrony_index` | 11× | 0.69 | 5/5 | **core** |
| `area_xcorr` (raw) | 18× | 0.75 | 5/5 | same information, confounded by firing rate |
| `fwhh_xcorr_norm` | 552× | 1.47 | **2/5** | large spread is mostly the resolution artefact |
| `kreuz_metric` | **1.2×** | **0.09** | 5/5 | near-constant across very different cultures |

Fewer, better-behaved numbers beat six that need caveats.

### The rest, on request

| `--sync_metrics` | Adds |
|---|---|
| `xcorr` | raw area and both FWHH values, with the resolution warning |
| `kreuz` | `kreuz_metric` = 1 − Kreuz SPIKE-distance |
| `index` | `synchrony_index` alone |
| `all` | everything — the full Axion Neural Metric Tool set |

`kreuz_metric` separates *different* cultures poorly, but it is parameter-free, which
makes it the right choice for tracking change **within one culture over time** — a
before/after stimulation comparison, for instance, where there is no window to be
accused of having chosen favourably.

```bash
# default (core): normalised correlogram area + Synchrony Index, 20 ms window
cascade --manufacturer [manufacturer] --input [recording] --output [dir]

# millisecond timescale
cascade --manufacturer [manufacturer] --input [recording] --sync_window 5 --output [dir]

# everything, including the O(N^2) spike-train measures
cascade --manufacturer [manufacturer] --input [recording] --sync_metrics all --output [dir]
```

Outputs: `synchrony_<well>.tsv`, plus `crosscorrelogram_<well>.png` and
`crosscorrelogram_norm_<well>.png` with the synchrony window marked by grey
vertical lines, as the Axion display draws them.

**The pooled correlogram** is computed over all unique electrode pairs in the
frequency domain, using the identity
`Σ_{i<j} conj(F_i)·F_j = (|Σ F_i|² − Σ|F_i|²)/2`, so the cost is O(N·T log T) rather
than O(N²·T log T). On the 341-active-channel MaxWell recording that avoids ~58,000
pairwise products.

**Why the timing parameters are absolute, not inferred.** `--sync_window`,
`--sync_bin`, `--sync_max_lag` and `--sync_sigma` are all absolute durations in
milliseconds. They are deliberately **not** derived from the sampling rate or the
mean ISI, unlike the avalanche bin width, which `infer_dt()` adapts to the data.
Sampling rates differ substantially between manufacturers (12.5 kHz on the Axion
sample, 19.75 kHz on the 3Brain sample), so a sampling-rate-derived correlogram bin
would make "area under the normalised cross-correlation at a 20 ms window" mean a
different physical quantity on each device. Keeping these absolute is what makes the
values comparable across platforms.

**Choosing `--sync_bin`: check that FWHH has converged.** Full width at half height
cannot resolve a peak narrower than the correlogram bin. Measured on the sample
recordings:

| Recording | 1/SR | FWHH @ 1.0 ms | @ 0.5 | @ 0.2 | @ 0.1 ms |
|---|---|---|---|---|---|
| Axion | 0.080 ms | 144.5 | 141.4 | 133.2 | 121.1 |
| 3Brain | 0.051 ms | 1.275 | 0.661 | 0.483 | 0.392 |
| MaxWell | 0.050 ms | 1.053 | 0.568 | 0.298 | 0.191 |

Axion's value barely moves across a 10× change in bin width — that is a converged
measurement of a genuinely broad (~130 ms) peak, and the 1 ms default is fine.
3Brain and MaxWell instead track the bin width (MaxWell at almost exactly 2× the bin
at every step), which means their true peaks are **narrower than the bin** and the
reported FWHH is a resolution limit, not a network property. For those recordings
FWHH is only a lower bound until it stops changing as `--sync_bin` shrinks, and at
0.1 ms it still has not converged — their synchrony is sub-100 µs, approaching the
0.05 ms sampling period, so spike-time quantisation may be the real floor.

Practical rule: halve `--sync_bin` until FWHH stops changing. The area metrics are
far less sensitive to bin width and are the better cross-recording comparison when
peaks are this sharp. Note that very fine bins are expensive — 0.1 ms over a 300 s
recording is 3 million bins per channel.

**Interpreting FWHH.** Full width at half height is only meaningful when the
correlogram has a genuine peak. For asynchronous activity the correlogram is flat
and the reported FWHH describes a noise fluctuation, not network structure — check
the correlogram plot before interpreting it. Verified against theory: for Gaussian
jitter of standard deviation *j* per train, the pooled peak has standard deviation
*j*·√2 and FWHH 2.355·*j*·√2, which the implementation reproduces across a 20× range
of *j*.

**Kreuz SPIKE-distance: install PySpike.** CASCADE defers to
[PySpike](https://github.com/mariomulansky/PySpike) when it is available — that is the
reference implementation, written by an author of Kreuz et al. (2013), and the values
then match it to machine precision:

```bash
pip install pyspike          # or: pip install "cascade-mea[synchrony]"
```

Without PySpike, CASCADE falls back to its own port. The fallback satisfies the
paper's guarantees (exactly 0 for identical trains, time-scale invariant to machine
precision, monotonic with mismatch) but differs from the reference by up to 2.7e-2 on
structured spike trains, because the paper specifies the dissimilarity profile
(Eqs 16-19) precisely while leaving the first and last interval underspecified beyond
"add auxiliary spikes at t=0 and t=T". The active path is recorded in the output as
`kreuz_backend`, so a fallback number never looks like a reference one. Use PySpike
for anything you intend to publish.

**Kreuz SPIKE-distance** is parameter-free and time-scale independent, so it ignores
`--sync_window` entirely. It uses the *improved* dissimilarity profile of Kreuz et
al. (2013) Eq 19 — the 2013 paper exists because the original definition produced
spuriously high instantaneous values for event-like firing patterns.

**Synchrony Index** is the Schreiber-style Gaussian-kernel correlation as presented
in Paiva et al. (2010) Eq 9, evaluated with the closed form of Eq 11, so no binning
is involved. Note that the derived dissimilarity 1 − *r* is a **semi-metric**: it is
symmetric and positive but does not satisfy the triangle inequality.

**Cost.** `kreuz` and `index` are O(N²) in channels and dominate the runtime on dense
arrays. The default `core` set includes the Synchrony Index, which is why
`--sync_max_channels` caps the O(N²) measures at 64 channels by default. Subsampling uses a fixed
seed and is recorded in the output as `kreuz_subsampled` / `index_subsampled`, so a
subsampled result never looks like a full one.

---

## Power-law normalisation and `--simple_PL`

Since v0.9.1 the avalanche size and duration exponents (τ, α) are fitted with the
discrete truncated power law normalised over **every integer** between the
cut-offs:

```
p(y) = y^-a / Σ_{k=xmin}^{xmax} k^-a
```

This is the convention of Clauset et al. (2009) and of `plmle.py` in the NCC
toolbox (Marshall et al. 2016).

CASCADE ≤ 0.9.0 summed only over the avalanche values that actually occurred in
the recording, the normalisation used by the Habibollahi et al. (2023) reference
`tplfit.py` (their paper's Eq. 1 states the full-range sum). That makes
the normalising constant depend on the sample, which biases the exponent and also
corrupts the KS-based xmin selection. The bias is largest when avalanches are few:

| avalanches | old fit | corrected fit | true exponent |
|-----------:|--------:|--------------:|--------------:|
| 421 | 1.295 | 1.502 | 1.5 |
| 898 | 1.372 | 1.498 | 1.5 |
| 9,197 | 1.484 | 1.496 | 1.5 |
| 27,288 | 1.499 | 1.500 | 1.5 |

Below roughly 1,000 avalanches the old exponent was off by 0.1–0.2; above ~10,000
the two agree.

**This changes τ, α, β_pred, DCC and SCE** relative to v0.9.0 and earlier. β_fit,
branching ratio, burst metrics, connectivity and entropy are unaffected.

To regenerate numbers published with v0.9.0 or earlier:

```bash
# Corrected fit (default) — use this for all new analyses
cascade --manufacturer [manufacturer] --input [input_recording] --output [output_dir]

# Reproduce pre-0.9.1 numbers exactly (prints a warning)
cascade --manufacturer [manufacturer] --input [input_recording] --simple_PL --output [output_dir]
```

`--simple_PL` exists only for reproducibility and provenance. It is not an
alternative model — both paths assume the same distribution, and the legacy path
simply evaluates its normalising constant incorrectly. Do not use it for new work.

---

## Reproducing the published v0.8.0 results

The published article used CASCADE v0.8.0. Two methods have been corrected since
(next section). To regenerate the v0.8.0 numbers exactly with the current version,
add both switches:

```bash
cascade --manufacturer [manufacturer] --input [input_recording] --simple_PL --simple_BR --output [output_dir]
```

`--simple_PL` restores the old power-law normalisation (τ, α, β_pred, DCC, SCE) and
`--simple_BR` the old branching ratio. Each also works on its own. Both print a
warning, because neither should be used for new analyses. The v0.8.0 release itself
remains available on GitHub and Zenodo.

---

## Differences from the published article

The article (Avila et al. 2026) describes v0.8.0. Since then:

| | v0.8.0 (article) | Current default | Old behaviour |
|---|---|---|---|
| **Power-law normalisation** (τ, α, β_pred, DCC, SCE) | Sum over the observed avalanche values only; biases the exponent below ~1,000 avalanches | Sum over every integer between the cut-offs (Clauset et al. 2009) | `--simple_PL` |
| **Branching ratio** | Single-lag regression of N(t+1) on N(t), biased toward 0 under subsampling. The article cites Wilting & Priesemann (2018), whose estimator v0.8.0 did not implement | Multi-lag MR estimator of Wilting & Priesemann (2018), validated against their `mrestimator` toolbox, so the software now matches the article's description | `--simple_BR` |
| **Upper-tail trimming** | None | None; opt-in `--pl_exclude` for comparison with the Habibollahi et al. (2023) pipeline | — |

Two further points on attribution:

- `--perc 0.25` is the default of the Habibollahi et al. (2023) reference code
  (`AV_analysis_BurstT`). Their paper reports a threshold of 40% of the median activity.
- The power-law fit follows the KS-minimising xmin search of that reference code, with
  the normalisation corrected as above; `--pl_exclude` reproduces its Python version.

Details and the measured effect of each change are in `CHANGELOG.md`.

---

## Project structure

```
environment.yml              conda environment
pyproject.toml               package metadata and CLI entry points
conda-recipe/meta.yaml       conda build recipe
CITATION.cff                 citation metadata (GitHub "Cite this repository")
CHANGELOG.md                 what changed in each version, and why
cascade/
├── run_analysis.py          cascade command
├── scanner.py               cascade-info command
├── ltpltd.py                cascade-ltpltd command
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
│   ├── synchrony.py         cross-correlogram, Kreuz, Synchrony Index
│   ├── entropy.py           spatial entropy
│   └── plots.py             figure generation
│
├── utils/
│   ├── common.py            MEARecording dataclass, helpers
│   ├── io.py                CSV/TSV writers, criticality report, run log
│   ├── citation.py          --citation text
│   └── inspect_file.py      HDF5 structure inspection
│
└── tests/
    └── self_test.py         the 14 checks run by cascade --test
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

Run `python -m cascade.utils.inspect_file dump [yourfile]` to explore file structure.

---

## Citing CASCADE

If you use CASCADE, please cite the article and state the CASCADE version you used
(`cascade --version`):

> Avila F, Kim J, Park J-C, Kang R, Lee H, Park S-H. CASCADE: criticality avalanche spike cross-platform analysis detection engine, a multi-manufacturer MEA bash analysis pipeline. *Bioinformatics*. 2026;42(8):btag570. https://doi.org/10.1093/bioinformatics/btag570

`cascade --citation` prints this with a BibTeX entry, and GitHub's **Cite this repository** button reads `CITATION.cff`.

---

## References

- Habibollahi, F., Kagan, B. J., Burkitt, A. N., & French, C. (2023). Critical dynamics arise during structured information presentation within embodied in vitro neuronal networks. *Nature Communications*, **14**(1), 5287. https://doi.org/10.1038/s41467-023-41020-3

- Beggs, J. M., & Plenz, D. (2003). Neuronal Avalanches in Neocortical Circuits. *The Journal of Neuroscience*, **23**(35), 11167–11177. https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003

- Wilting, J., & Priesemann, V. (2018). Inferring collective dynamical states from widely unobserved systems. *Nature Communications*, **9**(1), 2325. https://doi.org/10.1038/s41467-018-04725-4

- Friedman, N., Ito, S., Brinkman, B. A. W., Shimono, M., DeVille, R. E. L., Dahmen, K. A., Beggs, J. M., & Butler, T. C. (2012). Universal Critical Dynamics in High Resolution Neuronal Avalanche Data. *Physical Review Letters*, **108**(20), 208102. https://doi.org/10.1103/PhysRevLett.108.208102

- Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-Law Distributions in Empirical Data. *SIAM Review*, **51**(4), 661–703. https://doi.org/10.1137/070710111

- Marshall, N., Timme, N. M., Bennett, N., Ripp, M., Lautzenhiser, E., & Beggs, J. M. (2016). Analysis of Power Laws, Shape Collapses, and Neural Complexity: New Techniques and MATLAB Support via the NCC Toolbox. *Frontiers in Physiology*, **7**, 250. https://doi.org/10.3389/fphys.2016.00250

- Kreuz, T., Chicharro, D., Houghton, C., Andrzejak, R. G., & Mormann, F. (2013). Monitoring spike train synchrony. *Journal of Neurophysiology*, **109**(5), 1457–1472. https://doi.org/10.1152/jn.00873.2012

- Paiva, A. R. C., Park, I., & Príncipe, J. C. (2010). A comparison of binless spike train measures. *Neural Computing and Applications*, **19**(3), 405–419. https://doi.org/10.1007/s00521-009-0307-6

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
