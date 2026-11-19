# CASCADE v0.9.6

**Criticality Avalanche Spike Cross-platform Analysis Detection Engine**

<div align="center">
<img width="250" height="250" alt="CASCADE_Logo" src="https://github.com/user-attachments/assets/7a029bf6-93b7-4185-8859-888672b4bb65" />
</div>

CASCADE analyses multi-electrode array (MEA) spike recordings from six manufacturers with one command. Each file is converted into a common format, so firing rate, network bursts, neuronal avalanches, criticality, functional connectivity, synchrony and spatial entropy are computed the same way for every device. Figures are saved as PNG (300 dpi) and numbers as TSV.

Published in *Bioinformatics* (2026), [doi:10.1093/bioinformatics/btag570](https://doi.org/10.1093/bioinformatics/btag570).

See [Citing CASCADE](#citing-cascade).

> **This is v0.9.6.** The article describes v0.8.0. The pipeline is the same, with new features and two updated calculations; the article's methods remain available as options. See [What's new since the article](#whats-new-since-the-article).

## Contents

- [Tools](#tools)
- [Installation](#installation)
- [Quick start](#quick-start)
  - [Supported manufacturers](#supported-manufacturers)
  - [Recordings with several wells](#recordings-with-several-wells)
- [Customising the analysis](#customising-the-analysis)
  - [Analysing part of a recording](#analysing-part-of-a-recording)
  - [Removing electrodes with few spikes](#removing-electrodes-with-few-spikes)
- [Output files](#output-files)
- [All flags](#all-flags)
- [Key metrics](#key-metrics)
- [Noise filtering](#noise-filtering)
- [Synchrony metrics](#synchrony-metrics)
- [Differences from the published article](#differences-from-the-published-article)
- [Project structure](#project-structure)
- [Adding a new manufacturer](#adding-a-new-manufacturer)
- [What's new since the article](#whats-new-since-the-article)
- [Citing CASCADE](#citing-cascade)
- [References](#references)
- [Funding](#funding)
- [Acknowledgement](#acknowledgement)

---

## Tools

CASCADE has three tools.

| Command | What it does |
|---|---|
| `cascade-info` | Inspect a recording or a folder: duration, channels, active electrodes |
| `cascade` | Full pipeline |
| `cascade-ltpltd` | Classify electrodes as LTP/LTD from firing-rate change between two `cascade` runs |

`cascade-info` is used for scanning the overall information of an MEA recording file before the analysis. It prints the manufacturer, the wells in the file, the sampling rate, the recording duration, the analysed time window, the total and active electrodes (with `--ae`), and the number of stimulations (FinalSpark). Given a folder, it scans every file of that manufacturer and can save one row per file to `info_summary.tsv`.

`cascade` is the main tool. It preprocesses the spikes and calculates spike-based metrics such as firing rate, inter-spike interval (ISI) and refractory violations, network bursts, active area, functional connectivity, synchrony, neuronal avalanches, criticality (τ, α, β, DCC, branching ratio and shape collapse) and spatial entropy, together with rasters and heatmaps.

`cascade-ltpltd` is used to compare a control (`--baseline`) and a treated (`--measure`) MEA recording. It is the only tool that needs two inputs, which are the `rates_[well].tsv` files from two `cascade` runs, not the recordings themselves.

---

## Installation

### GitHub
```bash
git clone https://github.com/FA387/cascade-mea.git
cd cascade-mea
conda env create -f environment.yml
conda activate cascade-mea
pip install .
```

### conda
```bash
conda install fa387::cascade-mea
```

### pip
```bash
pip install cascade-mea
```

conda and pip install the latest release uploaded there, which can be older than this README. If `cascade --version` prints an older version, install from GitHub as above.

Then run the self-test:

```bash
cascade --test       # 16 self-tests on synthetic data
cascade --version    # prints the installed version
cascade --help       # lists every flag
```

Expected output of `cascade --test`:

```
CASCADE v0.9.6 — Built-in Self-Test
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
  Stimulation blanking (--stim_blank) .................... PASS
  FinalSpark loader (synthetic export) ................... PASS
===================================================================
  All 16 tests passed. CASCADE is installed correctly.
```

If you do not have your own MEA recordings, you can download the sample recordings from the Zenodo repository: https://doi.org/10.5281/zenodo.19161615. These are the same six recordings from five manufacturers used in the [Bioinformatics paper](https://doi.org/10.1093/bioinformatics/btag570), and the repository also contains the expected `cascade` output for each one. In the paper (Table 1, v0.8.0), each recording took 4.8–20.3 s to analyse with a peak memory of 249–1570 MB, including an 8.6 GB Cortical Labs file (12.4 s, 353 MB). Later versions run longer because of the new modules.

---

## Quick start

The minimum requirements to use `cascade-info` and `cascade` are the following:

1. The tool (`cascade-info` or `cascade`).
2. The MEA manufacturer, `--manufacturer [manufacturer]` (see [Supported manufacturers](#supported-manufacturers)).
3. The input file, the MEA recording, `--input [MEA_recording]`.
   - Several manufacturers use MEA plates with multiple wells. If the recording measured the whole plate, specify the well with `--well [well_id]` (see [Recordings with several wells](#recordings-with-several-wells)).

Although the output folder does not need to be specified, since the results are then written to `mea_out` in the current directory, we highly recommend using `--output [output_dir]` for better file management.

For `cascade-ltpltd`, the inputs are the `rates_[well].tsv` files generated by `cascade`: the baseline (control) `--baseline [baseline_dir]/rates_[well].tsv` and the measured (treated) `--measure [measure_dir]/rates_[well].tsv`.

```bash
# 1. cascade-info: inspect a recording before the analysis
cascade-info --manufacturer [manufacturer] --input [MEA_recording]
cascade-info --manufacturer [manufacturer] --input [folder] --ae 5 --output [output_dir]   # every file in a folder, saved to info_summary.tsv

# 2. cascade: full analysis
cascade --manufacturer [manufacturer] --input [MEA_recording] --output [output_dir]                     # minimum
cascade --manufacturer axion --input [MEA_recording] --well [well_id] --output [output_dir]             # select a well of a multi-well plate
cascade --manufacturer [manufacturer] --input [MEA_recording] --ae 5 --perc 0.25 --output [output_dir]  # settings used in the article

# 3. cascade-ltpltd: compare two cascade runs
cascade-ltpltd --baseline [baseline_dir]/rates_[well].tsv --measure [measure_dir]/rates_[well].tsv --output [output_dir]
```

Without `--well`, Axion merges all wells into one network; analyse wells separately for criticality.

### Supported manufacturers

Each MEA manufacturer saves its recordings in its own file format. Most of them are based on the Hierarchical Data Format (HDF5), with extensions such as `.h5`, `.hdf5` or `.bxr`. However, the data structure inside the file differs between manufacturers, and each manufacturer documents it in its own resources. Axion BioSystems is the exception, since its `.raw` and `.spk` files are a proprietary binary format that is read with Axion's MATLAB tools. Thus, for Axion, CASCADE reads the spike-list `.csv` file exported from the Axion software.

CASCADE v0.9.6 processes MEA recordings from **six** manufacturers. Each manufacturer has its own identifier for `--manufacturer`, which selects the correct spike loader. Since `.h5` is used by three devices, check which one recorded the file. A wrong choice stops with an error that names the right device when the file identifies it.

| Manufacturer | `--manufacturer` | File extension | What is read |
|---|---|---|---|
| 3Brain BrainWave | `3brain` | `.bxr` | SpikeTimes, SpikeChIdxs, SpikeUnits, SpikeForms |
| Axion BioSystems | `axion` | `.csv` | AxIS spike-list export |
| Cortical Labs | `cl` | `.h5` | `/spikes` table (timestamp, channel, waveform) |
| FinalSpark Neuroplatform | `finalspark` | `_package.hdf5` | `fsNNN_package.hdf5`: events table (time, electrode, peak voltage) and stimulation times |
| Maxwell Biosystems | `maxwell` | `.h5` | `data_store/<stream>/spikes` (frameno, channel, amplitude) |
| Multi Channel Systems | `mcs` | `.h5` | MCS HDF5 v3 SegmentStream spike cutouts |

#### Recordings with several wells

Several MEA devices use a chip or plate with multiple wells (separate chambers), so the samples and their media are separated from one another. By default, the MEA generates only a single recording file containing the data from all the wells, and each spike carries a label that indicates its well.

Thus, `cascade` has the flag `--well [well_id]`, which analyses only the specified well. Without `--well`, Axion and FinalSpark recordings are analysed with all wells together as one network, while Maxwell uses the first well. Thus, users must be careful when working with an MEA recording that has multiple wells and samples.

##### Axion BioSystems

Axion BioSystems offers many types of MEA plates, which vary in the number of wells.

For example, the data used in the development and publication of CASCADE was recorded on a 6-well plate (CytoView MEA 6), with 64 electrodes per well in an 8 × 8 grid.

Here are the names of the wells and electrodes of the 6-well plate:

| Row | Column | `--well` | Electrode IDs |
|---|---|---|---|
| A | 1 | `A1` | `A1_11` – `A1_88` |
| A | 2 | `A2` | `A2_11` – `A2_88` |
| A | 3 | `A3` | `A3_11` – `A3_88` |
| B | 1 | `B1` | `B1_11` – `B1_88` |
| B | 2 | `B2` | `B2_11` – `B2_88` |
| B | 3 | `B3` | `B3_11` – `B3_88` |

The two digits after the underscore give the electrode's position in the 8 × 8 grid (1–8 each).

| | 1 | 2 | 3 |
|---|---|---|---|
| **A** | A1 | A2 | A3 |
| **B** | B1 | B2 | B3 |

For example, a toxicology study is conducted by MEA measurement on a plate from Axion BioSystems (`axion`). In the beginning, a spontaneous (control) recording was made, which generated the control MEA recording `control.csv`. Then six drugs were added, each in a different well. After incubation, the plate was measured again, which generated `treatment.csv`.

1. Drug 1 was added in well A1
2. Drug 2 was added in well A2
3. Drug 3 was added in well A3
4. Drug 4 was added in well B1
5. Drug 5 was added in well B2
6. Drug 6 was added in well B3

The following shows how to analyse the effect of drug 1 using CASCADE:

```bash
# 1. Analyse well A1 before the addition of drug 1 (baseline/control)
cascade --manufacturer axion --input control.csv --well A1 --output control_drug1

# 2. Analyse well A1 after the addition of drug 1 (treatment/measure)
cascade --manufacturer axion --input treatment.csv --well A1 --output measure_drug1

# 3. Compare the control and treatment of drug 1
cascade-ltpltd --baseline control_drug1/rates_A1.tsv --measure measure_drug1/rates_A1.tsv --output ltpltd_drug1
```

Step 1 generates the results of the control, and step 2 generates the results of the treatment. `cascade-ltpltd` then compares the `rates_A1.tsv` files from both steps and classifies each electrode as LTP, LTD or no change. For drugs 2–6, change `--well` and the file names (for example `--well A2` and `rates_A2.tsv`).

##### FinalSpark

The FinalSpark Neuroplatform uses a custom MEA, in which each MEA has 4 sites, each holding one organoid recorded by 8 electrodes (32 electrodes in total).

Here are the names of the sites and electrodes in the FinalSpark recording:

| Organoid | `--well` | Electrode IDs |
|---|---|---|
| 1 | `0` | 0–7 |
| 2 | `1` | 8–15 |
| 3 | `2` | 16–23 |
| 4 | `3` | 24–31 |

These IDs are for the first MEA of the Neuroplatform; each following MEA adds 32 (for example, 32–39 for organoid 1 of the second MEA).

For example, a FinalSpark recording named `example_package.hdf5`:

```bash
# 1. Analyse site 0 (organoid 1) only
cascade --manufacturer finalspark --input example_package.hdf5 --well 0 --output site0_organoid1

# 2. Analyse all four sites (all 32 electrodes) together
cascade --manufacturer finalspark --input example_package.hdf5 --output all_sites
```

FinalSpark support is new in v0.9.6 and has so far been tested on one MEA recording sample (`fs437_package.hdf5`).

The sample recording was kindly provided by FinalSpark (see [Acknowledgement](#acknowledgement)).

---

## Customising the analysis

CASCADE has several flags that customise the spike preprocessing and analysis. This section shows the two most common ones with examples.

### Analysing part of a recording

By default, `cascade` analyses the MEA recording from the start (0 s) to the end of the recording.

However, `cascade-info` and `cascade` can analyse a particular part of the recording with `--start [time in s]` and `--dur [time in s]`. Note that `--dur` is the length of the window, not its end time.

For example, a 5-minute (300 s) recording from a Cortical Labs device (`cl`) named `time_example.h5`:

```bash
# 1. Analyse the whole recording
cascade --manufacturer cl --input time_example.h5 --output output_1

# 2. Analyse from 0 to 120 s
cascade --manufacturer cl --input time_example.h5 --start 0 --dur 120 --output output_2

# 3. Analyse from 120 to 180 s (start at 120 s, for 60 s)
cascade --manufacturer cl --input time_example.h5 --start 120 --dur 60 --output output_3
```

The time unit is seconds. A smaller time scale in milliseconds can be achieved by converting the milliseconds to seconds:

```bash
# 4. Analyse from 0 to 500 ms (0.5 s)
cascade --manufacturer cl --input time_example.h5 --start 0 --dur 0.5 --output output_4
```

A window this short usually has too few spikes for bursts, avalanches and criticality, so use it for rasters and firing rates only.

This feature is useful when an MEA recording has an electrical stimulation in the middle of the recording. For example, a 10-minute (600 s) recording from a Cortical Labs device (`cl`) named `stim_example.h5`, in which stimulation starts at 300 s:

```bash
# 1. Before stimulation (0 to 300 s)
cascade --manufacturer cl --input stim_example.h5 --start 0 --dur 300 --output before_stim

# 2. During stimulation (300 to 600 s)
cascade --manufacturer cl --input stim_example.h5 --start 300 --dur 300 --output during_stim

# 3. Or both in one run, with an extra figure comparing the two epochs
cascade --manufacturer cl --input stim_example.h5 --epochA 0,300 --epochB 300,300 --output stim_compare
```

The `rates_[well].tsv` files from runs 1 and 2 can then be compared with `cascade-ltpltd`.

### Removing electrodes with few spikes

There is an option to exclude weak or inactive electrodes, with the approach that an electrode which detected at least X spikes per minute is trustworthy, while an electrode below that value is not. Thus, the `--ae [spikes per minute]` flag removes every electrode that detected fewer spikes per minute (<) than the set value. The rate is calculated over the analysed window, so it changes with `--start` and `--dur`.

For example, a 2-minute `axion` recording named `ae_filter_example.csv`, in which electrode A1_11 has the following spikes:

| Minute | Spikes | Average spikes per minute |
|---|---|---|
| 1 | 4 | 4 / 1 = 4 |
| 2 | 6 | (4 + 6) / 2 = 5 |

```bash
# 1. A1_11 is included: whole recording, 5 spikes/min ≥ --ae 5
cascade --manufacturer axion --input ae_filter_example.csv --ae 5 --output output_1

# 2. A1_11 is excluded: first minute only, 4 spikes/min < --ae 5
cascade --manufacturer axion --input ae_filter_example.csv --ae 5 --dur 60 --output output_2

# 3. A1_11 is excluded: whole recording, 5 spikes/min < --ae 6
cascade --manufacturer axion --input ae_filter_example.csv --ae 6 --output output_3

# 4. A1_11 is included: the default keeps every electrode, the same as --ae 0
cascade --manufacturer axion --input ae_filter_example.csv --output output_4
cascade --manufacturer axion --input ae_filter_example.csv --ae 0 --output output_4
```

A higher `--ae` keeps only electrodes with many spikes, which removes silent and weak electrodes from the burst, avalanche and connectivity analyses. However, a low firing rate does not mean an electrode is untrustworthy, since some healthy neurons fire rarely. Removing them also changes the number of electrodes that every network metric is calculated on, so τ, α, the branching ratio and connectivity can change with `--ae`. Thus, there is no single correct value. The article used `--ae 5`. Check how many electrodes stay active at a few values with `cascade-info --ae`, choose one, use the same value for every recording you compare, and report it in your methods.

---

## Output files

### cascade-info output

`cascade-info` prints its results to the terminal. With `--output`, it also writes `[output_dir]/info_summary.tsv` with one row per file.

### cascade output

```
[output_dir]/
├── criticality_report_<well>.txt     τ, α, β, DCC, BR (with checks), SCE: the main result
├── cascade_run.log                   command, version, terminal output, SUCCESS/FAILED
│
├── rates_<well>.tsv                  per-channel firing rate per --bin (input for cascade-ltpltd)
├── isi_metrics_<well>.tsv            per-channel ISI stats and refractory violations
├── bursts_<well>.tsv                 network bursts: start, end, duration, spikes
├── avalanches_<well>.tsv             avalanches used for all criticality metrics
├── branching_mr_<well>.tsv           MR estimator: r_k per lag and both fits
├── branching_points_<well>.tsv       N(t), N(t+1) pairs for --simple_BR
├── corr_<well>.tsv                   Pearson correlation matrix (active channels)
├── degrees_<well>.tsv                connections per channel above --conn_thr
├── synchrony_<well>.tsv              synchrony metrics and the parameters used
├── active_area_<well>.tsv            fraction of active electrodes per window
├── entropy_<well>.tsv                spatial entropy per window
│
├── raster_max_<well>.png             raster, all electrode slots
├── raster_active_<well>.png          raster, active electrodes only
├── rate_heatmap_max_<well>.png       firing-rate heatmap, all slots
├── rate_heatmap_active_<well>.png    firing-rate heatmap, active only
├── population_rate_<well>.png        mean firing rate over time
├── bursts_timeline_<well>.png        network bursts over time
├── avalanches_<well>_size.png        avalanche size histogram
├── avalanches_<well>_duration.png    avalanche duration histogram
├── loglog_size_<well>.png            size distribution with power-law fit
├── loglog_duration_<well>.png        duration distribution with power-law fit
├── shape_collapse_<well>.png         avalanche shape collapse
├── branching_mr_<well>.png           r_k vs lag with fits (should decay to 0)
├── corr_heatmap_max_<well>.png       correlation matrix, all slots
├── corr_heatmap_active_<well>.png    correlation matrix, active only
├── degree_hist_<well>.png            node degree histogram
├── crosscorrelogram_<well>.png       pooled cross-correlogram
├── crosscorrelogram_norm_<well>.png  normalised cross-correlogram
├── entropy_curve_<well>.png          spatial entropy over time
├── active_area_curve_<well>.png      active electrode fraction over time
├── wave_templates_<well>.png         mean waveform per channel (3Brain with SpikeForms)
└── epoch_compare_<well>.png          epoch A vs B (with --epochA and --epochB)
```

Report the CASCADE version from `cascade_run.log` with your results, since criticality values differ between versions.

### cascade-ltpltd output

```
[output_dir]/
├── ltpltd_[label].tsv            per-electrode call: LTP / LTD / NoChange
├── ltpltd_[label]_bar.png        % change per electrode
├── ltpltd_[label]_pie.png        LTP / LTD / NoChange proportions
└── ltpltd_all_conditions.tsv     all conditions (several --measure files)
```

---

## All flags

Flags are grouped by workflow step, in the same order as `--help`.

| Step | Decides |
|---|---|
| **1. Setup** | Which file, which device, where results go |
| **2. Preprocessing** | Which part of the recording and which channels to trust |
| **3. Analysis** | Parameters of each calculation |
| **4. Output** | Extra figure formats and terminal output |
| **5. Misc** | Information only, no analysis |

Most runs need only step 1. Add step 2 filters for noisy recordings, change step 3 only with a reason, and report any non-default value in your methods.

### cascade-info

#### 1. Setup

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `axion`, `cl`, `finalspark`, `maxwell`, `mcs` |
| `--input` | required | Recording file, or a folder (scans that manufacturer's file type) |
| `--well` | first | Well or stream. Axion: all wells merged unless given. FinalSpark: site `0`–`3`, default all 32 electrodes |
| `--sr` | none | Sampling rate (Hz), only for a plain CSV with frame numbers |
| `--output` | off | Folder for `info_summary.tsv` |

#### 2. Preprocessing

| Flag | Default | Description |
|---|---|---|
| `--start` | `0` | Window start (s) |
| `--dur` | full | Window duration (s) |
| `--ae` | `0` | Count channels with ≥ N spikes/min as active (nothing removed) |

#### 5. Misc

| Flag | Description |
|---|---|
| `--help` / `--version` / `--citation` | Show flags / version / citation |

### cascade

#### 1. Setup

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `axion`, `cl`, `finalspark`, `maxwell`, `mcs` |
| `--input` | auto | Recording file. If omitted, the first `.bxr`/`.h5` in the current folder |
| `--well` | first | Well or stream. Axion: all wells merged unless given. FinalSpark: site `0`–`3`, default all 32 electrodes |
| `--output` | `mea_out` | Output folder |

#### 2. Preprocessing

| Flag | Default | Description |
|---|---|---|
| `--start` | `0` | Window start (s) |
| `--dur` | full | Window duration (s) |
| `--stim_blank` | off | Remove spikes within N ms after each stimulation, on all channels, before the other filters (try `20`). Needs stimulation times: FinalSpark |
| `--ae` | off | Keep channels with ≥ N spikes/min |
| `--max_rate` | off | Remove channels above N spikes/min (try `1000`) |
| `--max_rpv` | off | Remove channels with > N % of ISIs < 2 ms (try `1.0`) |
| `--max_amplitude` | off | Remove channels with any spike > N µV peak-to-peak (try `200`). Maxwell stores a single-peak value, if any |

#### 3. Analysis

**Firing rate and connectivity**

| Flag | Default | Description |
|---|---|---|
| `--bin` | `1.0` | Bin for rates, heatmaps and connectivity (s) |
| `--conn_thr` | `0.3` | Pearson \|r\| threshold for a connection |

**Network bursts**

| Flag | Default | Description |
|---|---|---|
| `--burst_min_spikes` | `5` | Minimum spikes per burst |
| `--burst_max_gap` | `0.1` | Maximum gap within a burst (s) |

**Synchrony** (all times in ms)

| Flag | Default | Description |
|---|---|---|
| `--sync_metrics` | `core` | `none`, `core`, `xcorr`, `kreuz`, `index`, `all`. See [Synchrony metrics](#synchrony-metrics) |
| `--sync_window` | `20.0` | Window around zero lag for the area metrics |
| `--sync_bin` | `1.0` | Correlogram bin |
| `--sync_max_lag` | `500.0` | Maximum correlogram lag; must contain the peak |
| `--sync_sigma` | = window | Gaussian σ for the Synchrony Index |
| `--sync_max_channels` | `64` | Channel cap for O(N²) metrics; `0` = no cap |

**Avalanche detection**

| Flag | Default | Description |
|---|---|---|
| `--ava_size` | `2` | Minimum avalanche size (spikes) |
| `--ava_dur` | `1` | Minimum avalanche duration (bins) |
| `--aval_dt` | auto | Bin width (s); auto = mean ISI, clipped to 1–50 ms |
| `--perc` | off | Population activity percentile threshold. `0.25` is the Habibollahi et al. (2023) reference-code default |

**Power-law fit (τ, α, DCC)**

| Flag | Default | Description |
|---|---|---|
| `--burst_min` | `12` | Upper limit of the xmin search, sizes |
| `--t_min` | `4` | Upper limit of the xmin search, durations |
| `--pl_xmin_S` | auto | Fix the size xmin; the automatic fit is still reported |
| `--pl_xmin_D` | auto | Fix the duration xmin |
| `--pl_exclude` | off | Trim the upper tail first (Habibollahi et al. 2023 Python code). For comparison with that pipeline only |
| `--skip_rare_sizes` | off | Drop sizes seen once from the plot line (fit unchanged) |
| `--skip_rare_durations` | off | Drop durations seen once from the plot line (fit unchanged) |

**Shape collapse, branching ratio, entropy**

| Flag | Default | Description |
|---|---|---|
| `--sc_Dmin` | `4` | Shortest duration in shape collapse (bins) |
| `--sc_Dmax` | `20` | Longest duration in shape collapse (bins) |
| `--br_tmax` | `2000` | Longest lag for the MR estimator (ms); the report warns if too short |
| `--entropy_win` | `10` | Window for entropy and active-area curves (s) |

**Epoch comparison** (needs both)

| Flag | Default | Description |
|---|---|---|
| `--epochA` | off | `start,duration` in s, e.g. `0,60` |
| `--epochB` | off | `start,duration` in s, e.g. `120,60` |

**Reproducing v0.8.0** (not for new analyses)

| Flag | Default | Description |
|---|---|---|
| `--simple_PL` | off | Old power-law normalisation (biased below ~1,000 avalanches) |
| `--simple_BR` | off | Old single-lag branching ratio (biased toward 0) |

#### 4. Output

| Flag | Default | Description |
|---|---|---|
| `--svg` / `--jpg` | off | Also save figures as SVG/JPG |
| `--verbose` | off | Print every step, including channels removed by each filter |

#### 5. Misc

| Flag | Description |
|---|---|
| `--help` | Show all flags |
| `--version` | Print the version |
| `--citation` | Print the citation (text and BibTeX) |
| `--test` | Run the 16 self-tests |

### cascade-ltpltd

#### 1. Setup

| Flag | Default | Description |
|---|---|---|
| `--baseline` | required | Baseline `rates_*.tsv` |
| `--measure` | required | One or more `rates_*.tsv` to compare |
| `--output` | `ltpltd_out` | Output folder |

#### 3. Analysis

| Flag | Default | Description |
|---|---|---|
| `--threshold` | `20` | % change from baseline for LTP/LTD |

#### 4. Output

| Flag | Default | Description |
|---|---|---|
| `--svg` / `--jpg` | off | Also save figures as SVG/JPG |
| `--verbose` | off | Print the per-electrode table |

#### 5. Misc

| Flag | Description |
|---|---|
| `--help` / `--version` / `--citation` | Show flags / version / citation |

---

## Key metrics

| Metric | Definition | Critical value |
|---|---|---|
| **τ** | Avalanche size exponent (truncated power-law MLE, xmin by KS) | ~1.5 |
| **α** | Avalanche duration exponent (same fit) | ~2.0 |
| **β_pred** | (α − 1) / (τ − 1) | |
| **β_fit** | Slope of log mean size vs log duration, xmin_D < D < xmin_D + 60 | |
| **DCC** | \|β_pred − β_fit\| | → 0 |
| **BR** | Branching ratio m, multi-lag MR estimator, r_k = b·m^k (Wilting & Priesemann 2018) | ~1 |
| **SCE** | Shape collapse error, \|β_opt − β_pred\| | → 0 |

Below ~50 avalanches the power-law fit is unreliable. Trust BR only when its three checks in the report pass.

---

## Noise filtering

Each filter removes whole channels. They run in this order: `--ae`, then `--max_rate` and `--max_rpv`, then `--max_amplitude`.

```bash
cascade --manufacturer [manufacturer] --input [MEA_recording] --output [output_dir] \
        --ae 5 --max_rate 1000 --max_rpv 1.0 --max_amplitude 200 --verbose
```

Add only the filters a first run shows you need, check how many channels each removed (`--verbose`), and use the same settings for every recording you compare. Amplitudes are in µV for all devices, but Maxwell stores a single-peak value, so calibrate the threshold per device.

**Stimulation.** If the file records stimulation times (FinalSpark), `cascade` reports how many spikes fall within 20 ms after a stimulation. Those spikes are stimulus artefacts or evoked responses, not spontaneous activity, and they can inflate avalanches. `--stim_blank N` removes every spike within N ms after each stimulation, on all channels, before the channel filters: `5` removes mainly artefacts, `20` also removes most evoked responses. In the fs437 export (802,657 stimulations), 20 % of all spikes fell within 20 ms after one, compared with 0.6 % expected by chance. `cascade-info` prints the number of stimulations in the window.

---

## Synchrony metrics

CASCADE computes the synchrony metrics of the Axion Neural Metric Tool for every manufacturer. Pearson connectivity uses 1 s bins; these metrics resolve millisecond timing.

| `--sync_metrics` | Reports |
|---|---|
| `core` (default) | `area_xcorr_norm` (normalised correlogram area) and `synchrony_index` (Paiva et al. 2010) |
| `xcorr` | Raw area and both full widths at half height (FWHH) |
| `kreuz` | `kreuz_metric` = 1 − SPIKE-distance (Kreuz et al. 2013) |
| `index` | `synchrony_index` only |
| `all` | Everything |
| `none` | Skip |

`core` keeps the two metrics that resolved on all five sample recordings and separated them best. `kreuz_metric` separates cultures poorly but is parameter-free, so it suits before/after comparisons within one culture.

- Timing flags are absolute ms, never inferred from the sampling rate, so values are comparable across devices.
- FWHH cannot resolve a peak narrower than `--sync_bin`. Halve the bin until FWHH stops changing, or compare the area metrics instead.
- For the reference SPIKE-distance, install PySpike (`pip install ".[synchrony]"`). Without it CASCADE uses its own port (up to 2.7e-2 off); `kreuz_backend` in the TSV says which was used.
- `kreuz` and `index` are O(N²) in channels and capped at `--sync_max_channels`; subsampling is recorded in the TSV.

---

## Differences from the published article

The pipeline is the same as in the article (v0.8.0). Spike loading, filtering, firing rates, bursts, connectivity and avalanche detection give the same results. Two criticality calculations were updated to follow their reference methods more closely, so τ, α, DCC, SCE and the branching ratio can differ slightly from v0.8.0.

| Calculation | v0.8.0 (article) | v0.9.6 (default) | v0.8.0 method |
|---|---|---|---|
| Power-law fit (τ, α, β_pred, DCC, SCE) | Normalised over the observed avalanche sizes | Normalised over every size in the fitted range (Clauset et al. 2009) | `--simple_PL` |
| Branching ratio (BR) | One time lag | Many time lags (Wilting & Priesemann 2018) | `--simple_BR` |

The difference is small for recordings with many avalanches (about 1% for the 27,288-avalanche Cortical Labs recording in the article) and larger for recordings with fewer than ~1,000 avalanches.

### Reproducing the published v0.8.0 results

```bash
cascade --manufacturer [manufacturer] --input [MEA_recording] --ae 5 --perc 0.25 --simple_PL --simple_BR --output [output_dir]
```

The two flags restore the v0.8.0 calculations, so the article's numbers can be regenerated exactly. For new analyses, use the defaults.

---

## Project structure

```
environment.yml             conda environment
pyproject.toml              package metadata, CLI entry points
conda-recipe/meta.yaml      conda build recipe
CITATION.cff                citation metadata
CHANGELOG.md                version history
cascade/
├── run_analysis.py         cascade
├── scanner.py              cascade-info
├── ltpltd.py               cascade-ltpltd
├── loaders/                one file per manufacturer → MEARecording
├── core/                   spikes, bursts, avalanche, criticality, connectivity, synchrony, entropy, plots
├── utils/                  MEARecording, writers, run log, citation, inspect_file
└── tests/self_test.py      the 16 checks run by cascade --test
```

## Adding a new manufacturer

1. Add `cascade/loaders/<name>_loader.py` subclassing `BaseLoader`.
2. Implement `load(path, well=None, start_s=0.0, duration_s=None) -> MEARecording` with spike times in **seconds**.
3. Register it in `LOADERS` and `EXTENSIONS` in `cascade/loaders/__init__.py`. `cascade` and `cascade-info` take their `--manufacturer` choices from `LOADERS`.

```python
from .base_loader import BaseLoader
from ..utils.common import MEARecording

class YourLoader(BaseLoader):
    @staticmethod
    def manufacturer_name() -> str:
        return "yourmanufacturer"

    def load(self, path, well=None, start_s=0.0, duration_s=None) -> MEARecording:
        # return MEARecording(sampling_rate, n_channels, channel_ids,
        #                     spike_times={ch: np.array(seconds)}, duration, metadata={...})
        ...
```

Explore an unknown file with `python -m cascade.utils.inspect_file dump [file]`.

---

## What's new since the article

### Added

- FinalSpark Neuroplatform support (`--manufacturer finalspark`). See [FinalSpark](#finalspark).
- `--stim_blank`: removes spikes right after each stimulation (artefacts and evoked responses). See [Noise filtering](#noise-filtering).
- [Synchrony metrics](#synchrony-metrics) for all manufacturers (`--sync_metrics`).
- `--max_amplitude`: removes artefact channels by spike amplitude. See [Noise filtering](#noise-filtering).
- `--pl_exclude`: opt-in upper-tail trimming from the Habibollahi et al. (2023) reference code.
- `--citation` in all three tools; the run log records the CASCADE version.
- `--input` / `--output` in all three tools; flags grouped by workflow step in `--help`.

### Updated

A default run can give slightly different criticality values than v0.8.0. See [Differences from the published article](#differences-from-the-published-article).

- Power-law fit normalisation: τ, α, β_pred, DCC, SCE.
- Branching ratio: multi-lag estimator (Wilting & Priesemann 2018).

### Fixed

- `--pl_xmin_S/D` now set the fit, not only the plot.
- MCS amplitudes are in µV.
- Single-channel recordings skip connectivity instead of failing.

The article's conclusions are not affected, because it reports no τ, α, DCC, BR or SCE values in its text or tables. Full history: [CHANGELOG.md](CHANGELOG.md).

---

## Citing CASCADE

Cite the article and the CASCADE version you used (`cascade --version`):

> Avila F, Kim J, Park J-C, Kang R, Lee H, Park S-H. CASCADE: criticality avalanche spike cross-platform analysis detection engine, a multi-manufacturer MEA bash analysis pipeline. *Bioinformatics*. 2026;42(8):btag570. https://doi.org/10.1093/bioinformatics/btag570

`cascade --citation` prints it with BibTeX. GitHub's **Cite this repository** reads `CITATION.cff`.

- Source code archive: https://doi.org/10.5281/zenodo.20840795
- Sample recordings with v0.8.0 output: https://doi.org/10.5281/zenodo.19161615

## References

- Beggs JM, Plenz D (2003). Neuronal avalanches in neocortical circuits. *J Neurosci* 23:11167–11177. https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003
- Clauset A, Shalizi CR, Newman MEJ (2009). Power-law distributions in empirical data. *SIAM Rev* 51:661–703. https://doi.org/10.1137/070710111
- Friedman N et al. (2012). Universal critical dynamics in high resolution neuronal avalanche data. *Phys Rev Lett* 108:208102. https://doi.org/10.1103/PhysRevLett.108.208102
- Habibollahi F, Kagan BJ, Burkitt AN, French C (2023). Critical dynamics arise during structured information presentation within embodied in vitro neuronal networks. *Nat Commun* 14:5287. https://doi.org/10.1038/s41467-023-41020-3
- Kreuz T et al. (2013). Monitoring spike train synchrony. *J Neurophysiol* 109:1457–1472. https://doi.org/10.1152/jn.00873.2012
- Marshall N et al. (2016). Analysis of power laws, shape collapses, and neural complexity: the NCC toolbox. *Front Physiol* 7:250. https://doi.org/10.3389/fphys.2016.00250
- Paiva ARC, Park I, Príncipe JC (2010). A comparison of binless spike train measures. *Neural Comput Appl* 19:405–419. https://doi.org/10.1007/s00521-009-0307-6
- Poil S-S et al. (2012). Critical-state dynamics of avalanches and oscillations jointly emerge from balanced excitation/inhibition in neuronal networks. *J Neurosci* 32:9817–9823. https://doi.org/10.1523/JNEUROSCI.5990-11.2012
- Wilting J, Priesemann V (2018). Inferring collective dynamical states from widely unobserved systems. *Nat Commun* 9:2325. https://doi.org/10.1038/s41467-018-04725-4

## Funding

- Korea Institute of Toxicology (KIT), Project No. KK-2607
- Korea Ministry of Food and Drug Safety, RS-2024-00331678
- National Research Foundation of Korea (NRF), funded by the Korea government (MSIT), RS-2026-25514372
- Research Institute for Convergence of Biomedical Science and Technology, Pusan National University Yangsan Hospital, 30-2026-005
- National Research Foundation of Korea (NRF), funded by the Korea government (MSIT), RS-2026-25520892
- University of Science and Technology (UST) Korea, 2026 2nd International Presentation Support Program

## Acknowledgement

The developers of CASCADE thank FinalSpark, especially Dr. Ewelina Kurtys, for providing the FinalSpark MEA recording sample.

## Contributing

New manufacturer loaders are welcome. Open an issue or pull request on GitHub.

## Licence

Apache 2.0
