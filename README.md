# CASCADE v0.9.5

**Criticality Avalanche Spike Cross-platform Analysis Detection Engine**

<div align="center">
<img width="250" height="250" alt="CASCADE_Logo" src="https://github.com/user-attachments/assets/7a029bf6-93b7-4185-8859-888672b4bb65" />
</div>

CASCADE analyses multi-electrode array (MEA) spike recordings from five manufacturers with one command. Each file is converted into a common format, so firing rate, network bursts, neuronal avalanches, criticality, functional connectivity, synchrony and spatial entropy are computed the same way for every device. Figures are saved as PNG (300 dpi) and numbers as TSV.

Published in *Bioinformatics* (2026), [doi:10.1093/bioinformatics/btag570](https://doi.org/10.1093/bioinformatics/btag570). See [Citing CASCADE](#citing-cascade).

> **This is v0.9.5.** The article describes v0.8.0. The pipeline is the same, with new features and two corrected calculations. See [What's new since the article](#whats-new-since-the-article).

---

## What's new since the article

**Added**

- [Synchrony metrics](#synchrony-metrics) for all five manufacturers (`--sync_metrics`).
- `--max_amplitude`: removes artefact channels by spike amplitude. See [Noise filtering](#noise-filtering).
- `--pl_exclude`: opt-in upper-tail trimming from the Habibollahi et al. (2023) reference code.
- `--citation` in all three tools; the run log records the CASCADE version.
- `--input` / `--output` in all three tools; flags grouped by workflow step in `--help`.

**Corrected** (a default run can give different criticality values than v0.8.0)

- Power-law fit normalisation: τ, α, β_pred, DCC, SCE.
- Branching ratio: now the Wilting & Priesemann (2018) estimator cited in the article's Methods.

**Fixed:** `--pl_xmin_S/D` now set the fit, not only the plot; MCS amplitudes are in µV; single-channel recordings skip connectivity.

The article's conclusions are not affected, since it reports no τ, α, DCC, BR or SCE values in its text or tables. See [Effect on the published article](#effect-on-the-published-article). Full history: [CHANGELOG.md](CHANGELOG.md).

---

## Tools

| Command | What it does |
|---|---|
| `cascade-info` | Inspect a recording or a folder: duration, channels, active electrodes |
| `cascade` | Full pipeline |
| `cascade-ltpltd` | Classify electrodes as LTP/LTD from firing-rate change between two `cascade` runs |

## Supported manufacturers

`.h5` is used by three devices, so check which one recorded the file. A wrong choice stops with an error that names the right device when the file identifies it.

| Manufacturer | `--manufacturer` | File | What is read |
|---|---|---|---|
| 3Brain BrainWave | `3brain` | `.bxr` | SpikeTimes, SpikeChIdxs, SpikeUnits, SpikeForms |
| Maxwell Biosystems | `maxwell` | `.h5` | `data_store/<stream>/spikes` (frameno, channel, amplitude) |
| Cortical Labs | `cl` | `.h5` | `/spikes` table (timestamp, channel, waveform) |
| Multi Channel Systems | `mcs` | `.h5` | MCS HDF5 v3 SegmentStream spike cutouts |
| Axion BioSystems | `axion` | `.csv` | AxIS spike-list export |

---

## Installation

```bash
git clone https://github.com/FA387/cascade-mea.git
cd cascade-mea
conda env create -f environment.yml
conda activate cascade-mea
pip install .
```

Check it:

```bash
cascade --test       # 14 self-tests on synthetic data; ends with "All 14 tests passed."
cascade --version
```

`pip install cascade-mea` installs the latest PyPI release. If `cascade --version` is older than this README, install from GitHub as above.

Sample recordings from the article: https://doi.org/10.5281/zenodo.19161615

---

## Quick start

Replace `[manufacturer]` with a value from [Supported manufacturers](#supported-manufacturers).

```bash
# 1. Inspect a recording (or a whole folder; --output saves info_summary.tsv)
cascade-info --manufacturer [manufacturer] --input [recording] --ae 5
cascade-info --manufacturer [manufacturer] --input [folder] --ae 5 --output [output_dir]

# 2. Full analysis
cascade --manufacturer [manufacturer] --input [recording] --output [output_dir]
cascade --manufacturer [manufacturer] --input [recording] --ae 5 --perc 0.25 --output [output_dir]   # article settings
cascade --manufacturer axion --input [recording] --well B1 --output [output_dir]                    # one Axion well
cascade --manufacturer [manufacturer] --input [recording] --epochA 0,60 --epochB 120,60 --output [output_dir]

# 3. LTP/LTD between runs (one or more --measure files)
cascade-ltpltd --baseline [baseline]/rates_[well].tsv --measure [stim]/rates_[well].tsv --output [output_dir]
```

Without `--well`, Axion merges all wells into one network; analyse wells separately for criticality.

---

## Output files

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

### cascade-info output

Printed to the terminal. With `--output`, also saved as `[output_dir]/info_summary.tsv`, one row per file.

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

Grouped by workflow step, in the same order as `--help`:

| Step | Decides |
|---|---|
| **1. Setup** | Which file, which device, where results go |
| **2. Preprocessing** | Which part of the recording and which channels to trust |
| **3. Analysis** | Parameters of each calculation |
| **4. Output** | Extra figure formats and terminal output |
| **5. Misc** | Information only, no analysis |

Most runs need only step 1. Add step 2 filters for noisy recordings, change step 3 only with a reason, and report any non-default value in your methods.

### cascade

#### 1. Setup

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `maxwell`, `cl`, `mcs`, `axion` |
| `--input` | auto | Recording file. If omitted, the first `.bxr`/`.h5` in the current folder |
| `--well` | first | Well or stream. Axion: all wells merged unless given |
| `--output` | `mea_out` | Output folder |

#### 2. Preprocessing

| Flag | Default | Description |
|---|---|---|
| `--start` | `0` | Window start (s) |
| `--dur` | full | Window duration (s) |
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
| `--test` | Run the 14 self-tests |

### cascade-info

#### 1. Setup

| Flag | Default | Description |
|---|---|---|
| `--manufacturer` | required | `3brain`, `maxwell`, `cl`, `mcs`, `axion` |
| `--input` | required | Recording file, or a folder (scans that manufacturer's file type) |
| `--well` | first | Well or stream. Axion: all wells merged unless given |
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

Filters remove whole channels, in this order: `--ae`, then `--max_rate` and `--max_rpv`, then `--max_amplitude`.

```bash
cascade --manufacturer [manufacturer] --input [recording] --output [output_dir] \
        --ae 5 --max_rate 1000 --max_rpv 1.0 --max_amplitude 200 --verbose
```

Add only the filters a first run shows you need, check how many channels each removed (`--verbose`), and use the same settings for every recording you compare. Amplitudes are in µV for all devices, but Maxwell stores a single-peak value, so calibrate the threshold per device.

---

## Synchrony metrics

Spike-timing synchrony from the Axion Neural Metric Tool set, computed for every manufacturer. Unlike Pearson connectivity (1 s bins), it resolves millisecond timing.

| `--sync_metrics` | Reports |
|---|---|
| `core` (default) | `area_xcorr_norm` (normalised correlogram area) and `synchrony_index` (Paiva et al. 2010) |
| `xcorr` | Raw area and both full widths at half height (FWHH) |
| `kreuz` | `kreuz_metric` = 1 − SPIKE-distance (Kreuz et al. 2013) |
| `index` | `synchrony_index` only |
| `all` | Everything |
| `none` | Skip |

`core` keeps the two metrics that resolved on all five sample recordings and separated them best. `kreuz_metric` separates cultures poorly but is parameter-free, so it suits before/after comparisons within one culture.

Notes:

- Timing flags are absolute ms, never inferred from the sampling rate, so values are comparable across devices.
- FWHH cannot resolve a peak narrower than `--sync_bin`. Halve the bin until FWHH stops changing, or compare the area metrics instead.
- For the reference SPIKE-distance, install PySpike (`pip install ".[synchrony]"`). Without it CASCADE uses its own port (up to 2.7e-2 off); `kreuz_backend` in the TSV says which was used.
- `kreuz` and `index` are O(N²) in channels and capped at `--sync_max_channels`; subsampling is recorded in the TSV.

---

## Differences from the published article

| | v0.8.0 (article) | v0.9.5 default | Old behaviour |
|---|---|---|---|
| **Power-law normalisation** (τ, α, β_pred, DCC, SCE) | Sum over observed avalanche values only; biased below ~1,000 avalanches | Sum over every integer between the cut-offs (Clauset et al. 2009) | `--simple_PL` |
| **Branching ratio** | Single-lag regression, biased toward 0 under subsampling | Multi-lag MR estimator (Wilting & Priesemann 2018), validated against `mrestimator` | `--simple_BR` |
| **Upper-tail trimming** | None | None; opt-in `--pl_exclude` | |

The old normalisation came from the Habibollahi et al. (2023) reference `tplfit.py`. On synthetic data with a true exponent of 1.5, it gave 1.295 at 421 avalanches and 1.499 at 27,288, while the corrected fit gives 1.50 at both.

### Reproducing the published v0.8.0 results

```bash
cascade --manufacturer [manufacturer] --input [recording] --ae 5 --perc 0.25 --simple_PL --simple_BR --output [output_dir]
```

Both flags print a warning. Use them only to regenerate old numbers.

### Effect on the published article

| Article item | Affected? |
|---|---|
| Table 1: active electrodes, spikes, avalanches | No. Filtering and avalanche detection are unchanged |
| Table 1: run time, peak memory | Longer in v0.9.5 because of new modules, not because of the corrections |
| Fig. 1: rasters vs proprietary software | No |
| Methods: BR "following Wilting and Priesemann (2018)" | v0.9.4+ implements this estimator, so the software now matches the text |
| Supplementary Fig. 1C–D (CL1 game, 27,288 avalanches) | τ 1.752 → 1.770, α 1.974 → 1.990 (about 1%) |
| Supplementary Table S1 | No |

Two clarifications on the article's wording: `cascade-ltpltd` is Python (not R), ships with CASCADE since v0.8.0, and classifies electrodes by firing-rate change. `--perc 0.25` is the reference-code default; the Habibollahi et al. (2023) paper reports 40% of the median.

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
└── tests/self_test.py      the 14 checks run by cascade --test
```

## Adding a new manufacturer

1. Add `cascade/loaders/<name>_loader.py` subclassing `BaseLoader`.
2. Implement `load(path, well=None, start_s=0.0, duration_s=None) -> MEARecording` with spike times in **seconds**.
3. Register it in `run_analysis.py:get_loader()` and the `--manufacturer` choices.

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

Korea Institute of Toxicology (KIT, Project No. KK-2607) and the Korea Ministry of Food and Drug Safety (RS-2024-00331678).

## Contributing

New manufacturer loaders are welcome. Open an issue or pull request on GitHub.

## Licence

Apache 2.0
