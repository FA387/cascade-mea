# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/spec/v2.0.0.html).

## [0.9.5] - 2026-11-12

Interface and documentation release. All TSV outputs are byte-identical to 0.9.4; only header wording in `criticality_report_<well>.txt` changed.

### Added
- `--citation` in all three tools (text and BibTeX). `CITATION.cff` gains `preferred-citation`; `cascade_run.log` records the version.
- A wrong `--manufacturer` now prints a one-line error instead of a traceback, and names the right device when the file identifies it (`detect_manufacturer`, `read_error_message` in `cascade/loaders/__init__.py`).
- Warning when only one of `--epochA` / `--epochB` is given. Previously the comparison was skipped silently.
- Rian Kang added as author (LICENSE, `CITATION.cff`, `__author__`).

### Changed
- `--output` is a folder in all three tools.
  - `cascade-ltpltd --out` → `--output`; `--out` still works.
  - **Breaking:** `cascade-info --output DIR` writes `DIR/info_summary.tsv`, for a file or a folder. `--output scan.tsv` now creates a folder of that name.
- `cascade-info --folder` folded into `--input`, which takes a file or a folder and scans only the manufacturer's file type, so Axion `.csv` folders now work. `--folder` is kept as a hidden alias.
- Flags grouped by workflow step in `--help` and the README: 1 Setup, 2 Preprocessing, 3 Analysis, 4 Output, 5 Misc. Names, defaults and behaviour unchanged.
- Cortical Labs is `cl` everywhere: `cl1_loader.py` / `CL1Loader` → `cl_loader.py` / `CLLoader`.
- `--legacy_powerlaw` (never released) → `--simple_PL`, pairing with `--simple_BR`.
- Method wording now credits the Habibollahi et al. reference **code** where it differs from their paper (`--perc 0.25`, `--pl_exclude`, β_fit range). No numbers changed. Citation year 2022 → 2023.
- README rewritten: shorter, single manufacturer table, *What's new since the article* and *Effect on the published article*.

### Removed
- `-h`; use `--help`.

### Fixed
- `environment.yml` pip-installed 0.8.0 from PyPI. It now installs only dependencies; CASCADE comes from `pip install .`.
- Version string read from `cascade.__version__` instead of four hard-coded copies.
- `--help`, `--version` and `--citation` no longer create `mea_out/cascade_run.log`.
- Help text: `--well` (Axion merges all wells), `--max_amplitude` (Maxwell single-peak), `--input` (lists `.csv`), `rates_*.tsv` not `.csv`.

## [0.9.4] - 2026-11-12

### Changed
- Branching ratio is now the multi-lag MR estimator of Wilting & Priesemann (2018). v0.8.0 to 0.9.3 used a single-lag regression, which is biased toward 0 under subsampling and shifted when dt was rounded. The new code matches `mrestimator` (r_k to 1e-16, m to 6 decimals) and adds no dependency. **BR values change**; all other outputs are identical.

  | Sample | Old BR | MR m | Checks |
  |---|---:|---:|---|
  | 3Brain | 0.534 | 0.505 | pass |
  | Axion | 0.962 | 0.848 | pass |
  | CL1 | 0.409 | 0.264 | pass |
  | MCS | 0.658 | 0.972 | offset fit disagrees |
  | MaxWell | −0.0004 | 0.988 | offset fit disagrees |

### Added
- Report prints b, autocorrelation time and three checks (r_1 > 0, exp vs exp+offset fit, window ≥ 3τ), each also a `[WARN]` on failure.
- `branching_mr_<well>.tsv` / `.png`.
- `--simple_BR` (old value) and `--br_tmax` (longest lag, default 2000 ms).
- Self-test 12: branching process at 5% subsampling. 14 tests in total.

## [0.9.3] - 2026-11-12

### Added
- `--pl_exclude`: opt-in upper-tail trimming (`EXCLUDE` step of the Habibollahi et al. Python code). With `--simple_PL` it matches their code exactly on all five samples. The trimmed range can shrink to a few values, so the untrimmed fit is always reported. It cannot be combined with `--pl_xmin_S/D`. Default output is unchanged.

## [0.9.2] - unreleased

### Added
- Synchrony metrics (`core/synchrony.py`) for all five manufacturers, from the Axion Neural Metric Tool set: correlogram area and FWHH (raw and normalised), Kreuz SPIKE-distance, Synchrony Index.
  - Default `--sync_metrics core`: `area_xcorr_norm` and `synchrony_index`, the two that resolved on all five samples and separated them best.
  - Flags: `--sync_metrics`, `--sync_window`, `--sync_bin`, `--sync_max_lag`, `--sync_sigma`, `--sync_max_channels`.
  - Outputs: `synchrony_<well>.tsv`, `crosscorrelogram_<well>.png`, `crosscorrelogram_norm_<well>.png`.
  - Timing parameters are absolute ms, never inferred from the sampling rate, so values compare across devices.
  - SPIKE-distance uses PySpike when it is installed (`.[synchrony]`) and a local port otherwise (up to 2.7e-2 off); `kreuz_backend` records which one ran.
  - `fwhh_resolution_limited` flags an FWHH limited by the bin or the lag window.

### Fixed
- `--pl_xmin_S/D` now fix the fit, not only the plot line. The automatic fit is still reported. Runs that used these flags get different τ/α/DCC/SCE.
- `cascade_run.log` said `FAILED` after a clean `sys.exit(0)`.

## [0.9.1] - unreleased

### Fixed
- Power-law normalisation. `tplfit()` summed over the observed avalanche values (`np.unique`) instead of every integer between the cut-offs, which biased the exponent and the KS xmin search. The code came from the Habibollahi et al. reference `tplfit.py`; `plmle.py` in the same package uses the full range. **τ, α, β_pred, DCC and SCE change.** On synthetic data with a true exponent of 1.5, the old fit gave 1.295 at 421 avalanches and 1.499 at 27,288; the corrected fit gives 1.50 at both.

### Added
- `--simple_PL` reproduces the old normalisation byte-for-byte. Self-test covers both paths.

## [0.9.0] - 2026-10-01

### Added
- `--max_amplitude`: removes channels with any spike above N µV peak-to-peak. Off by default; runs after `--ae`, `--max_rate`, `--max_rpv`.
- Per-spike amplitudes in all five loaders (`MEARecording.spike_amplitudes`). Maxwell stores a single-peak value; files without amplitudes skip the filter with a warning.

### Fixed
- Single-channel recordings now skip connectivity. The old guard counted DataFrame columns, so it never fired.
- MCS amplitudes were ADC units, not µV (~16.8× too large); now converted with `ConversionFactor` × 10^`Exponent`.
- MCS and Axion spike times and amplitudes now stay aligned after sorting.
- `CITATION.cff` release date, and the missing `conda-recipe/meta.yaml`.

## [0.8.0] - 2026-04-21

Peer-reviewed release (Avila et al. 2026, *Bioinformatics* 42(8):btag570). Five loaders, spike filtering, bursts, avalanches, criticality (τ, α, DCC, BR, shape collapse), connectivity, spatial entropy, epoch comparison, LTP/LTD.

[0.9.5]: https://github.com/FA387/cascade-mea/releases/tag/v0.9.5
[0.9.4]: https://github.com/FA387/cascade-mea/tree/v0.9.4
[0.9.3]: https://github.com/FA387/cascade-mea/tree/v0.9.3
[0.9.2]: https://github.com/FA387/cascade-mea/tree/v0.9.2
[0.8.0]: https://github.com/FA387/cascade-mea/releases/tag/v0.8.0
