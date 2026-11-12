# Changelog

All notable changes to CASCADE are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
CASCADE adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.9.5] — 2026-11-12

Documentation and interface release: no analysis result changes. Every TSV output is
byte-identical to 0.9.4; only the header lines of `criticality_report_<well>.txt`
changed wording.

### Changed — `cascade-ltpltd --out` is now `--output`

All three tools now name their output `--output`. The old `--out` still works
(argparse accepts it as an abbreviation of `--output`), so existing scripts are
unaffected. Its help text also now says `rates_*.tsv`, which is what `cascade`
writes, instead of `rates_*.csv`.

### Changed — `cascade-info`: `--input` takes a folder, `--output` is a folder

`--output` now means an output folder in all three tools.

- `--input` accepts a recording file **or a folder**. A folder is scanned for the
  `--manufacturer`'s own file type (3brain `.bxr`; maxwell/cl/mcs `.h5`; axion `.csv`),
  so Axion folders can now be scanned and other manufacturers' files are skipped.
  `--folder` is no longer documented; it still works as a hidden alias of `--input`.
- `--output DIR` saves the summary table as `DIR/info_summary.tsv`, for a single file
  as well as a folder. Previously it took a file path, worked in folder mode only, and
  was silently ignored with `--input`. **Scripts that passed `--output scan.tsv` now get
  a folder named `scan.tsv`** containing `info_summary.tsv`.
- Without `--output` nothing is written, as before.
- `cascade --input` help now lists `.csv` and says the auto-search uses the first
  `.bxr`/`.h5` file in the current folder.

### Changed — flags grouped by workflow step

`--help` of all three tools and the README flag tables now group the flags by step:
**1 Setup** (input, output, manufacturer), **2 Preprocessing** (time window, spike
filters), **3 Analysis** (parameters of each calculation), **4 Output** (figure formats,
terminal), **5 Misc** (`--help`, `--version`, `--citation`, `--test`). `--bin` moved out
of "Spike filtering" (it sets the firing-rate bin, it filters nothing) and
`--simple_PL`/`--simple_BR` share a "reproducing v0.8.0" group. Layout only: every
flag's name, default and behaviour is unchanged. The short `-h` is removed; use `--help`.

### Fixed — `environment.yml` installed the old release

`environment.yml` also pip-installed `cascade-mea` from PyPI, which serves 0.8.0, so
`conda env create` alone gave the old methods (the README's following `pip install .`
replaced it). The PyPI line is removed; CASCADE is installed from the folder with
`pip install .` as the README says.

### Changed — clear message for a file the chosen manufacturer can't read

A wrong `--manufacturer` (easy with `.h5`, which three devices share) used to end in a
Python traceback. `cascade` and `cascade-info` now print one plain message with the
original error and, when the file's own attributes identify the device, the right value:
`This file looks like a 'cl' recording - try --manufacturer cl.` Detection lives in
`cascade/loaders/__init__.py` (`detect_manufacturer`, `read_error_message`), shared by both
tools. Results for readable files are unchanged.

### Changed — Cortical Labs is `cl` everywhere

The loader is now `cascade/loaders/cl_loader.py` with class `CLLoader` (was
`cl1_loader.py` / `CL1Loader`), and it reports itself as "Cortical Labs", matching the
`--manufacturer cl` value. Code that imported `CL1Loader` must use `CLLoader`. Output
files are byte-identical. (Detection still matches the `CL1-…` `system_id` the device
writes into its files.)

### Fixed — help text and a silent skip

- `--epochA` without `--epochB` (or the reverse) was silently ignored; it now prints a
  warning that the epoch comparison was skipped.
- `--well` help (cascade, cascade-info) now says that Axion merges all wells when no
  well is given; `--max_amplitude` help notes that Maxwell stores a single-peak value.

### Fixed — README

- `cascade-ltpltd` examples now work (see above).
- `--manufacturer` value for Cortical Labs is `cl`, not `cl1`.
- `--sync_metrics` default is `core` (normalised correlogram area + Synchrony
  Index), not `xcorr`.
- β_fit range is xmin_D < D < xmin_D + 60 (fitted duration cut-off), not
  `[t_min, t_min + 60]`.
- New sections: *What's new since the article* (near the top), *Citing CASCADE*,
  *Reproducing the published v0.8.0 results* (`--simple_PL --simple_BR`),
  *Differences from the published article*.
- The two manufacturer tables are merged into one under *Supported manufacturers*,
  now including the `--manufacturer` value.
- Logo restored; self-test output, output-file list (`wave_templates_<well>.png`),
  report and run-log descriptions, project structure and the `inspect_file`
  command brought up to date; Clauset 2009, Marshall 2016, Kreuz 2013 and
  Paiva 2010 added to the references.
- `--test`, `--version` and `--citation` documented in the flag tables.

### Added — `--citation`

- `cascade`, `cascade-info` and `cascade-ltpltd` accept `--citation`: prints the
  published article (Avila et al. 2026, *Bioinformatics* 42(8):btag570,
  doi:10.1093/bioinformatics/btag570) as text and BibTeX, plus the running version.
- `CITATION.cff` gains a `preferred-citation` block for the article.
- `cascade_run.log` now records the CASCADE version that produced the results.
- Rian Kang added as an author (LICENSE, `CITATION.cff`, `__author__`), matching the article.

### Fixed — version string had four hand-typed copies

`--version` in all three tools and the run banner now read `cascade.__version__`
instead of a hard-coded "v0.9.4", so they cannot drift at the next bump.
`cascade --help`, `--version` and `--citation` no longer create
`mea_out/cascade_run.log` in the current folder.

### Changed — `--legacy_powerlaw` renamed to `--simple_PL`

Pairs with `--simple_BR`: each flag switches one metric back to the simpler v0.8.0
method (`--simple_PL` the power-law normalisation, `--simple_BR` the branching ratio).
Use both to reproduce v0.8.0 exactly. The old name was never in a released version,
so there is no alias; CHANGELOG entries below use the new name.

### Fixed — method attribution wording (no numerical change)

The report, CLI help and docstrings said "Habibollahi 2022 method/default" where
CASCADE actually follows the Habibollahi et al. reference **code**, which differs from
their published paper in places. Wording only: every number, default, flag and
output file is unchanged. The one visible change in outputs is the power-law header
line of `criticality_report_<well>.txt`.

- Report header: "truncated MLE, Habibollahi 2022 method" → "discrete truncated
  power-law MLE, xmin by KS minimisation".
- `--perc 0.25` is the default of the reference code (`AV_analysis_BurstT`); the
  paper reports a threshold of 40 % of the median.
- `--burst_min 12` / `--t_min 4` are CASCADE's defaults; the reference example uses
  bm = 10, tm = 2.
- β_fit range is t_min < D < t_min + 60, strict at both ends (Python reference; the
  MATLAB original uses + 30), not `[t_min, t_min + 60]`.
- `--simple_PL` reproduces CASCADE ≤ 0.9.0 exactly, with the normalisation of
  the reference `tplfit.py`. It does not reproduce `tplfit.py` exactly, because
  CASCADE clips the fitted CDF to [0, 1] and the reference does not.
- `--pl_exclude` reproduces the **Python** reference
  (`AV_analysis_ExponentErrorComments.py`); the MATLAB original differs (final fit
  with `plmle`, no duration pre-cut).
- Citation year 2022 → 2023 (Nat. Commun. 14:5287), matching the README and the
  CASCADE paper.

---

## [0.9.4] — 2026-11-12

Replaces the branching ratio with the estimator the CASCADE paper cites. **BR values
change.** No other output changes: every TSV, and every line of the criticality
report outside the branching-ratio block, is identical to v0.9.3.

### Changed — branching ratio is now the multi-lag (MR) estimator of Wilting & Priesemann (2018)

v0.8.0–v0.9.3 labelled BR "Wilting & Priesemann 2018" but computed a simple
single-lag regression of mean N(t+1) on N(t), one point per distinct N(t),
unweighted. That is the conventional estimator Wilting & Priesemann show is biased
toward 0 under subsampling — and an MEA always subsamples the network. It was also
unstable: rounding dt to 1 µs moved 3Brain from 0.534 to 0.607.

BR is now m from the MR estimator, as described in Habibollahi et al. (2023, Eq. 7):
regress a(t+k) on a(t) for k = 1 … kmax, fit r_k = b·m^k. Validated against the
`mrestimator` toolbox (Spitzner et al. 2021) on all five sample recordings: r_k agree
to 1e-16 and m, b and the offset-fit m to 6 decimals. No new dependency.

| sample | v0.9.3 BR (simple) | v0.9.4 BR (MR m) | b | τ (ms) | checks |
|---|---:|---:|---:|---:|---|
| 3Brain | 0.534 | 0.505 | 1.19 | 9.1 | pass |
| Axion | 0.962 | 0.848 | 1.17 | 89.4 | pass |
| CL1 | 0.409 | 0.264 | 1.74 | 37.6 | pass |
| MCS | 0.658 | 0.972 | 1.08 | 481 | offset fit disagrees |
| MaxWell | −0.0004 | 0.988 | 0.04 | 302 | offset fit disagrees |

The report also gives b (subsampling factor), τ = −dt/ln m (autocorrelation time,
independent of dt), the lag window, and three checks, each also printed as `[WARN]`
when it fails:
- **r_1 > 0**: lag-1 slope significant against uncorrelated activity
  (z = r_1·√(T−1), p < 0.05); fails on white noise.
- **exponential vs exponential+offset**: τ within ±50 %. They disagree under slow
  drift or recurring bursts, as on MCS and MaxWell, where r_k turns negative at long
  lags (see the plot). Treat m ≈ 1 there with caution.
- **window ≥ 3 τ**: the lag window must span the decay; raise `--br_tmax` if not.

New outputs: `branching_mr_<well>.tsv` (r_k per lag with both fits) and
`branching_mr_<well>.png`. `branching_points_<well>.tsv` is still written.

New flags:
- `--simple_BR` reports the old single-lag value. With `--simple_PL` it
  reproduces v0.8.0 numbers exactly. Use only for previously published results.
- `--br_tmax` (default 2000 ms) sets the longest lag.

Self-test: new test 12 simulates a branching process (m = 0.9) recorded at 5 %
subsampling; MR must recover m within 0.05 and the simple estimate must stay below
0.5 (it gives ≈ 0.2). 14 tests.

---

## [0.9.3] — 2026-11-12

Adds Habibollahi's upper-tail trimming as an opt-in flag. No default output changes: without `--pl_exclude`, every output is identical to v0.9.2.

### Added — `--pl_exclude`: Habibollahi's upper-tail trimming, opt-in

CASCADE has always copied only the xmin search of Habibollahi et al. (2022)
(`tplfit`), fitting up to the largest avalanche. Their pipeline first trims the
upper tail with `EXCLUDE()`: pre-cut at max^0.8, then drop the largest value until
KS < min(1/√n, 0.1), cap duration xmax at max^0.9 when it ends below 40, and fit
on [xmin, xmax]. `--pl_exclude` adds that step. Off by default; nothing changes
without it.

With `--simple_PL` it reproduces their code exactly: τ, α, both fit ranges,
β_pred, β_fit and DCC match `AV_analysis_ExponentErrorComments` to 0.0 on all five
sample recordings, at bm/tm = 12/4 and 20/3. That needed `tplfit(clip_cdf=False)`:
CASCADE has clipped the fitted CDF to [0, 1] since 0.8.0, but the legacy CDF can
exceed 1 and the reference does not clip, so KS and the trimming point differed
(MCS stopped at 180 instead of 132). Clipping stays the default, so 0.8.0 output
is still reproduced byte-for-byte.

Caveat, stated in the help text too: the stopping rule chases KS, so the trimmed
range passes the fit by construction and can shrink to very few values (3Brain
durations [4, 6], 38 of 447 avalanches). It moves DCC in both directions on the
samples (3Brain 8.61 → 0.35, MaxWell 4.54 → 21.99). The untrimmed fit is always
written to the report, with the fit ranges and avalanches used. Cannot be combined
with `--pl_xmin_S/D`.

---

## [0.9.2] — unreleased

Adds spike-timing synchrony analysis. No existing metric changes: every output
produced by v0.9.1 is unaffected, and the new metrics are written to their own
files.

### Fixed — `--pl_xmin_S` / `--pl_xmin_D` now fix the fit, not just the plot

Since v0.8.0 these flags only moved where the log-log fit line started; τ, α,
β_pred, DCC and SCE still came from the automatically selected xmin, although the
help text said "override the xmin for the fit". The plot could therefore draw a τ
fitted from one xmin starting at another. The flags now fix the lower cut-off of
the fit itself (`tplfit(..., xmin=k)`, the `bm` floor of Habibollahi's
`EX_burst=0` path), and every downstream value and the plot use it.

A hand-picked xmin is a free parameter, and near τ ≈ 1 it moves DCC a lot (3Brain
sample: DCC ≈ 1.0 at xmin_S = 2 vs 11.1 at the automatic 4). So the automatic fit is
always computed too, printed as a `[WARN]` and written to the criticality report
beside the manual one. Runs without the flags are unchanged byte-for-byte. Runs that
passed them before now get different τ/α/DCC/SCE. The flags now take integers.
Self-test: the power-law test also asserts that a fixed xmin reproduces the search's
fit at that xmin exactly.

### Fixed — run log said `FAILED` after a clean exit

`cascade --test` and `--help` end with `sys.exit(0)`, which `_run_with_log` counted
as a failure, so `cascade_run.log` read `Result: FAILED` under "All 13 tests
passed". Exit code 0 / None is now logged as `SUCCESS`.

### Added — synchrony metrics (new `core/synchrony.py`)

The synchrony metric set of the Axion Neural Metric Tool, applied uniformly to
**all five** supported manufacturers rather than only Axion. These are spike-timing
measures and are distinct from the existing Pearson connectivity, which correlates
firing rates in 1 s bins and therefore cannot resolve millisecond structure.

| Metric | Output key | Source |
|---|---|---|
| Area under cross-correlation | `area_xcorr` | Axion NMT v4.1 §4.4 |
| Area under normalised cross-correlation | `area_xcorr_norm` | Axion NMT v4.1 §4.4 |
| Full width at half height | `fwhh_xcorr` | Axion NMT v4.1 §4.4 |
| FWHH of normalised cross-correlation | `fwhh_xcorr_norm` | Axion NMT v4.1 §4.4 |
| Kreuz metric = 1 − SPIKE-distance | `kreuz_metric` | Kreuz et al. 2013, Eq 19 |
| Synchrony Index | `synchrony_index` | Paiva et al. 2010, Eqs 9/11 |

**Default is `core`: two metrics, not six.** Measured across the five published
sample recordings, only `area_xcorr_norm` (26× spread, CV 0.78) and
`synchrony_index` (11×, CV 0.69) both resolve on every manufacturer *and* separate
recordings. `area_xcorr` raw carries the same information confounded by firing rate;
`fwhh_xcorr_norm` resolves on only 2 of 5 and its apparent 552× spread is largely the
resolution artefact; `kreuz_metric` varies by just 1.2× (CV 0.09) across very
different cultures. The others stay available explicitly. Kreuz is retained because
it is parameter-free, which suits tracking change within a single culture over time
even though it discriminates between cultures poorly.

**FWHH resolution flag.** Full width at half height cannot resolve a peak narrower
than the bin or wider than the lag window, yet still returns a finite-looking number.
The value is reported exactly as computed and never altered; `fwhh_resolution_limited`
and `fwhh_limit_reason` now accompany it, and a `[WARN]` names the flag that would fix
it. Verified: fires on 3Brain (bin-limited) and MCS (lag-truncated at 584 ms vs a
500 ms window), silent on Axion where the value is genuine.

New flags: `--sync_metrics {none,core,xcorr,kreuz,index,all}` (default `core`),
`--sync_window` (ms, default 20), `--sync_bin` (ms, default 1), `--sync_max_lag`
(ms, default 500), `--sync_sigma` (ms, defaults to the window),
`--sync_max_channels` (default 64).

New outputs: `synchrony_<well>.tsv`, `crosscorrelogram_<well>.png` and
`crosscorrelogram_norm_<well>.png`, the latter two with the synchrony window drawn
as grey vertical lines as the Axion display does.

**Implementation notes**

- The **Kreuz SPIKE-distance** defers to **PySpike** when it is installed (new
  optional dependency, `pip install "cascade-mea[synchrony]"`). PySpike is the
  reference implementation by an author of the 2013 paper, and CASCADE then matches
  it to machine precision — verified across identical, offset, rate-mismatched and
  random Poisson trains, worst absolute difference 0.0e+00. A cross-validation
  assertion runs in the self-test whenever PySpike is present.

  Without PySpike, CASCADE uses its own port. That port satisfies every guarantee the
  paper states — exactly 0 for identical trains, time-scale invariant to 5.6e-17,
  monotonic with mismatch — but differs from the reference by up to 2.7e-2 on
  structured trains. The cause is the edge intervals: both implementations agree
  algebraically on the profile (PySpike's `dist_at_t` is Eq 19), but the paper says
  only "add auxiliary leading spikes at t=0 and trailing spikes at t=T", while
  PySpike clamps the edge ISI against the adjacent real ISI in a way the paper does
  not state. Rather than ship a near-miss of a published measure, CASCADE prefers the
  reference and records which path produced each number as `kreuz_backend`.

- The local port uses the *improved* dissimilarity profile of Eq 19,
  not the original Eq 13. The 2013 paper exists precisely because the original gave
  spuriously high instantaneous values for event-like firing. Eq 19 rather than the
  intermediate Eq 18 matters too: Eq 19 adds the cross-train ISI weighting that the
  paper says makes the distance "increase with the difference in firing rate".
  Auxiliary spikes are inserted at t=0 and t=T as the paper requires. S(t) is
  piecewise linear on the merged spike grid, so the time integral is exact and needs
  no sampling.
- The **Synchrony Index** is evaluated with the closed form of Paiva Eq 11, so no
  binning is involved. The double sums are windowed to ±5σ via `searchsorted`,
  which is exact to ~4e-6 of the peak and turns O(N_i·N_j) into roughly O(N_i·k).
- The **pooled correlogram** is computed over all unique pairs in the frequency
  domain using `Σ_{i<j} conj(F_i)·F_j = (|Σ F_i|² − Σ|F_i|²)/2`, giving
  O(N·T log T) instead of O(N²·T log T) — ~58,000 pairwise products avoided on the
  341-active-channel MaxWell recording. The per-channel FFTs are accumulated in a
  streaming loop rather than materialised as one array: holding all of them would
  cost ~1.4 GB for MaxWell at 1 ms bins, several times the pipeline's entire
  reported peak memory.

**Cross-manufacturer design decision.** All synchrony timing parameters are absolute
durations in milliseconds and are deliberately **not** inferred from the sampling
rate or mean ISI, unlike the avalanche bin width which `infer_dt()` adapts to the
data. Sampling rates differ substantially between devices, so a sampling-rate-derived
correlogram bin would make "area under the normalised cross-correlation at a 20 ms
window" a different physical quantity per manufacturer, defeating the cross-platform
comparability that is CASCADE's purpose. This asymmetry with the avalanche pipeline
is intentional and documented in the module so it is not "fixed" later.

**Known limitation — FWHH resolution.** Full width at half height cannot resolve a
correlogram peak narrower than `--sync_bin`. On the sample recordings Axion's FWHH is
converged at the 1 ms default (144.5 / 141.4 / 133.2 / 121.1 ms at bins of
1.0 / 0.5 / 0.2 / 0.1 ms — a real ~130 ms peak), whereas 3Brain and MaxWell track the
bin width (MaxWell ~2× the bin at every step) and are therefore resolution-limited;
at 0.1 ms they have still not converged, putting their synchrony below 100 µs, close
to the ~0.05 ms sampling period. For such recordings FWHH is a lower bound and the
area metrics are the better comparison. Halve `--sync_bin` until FWHH stops changing.

**Validated** against cases with known answers: SPIKE-distance is exactly 0 for
identical trains, grows monotonically with mismatch, and is time-scale invariant to
machine precision under a 3× stretch; the Synchrony Index is exactly 1 for identical
trains and decreases with offset; and the correlogram FWHH reproduces the analytic
prediction 2.355·j·√2 for Gaussian jitter of standard deviation j across a 20× range
of j. Self-test count 12 → 13.

---

## [0.9.1] — unreleased

### Fixed — power-law normalisation (changes published τ, α, β_pred, DCC, SCE)

`tplfit()` normalised the discrete truncated power-law likelihood over only the
**observed unique** avalanche values instead of **every integer** between the
cut-offs:

```python
# before                      # after
s = np.unique(subset)         s = np.arange(x0, xmax + 1)
```

The correct form is `p(y) = y^-a / Σ_{k=xmin}^{xmax} k^-a` (Clauset et al. 2009).
Restricting the sum to observed values makes the normalising constant depend on
the sample, which is not a valid likelihood: the MLE matches observed mean-log to
model mean-log, and dropping unobserved large integers lowers the model's
mean-log, so the exponent falls to compensate. The same `s` was also reused for
`A` in the KS block while `fit_cdf` spanned the full integer grid, so the fitted
CDF did not integrate to 1 and **xmin selection** was corrupted too. One change
fixes both.

**Provenance.** This was inherited, not introduced. `tplfit.py` and `tplfit.m`
from Habibollahi et al. (2022) both use `unique(...)`, and CASCADE replicated
them faithfully. Notably, `plmle.py` in that same reference package — the NCC
toolbox fitter of Marshall et al. (2016) — computes `unqZ = np.unique(z)` and
then **never uses it**, normalising over `allZ = np.arange(xmin, xmax+1)`
instead. CASCADE now follows `plmle`, not `tplfit`.

**Verified** against the real function on synthetic data with known exponent
(true α = 1.5):

| avalanches | legacy | corrected |
|-----------:|-------:|----------:|
| 421 | 1.295 | 1.502 |
| 520 | 1.308 | 1.496 |
| 705 | 1.338 | 1.500 |
| 898 | 1.372 | 1.498 |
| 9,197 | 1.484 | 1.496 |
| 27,288 | 1.499 | 1.500 |

Bias shrinks as counts grow, because the observed values fill in the range. It
also grows with xmin: on the real 3Brain sample the legacy path returns τ < 1 for
xmin ≥ 3, which is physically implausible and is the diagnostic tell.

**Impact on the published sample dataset** — the direction is *not* uniform,
because KS-based xmin selection shifts as well, so every recording must be refit
rather than offset by a constant:

| | 3Brain (421) | Axion (898) | MCS (705) | CL1 (520) |
|---|---:|---:|---:|---:|
| τ legacy | 1.570 | 2.221 | 1.622 | 2.115 |
| τ corrected | 1.193 | 1.768 | 1.889 | 2.536 |
| DCC legacy | 2.561 | 0.655 | 0.036 | 0.300 |
| DCC corrected | 11.127 | 0.456 | 0.157 | 0.319 |
| KS(size) | 0.216 → 0.117 | 0.063 → 0.043 | 0.116 → 0.064 | 0.025 → 0.026 |

KS improves in 3 of 4, independent evidence the correction is right. **3Brain's
corrected DCC of 11.1 should not be read as a measurement**: with 421 avalanches
and KS still 0.117 after the fix, that recording does not support a clean
power-law fit, and τ → 1.19 makes β_pred = (α−1)/(τ−1) diverge as the denominator
approaches zero. The bias was masking a weak fit.

### Added

- **`--simple_PL`** reproduces the pre-0.9.1 normalisation exactly, so
  previously published numbers stay reproducible. Verified: with this flag, all
  four manufacturers reproduce the bundled v0.8.0 reference output byte-for-byte.
  It prints a warning and must not be used for new analyses.
- **Self-test 12 tests** (was 11). The new test fits data drawn from a known
  truncated power law and asserts the corrected path recovers the exponent to
  within 0.10, *and* that `legacy_normalization=True` still reproduces the old
  bias, so the reproducibility path cannot silently rot.

### Action required

- τ, α, β_pred, DCC and SCE in the v0.8.0 manuscript need regenerating. β_fit,
  branching ratio, bursts, connectivity and entropy are unaffected.
- The correction should state that the convention came from the cited reference
  implementation, since that is where it originates.
- Version strings are already set to 0.9.1 throughout (`pyproject.toml`,
  `CITATION.cff`, `cascade/__init__.py`, the four sub-package docstrings, the
  banner, the three `--version` outputs and `conda-recipe/meta.yaml`). Tag the
  release as `v0.9.1` when cutting it.

---

## [0.9.0] — 2026-10-01

First release after the peer review of v0.8.0. This is a **refinement and
addition** release: one new opt-in filter, plus bug fixes. No existing analysis
mathematics was changed.

> **Reproducibility — verified on real recordings.** Every module behind the
> published metrics is byte-identical to v0.8.0: `core/criticality.py`,
> `core/avalanche.py`, `core/connectivity.py`, `core/bursts.py`, `core/entropy.py`
> and `core/plots.py`.
>
> This was confirmed empirically against the published sample dataset
> ([10.5281/zenodo.19161615](https://doi.org/10.5281/zenodo.19161615)) and its
> bundled v0.8.0 reference output. Re-running the documented command line
>
> ```
> --ae 5 --perc 0.25 --ava_size 2 --ava_dur 1
> ```
>
> on v0.9.0 reproduced **every** `.tsv` and the `criticality_report` byte-for-byte
> for 3Brain and Axion, including τ, α, DCC, branching ratio and shape collapse —
> on a different OS and kernel from the one that produced the reference. The new
> filter is off by default, so existing command lines are unaffected.

### Added

- **`--max_amplitude` — amplitude-based artefact rejection.** Removes any channel
  whose largest single spike exceeds the given peak-to-peak amplitude in µV.
  Electrodes recording mechanical knocks or electrical interference produce spikes
  far larger than real action potentials, and those channels previously had to be
  excluded by hand. Default `0.0` (off); recommended starting value `200`.

  ```bash
  cascade --manufacturer 3brain --input rec.bxr --output mea_out \
          --ae 5 --max_rate 1000 --max_rpv 1.0 --max_amplitude 200
  ```

  It is the fourth spike-level filter and runs last, after `--ae`, `--max_rate`
  and `--max_rpv`. With `--verbose` the removed channel IDs are listed.

- **`filter_amplitude_electrodes()`** in `core/spikes.py`. Returns
  `(filtered, n_removed, removed_ch)`, matching the existing filter functions.
  A threshold of `0`, or a recording without amplitude data, is a no-op.

- **`MEARecording.spike_amplitudes`** — new optional field holding
  `{channel_id: array of per-spike peak-to-peak amplitudes in µV}`. Defaults to
  `None`, so existing code constructing an `MEARecording` is unaffected.

- **Per-spike amplitude extraction in all five loaders.** v0.8.0 read spike
  *times* only and discarded the waveform information the files already carried.

  | Loader  | Source | Conversion to µV |
  |---------|--------|------------------|
  | CL1     | `spikes["samples"]` (75 samples/spike) | already µV |
  | 3Brain  | `<well>/SpikeForms`, reshaped by the `Wavelength` attribute | already µV |
  | Maxwell | `amplitude` field of the spike table | `abs()` of the stored value |
  | MCS     | `SegmentData_<N>` waveform cutouts | `ConversionFactor × 10^Exponent × 1e6` |
  | Axion   | amplitude column of the spike-list CSV | `abs()`, mV → µV |

  Where a file carries no waveform data — a Maxwell export without the optional
  amplitude field, a 3Brain file without `SpikeForms` — `spike_amplitudes` is
  `None`, and `--max_amplitude` prints a warning and continues instead of failing.

- **Self-test coverage for the new filter.** `cascade --test` now runs 12 tests
  (was 10; see also the power-law regression test above). The synthetic recording carries per-spike amplitudes around 50 µV with
  three channels seeded with a 500 µV artefact spike, and the new test asserts that
  a 200 µV threshold removes exactly those three, that `0` disables the filter, and
  that a recording without amplitudes passes through untouched.

- **`CLAUDE.md`** — repository guidance for AI coding assistants.

### Fixed

- **Connectivity was never correctly skipped for single-channel recordings.**
  The guard in `run_analysis.py` tested `rate_df.shape[1] < 2`, which counts
  DataFrame *columns*. `rate_df` is long-format with a fixed schema, so that value
  is a constant ≥ 2 and the condition could never be true. A recording with one
  active channel therefore passed the guard and failed inside the Pearson
  correlation instead of reporting the intended warning. The guard now tests
  `rate_df.empty or rate_df["ch_idx"].nunique() < 2`. The correlation mathematics
  in `core/connectivity.py` is unchanged.

- **MCS amplitudes were expressed in ADC units, not µV.** The MCS loader used raw
  integer waveform values, so a `--max_amplitude 200` meant something different for
  MCS files than for the other four manufacturers. Measured on the published MCS
  sample, the as-shipped code produced values of 98–2733 (median 775) where the
  true amplitudes are 5.8–162.9 µV (median 46.2) — a scale factor of ~16.8, i.e.
  ~0.0595 µV per ADC step. `--max_amplitude 200` would therefore have discarded
  **57 of 58 channels**, leaving a single channel, which would then have hit the
  broken single-channel connectivity guard fixed below. The two defects compounded. The loader now converts using
  `ConversionFactor` and `Exponent` from `InfoChannel`
  (`value_V = raw × ConversionFactor × 10^Exponent`, then ×1e6 for µV), and records
  `metadata["amplitude_unit"]`. Peak-to-peak is a difference, so the `ADZero`
  offset cancels and is not required. If a file lacks these fields, or reports a
  unit other than volts, the loader returns no amplitudes rather than raw counts —
  a µV threshold is never silently compared against ADC units.

- **MCS lost time/amplitude alignment across segments.** When one source channel
  appeared in more than one `SegmentStream` segment, the merged spike times were
  re-sorted while the amplitudes were only concatenated, so spike *i* no longer
  paired with amplitude *i*. Both arrays now carry the same sort permutation. This
  did not affect `--max_amplitude`, which takes a per-channel maximum and is
  order-independent, but it would have corrupted any per-spike use of amplitudes.

- **Axion per-channel sorting** now uses an explicit `argsort` index rather than
  `np.sort`, keeping spike times and amplitudes aligned. Output ordering is
  unchanged.

- **`CITATION.cff`** carried `date-released: 2026-04-21`, the v0.8.0 release date.

- **`conda-recipe/meta.yaml`** was missing from the v0.9.0 tree and has been
  restored at version 0.9.0.

### Changed

- Version strings bumped to `0.9.0` in `pyproject.toml`, `CITATION.cff`,
  `cascade/__init__.py`, the four sub-package docstrings, the startup banner and
  the `--version` output of all three entry points.
- `README.md` documents `--max_amplitude` in the flag table and the filtering
  examples, and the self-test block shows 11 tests.

### Migration notes

- No action required. `--max_amplitude` defaults to off and nothing else in the
  pipeline changed.
- Adopting `--max_amplitude` **will** change downstream results, because removing
  channels changes the spike pool feeding avalanche detection. Measured on the
  published 3Brain sample (`--max_amplitude 200`, 3 of 142 active channels removed,
  about 2%):

  | Metric | Without filter | With filter |
  |--------|---------------:|------------:|
  | τ (size exponent)   | 1.570 | 1.517 |
  | α (duration exponent) | 3.499 | 3.329 |
  | DCC                 | 2.561 | 2.400 |
  | β_opt (shape collapse) | 0.875 | 2.325 |
  | avalanche bin dt     | 6.30 ms | 6.35 ms |

  Removing 2% of channels moved τ by 0.05 — toward the theoretical 1.5 — and
  lowered DCC. β_opt moved substantially, so re-inspect `shape_collapse_<well>.png`
  whenever the filter is enabled rather than trusting the metric alone. These are
  real effects of excluding artefact channels, not a regression; report the
  threshold alongside the other filter settings in any methods section.
- On the published Axion sample the filter removed 0 of 70 channels and every
  output file was byte-identical to the unfiltered run, confirming the no-op path.
- Inspect a recording with `cascade-info` before choosing a threshold. 200 µV is a
  starting point, not a default suited to every culture or device.

### Known limitations

- Maxwell's `amplitude` field is a single-peak magnitude, not a true peak-to-peak
  measurement, so its values are not strictly comparable to the four loaders that
  compute `max − min` over a waveform. One threshold is not guaranteed to be
  equivalent across manufacturers; calibrate per device.
- The amplitude filter is all-or-nothing per channel: one artefact spike discards
  the whole channel for the entire recording. There is no per-spike rejection and
  no time-windowed variant.

---

## [0.8.0] — 2026-04-21

Peer-reviewed release. See `README.md` for the full feature set: five
manufacturer loaders (3Brain, Maxwell, Cortical Labs, MCS, Axion), spike
filtering, burst detection, avalanche detection, criticality metrics
(τ, α, DCC, branching ratio, shape collapse), functional connectivity, spatial
entropy, epoch comparison and LTP/LTD classification.

[0.9.5]: https://github.com/FA387/cascade-mea/releases/tag/v0.9.5
[0.9.4]: https://github.com/FA387/cascade-mea/tree/v0.9.4
[0.9.3]: https://github.com/FA387/cascade-mea/tree/v0.9.3
[0.9.2]: https://github.com/FA387/cascade-mea/tree/v0.9.2
[0.8.0]: https://github.com/FA387/cascade-mea/releases/tag/v0.8.0
