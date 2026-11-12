"""
core/criticality.py
-------------------
Criticality metrics: power-law exponents, DCC, branching ratio, shape collapse.

References
----------
Habibollahi et al. (2022)  – tplfit MLE method, DCC definition
Beggs & Plenz (2003)       – neuronal avalanche framework
Wilting & Priesemann (2018)– branching ratio estimation
Sethna (2001), Friedman et al. (2012) – shape collapse
"""

import numpy as np
import pandas as pd
from scipy.optimize import fmin
from .avalanche import avalanche_segments, counts_per_dt


# ── Truncated power-law fit (Habibollahi 2022 method) ─────────────────────────
def tplfit(data: np.ndarray, limit: int = 20, legacy_normalization: bool = False,
           xmin: int = None, clip_cdf: bool = True):
    """
    Fit a truncated power-law by finding the xmin that minimises the
    Kolmogorov–Smirnov (KS) statistic between the empirical CDF and the
    fitted truncated power-law CDF.

    Normalisation
    -------------
    The discrete truncated power law is normalised over EVERY integer between
    the cut-offs,

        p(y) = y^-a / sum_{k=xmin}^{xmax} k^-a

    which is the convention of Clauset et al. (2009) and of `plmle.py` in the
    NCC toolbox (Marshall et al. 2016).

    `legacy_normalization=True` sums only over the *observed unique* values,
    reproducing `tplfit.py` / `tplfit.m` from Habibollahi et al. (2022) and
    CASCADE <= 0.9.0. That convention makes the normalising constant depend on
    which values happen to appear in the sample, which biases the exponent
    (severely when avalanches are few) and leaves `A` inconsistent with the
    full-range KS grid below. It exists only to reproduce previously published
    numbers — do not use it for new analyses.

    Parameters
    ----------
    data  : 1-D array of positive integers (avalanche sizes or durations)
    limit : maximum xmin to try (search range: 1 … limit)
            Recommended: 8–20 for sizes, 3–6 for durations (Habibollahi 2022)
    legacy_normalization : reproduce the biased pre-0.9.1 behaviour (see above)
    xmin  : fit at this lower cut-off only instead of searching 1 … limit
            (the `bm` floor of Habibollahi's EX_burst=0 path). Default: search.
    clip_cdf : clip the fitted CDF to [0, 1] before the KS (CASCADE since 0.8.0).
            With legacy normalisation the fitted CDF can exceed 1; Habibollahi's
            tplfit.py does not clip, so KS and xmin can differ. The --pl_exclude
            path passes False to reproduce the reference. A no-op with the
            corrected normalisation, whose CDF never exceeds 1.

    Returns
    -------
    alpha   : float – fitted power-law exponent  P(x) ~ x^(-alpha)
    xmin    : int   – optimal lower cut-off (minimises KS)
    ks      : float – minimum KS statistic achieved
    loglike : float – log-likelihood at the optimum
    """
    data = np.asarray(data, dtype=np.float64)
    data = data[data > 0]
    if data.size < 5:
        return np.nan, 1, np.nan, np.nan

    xmax = int(data.max())

    if xmin is not None and int(xmin) < 1:
        raise ValueError(f"xmin must be >= 1, got {xmin}")
    candidates = [int(xmin)] if xmin is not None else list(range(1, int(limit) + 1))

    KS_list      = []
    alpha_list   = []
    loglike_list = []

    for x0 in candidates:
        idx = np.where((data >= x0) & (data <= xmax))[0]
        if idx.size < 3:
            KS_list.append(np.inf)
            alpha_list.append(np.nan)
            loglike_list.append(np.nan)
            continue

        subset = data[idx]
        n      = float(idx.size)
        # Support of the truncated distribution: every integer in [x0, xmax].
        # legacy_normalization restricts it to the observed values instead,
        # which is what Habibollahi 2022 tplfit.py does — see the docstring.
        s = (np.unique(subset) if legacy_normalization
             else np.arange(x0, xmax + 1, dtype=np.float64))

        # MLE: minimise negative log-likelihood of truncated power law
        # L(a) = a·Σlog(xi) - n·log(1/Σ s_k^(-a))
        def neg_ll(a):
            a = float(a[0]) if hasattr(a, '__len__') else float(a)
            denom = np.sum(np.power(s, -a))
            if denom <= 0:
                return 1e18
            return a * np.sum(np.log(subset)) - n * np.log(1.0 / denom)

        result = fmin(neg_ll, x0=[2.3], disp=False)
        a      = float(result[0])
        loglike_list.append(float(-neg_ll([a])))
        alpha_list.append(a)

        # KS statistic between empirical and fitted CDF
        z       = data[data >= x0]
        n_z     = float(z.size)
        bins    = np.arange(x0, xmax + 2)
        cdf_emp = np.cumsum(np.histogram(z, bins=bins)[0] / n_z)
        A       = 1.0 / np.sum(np.power(s, -a))
        fit_cdf = np.cumsum(A * np.power(np.arange(x0, xmax + 1, dtype=np.float64), -a))
        # clip fit_cdf to [0,1] to handle floating-point edge cases
        if clip_cdf:
            fit_cdf = np.clip(fit_cdf, 0.0, 1.0)
        KS_list.append(float(np.max(np.abs(cdf_emp - fit_cdf))))

    best     = int(np.argmin(KS_list))
    xmin_out = candidates[best]
    return (float(alpha_list[best]) if not np.isnan(alpha_list[best]) else np.nan,
            int(xmin_out),
            float(KS_list[best]),
            float(loglike_list[best]) if not np.isnan(loglike_list[best]) else np.nan)


# ── Upper-tail trimming (Habibollahi 2022 EXCLUDE, opt-in via --pl_exclude) ────
def exclude_tail(data: np.ndarray, limit: int, legacy_normalization: bool = False,
                 num: float = 1.0):
    """
    Port of `EXCLUDE.py` (Ma 2016, used by Habibollahi et al. 2022).

    Repeatedly fits tplfit(xmin search 1 … limit) and drops every copy of the
    current largest value until KS < min(num/sqrt(n), 0.1) or KS stops changing
    (|dKS| <= 0.0005). The stopping rule chases KS, so the trimmed range passes
    the fit by construction — see CHANGELOG before using it for new claims.

    Returns (xmax, xmin, alpha). xmax is the largest value left AFTER the last
    drop, exactly as the reference returns it (callers then keep <= xmax + 1).
    Raises ValueError when trimming empties the data, as the reference does.
    """
    burst = np.asarray(data, dtype=np.float64)
    KS, dKS, xmin, alpha, xmax = 1.0, 1.0, 1, np.nan, np.nan
    while True:
        n_above = np.sum(burst > xmin)
        thr = min(num / np.sqrt(n_above), 0.1) if n_above > 0 else 0.1
        if not (KS > thr and dKS > 0.0005):
            break
        if burst.size == 0:
            raise ValueError("EXCLUDE trimmed away all data")
        alpha, xmin, ks, _ = tplfit(burst, limit=limit, clip_cdf=False,
                                    legacy_normalization=legacy_normalization)
        if not np.isfinite(ks):
            raise ValueError("EXCLUDE: too few values left to fit")
        KS_old, KS = KS, ks
        dKS = abs(KS_old - KS)
        burst = burst[burst < burst.max()]
        if burst.size == 0:
            raise ValueError("EXCLUDE trimmed away all data")
        xmax = burst.max()
    return float(xmax), int(xmin), float(alpha)


def _fit_excluded(x: np.ndarray, floor: int, legacy_normalization: bool,
                  is_duration: bool):
    """
    One branch of Habibollahi's AV_analysis_ExponentErrorComments (EX_burst=1 /
    EX_time=1): pre-cut at max^0.8, EXCLUDE, then tplfit on [xmin, xmax(+1)].
    Returns (exponent, xmin, xmax, n_used, fell_back).
    """
    x = np.asarray(x, dtype=np.float64)
    top = np.max(x)
    try:
        xmax, xmin, _ = exclude_tail(x[x < np.power(top, 0.8)], limit=floor,
                                     legacy_normalization=legacy_normalization)
        if is_duration and xmax < 40:
            xmax = float(np.floor(np.power(top, 0.9)))
        sel = x[(x >= xmin) & (x <= xmax + 1)]
        fell_back = False
    except ValueError:
        # reference fallback: "could not minimize KS ... wants to exclude all data"
        xmax, xmin = float(np.power(top, 0.8)), int(floor)
        sel = x[(x >= xmin) & (x < xmax)]
        fell_back = True
    # final fit: the reference searches xmin 1 … xmin over data already >= xmin,
    # so every candidate sees the same subset and gives the same exponent
    expo, _, _, _ = tplfit(sel, limit=max(int(xmin), 1), clip_cdf=False,
                           legacy_normalization=legacy_normalization)
    return expo, int(xmin), float(xmax), int(sel.size), fell_back


# ── Compute all exponents + DCC ────────────────────────────────────────────────
def compute_exponents(sizes: np.ndarray,
                      durations: np.ndarray,
                      burst_min: int = 12,
                      t_min: int = 4,
                      legacy_normalization: bool = False,
                      xmin_s: int = None,
                      xmin_d: int = None,
                      exclude: bool = False) -> dict:
    """
    Fit power-law exponents for size (τ) and duration (α), then compute
    the predicted scaling exponent β_pred and DCC.

    Parameters
    ----------
    sizes     : 1-D array of avalanche sizes
    durations : 1-D array of avalanche durations (bins)
    burst_min : upper limit for xmin search in size distribution
                (Habibollahi 2022 default: 12; try 8–20)
    t_min     : upper limit for xmin search in duration distribution
                (Habibollahi 2022 default: 4; try 3–6)
    legacy_normalization : passed to tplfit(); reproduces the biased pre-0.9.1
                exponents. See tplfit() for why you should not use it.
    xmin_s, xmin_d : fix the lower cut-off instead of searching for it
                (--pl_xmin_S / --pl_xmin_D). Default None = automatic.
    exclude   : also trim the upper tail with Habibollahi's EXCLUDE (--pl_exclude).
                burst_min / t_min become EXCLUDE's bm / tm. Cannot be combined
                with xmin_s / xmin_d.

    Returns
    -------
    dict with keys: tau, xmin_s, alpha, xmin_d, beta_pred, beta_fit, dcc;
    with exclude also xmax_s, xmax_d, n_fit_s, n_fit_d, n_total,
    fallback_s, fallback_d
    """
    result = dict(tau=np.nan, xmin_s=np.nan, alpha=np.nan, xmin_d=np.nan,
                  beta_pred=np.nan, beta_fit=np.nan, dcc=np.nan)

    if sizes.size < 5 or durations.size < 5:
        return result

    if exclude:
        if xmin_s is not None or xmin_d is not None:
            raise ValueError("exclude cannot be combined with a fixed xmin")
        tau, xmin_s, xmax_s, n_s, fb_s = _fit_excluded(
            sizes, burst_min, legacy_normalization, is_duration=False)
        alpha, xmin_d, xmax_d, n_d, fb_d = _fit_excluded(
            durations, t_min, legacy_normalization, is_duration=True)
        result.update(tau=tau, xmin_s=xmin_s, alpha=alpha, xmin_d=xmin_d,
                      xmax_s=xmax_s, xmax_d=xmax_d, n_fit_s=n_s, n_fit_d=n_d,
                      n_total=int(sizes.size), fallback_s=fb_s, fallback_d=fb_d)
    else:
        # τ – size exponent
        tau, xmin_s, _, _ = tplfit(sizes, limit=burst_min,
                                   legacy_normalization=legacy_normalization,
                                   xmin=xmin_s)
        result["tau"]    = tau
        result["xmin_s"] = xmin_s

        # α – duration exponent
        alpha, xmin_d, _, _ = tplfit(durations, limit=t_min,
                                     legacy_normalization=legacy_normalization,
                                     xmin=xmin_d)
        result["alpha"]  = alpha
        result["xmin_d"] = xmin_d

    # β_pred = (α−1)/(τ−1)
    if np.isfinite(tau) and np.isfinite(alpha) and (tau != 1.0):
        result["beta_pred"] = (alpha - 1.0) / (tau - 1.0)

    # β_fit from <S>(T) ~ T^β  — fit range [xmin_d, xmin_d + 60] (Habibollahi 2022)
    beta_fit = _compute_beta_fit(sizes, durations, xmin_d if np.isfinite(xmin_d) else 1)
    result["beta_fit"] = beta_fit

    if np.isfinite(result["beta_pred"]) and np.isfinite(beta_fit):
        result["dcc"] = float(abs(result["beta_pred"] - beta_fit))

    return result


def _compute_beta_fit(sizes: np.ndarray,
                      durations: np.ndarray,
                      t_min: int) -> float:
    """
    Empirical scaling exponent β from  <S>(T) ~ T^β.

    Fit range is restricted to [t_min, t_min+60] following Habibollahi 2022.
    """
    TT  = np.arange(1, int(durations.max()) + 1)
    Sm  = np.array([
        np.mean(sizes[durations == t]) if np.any(durations == t) else np.nan
        for t in TT
    ])

    loc = np.where(Sm > 0)[0]
    if loc.size == 0:
        return np.nan

    TT = TT[loc]
    Sm = Sm[loc]

    # Habibollahi 2022: restrict fit to [tMin, tMin+60]
    fit_mask = (TT > int(t_min)) & (TT < int(t_min) + 60)
    if np.sum(fit_mask) < 2:
        # fall back to all available points if range is too narrow
        fit_mask = np.ones(TT.size, dtype=bool)
        if np.sum(fit_mask) < 2:
            return np.nan

    x = np.log(TT[fit_mask].astype(np.float64))
    y = np.log(Sm[fit_mask].astype(np.float64))
    beta_fit, _ = np.polyfit(x, y, 1)
    return float(beta_fit)


# ── Branching ratio ────────────────────────────────────────────────────────────
def compute_branching_ratio(counts: np.ndarray):
    """
    Estimate the branching ratio σ (Wilting & Priesemann 2018).

    Method: linear regression of E[N(t+1) | N(t)] vs N(t).
    Slope = branching ratio.  σ≈1 → critical; σ<1 → subcritical; σ>1 → supercritical.

    Parameters
    ----------
    counts : 1-D int array of population spike counts per time bin

    Returns
    -------
    BR        : float – branching ratio (slope of the regression)
    intercept : float – intercept of the regression line
    pts_df    : DataFrame – {Nt, E_Nt1_given_Nt, n_obs} for plotting
    """
    empty_df = pd.DataFrame(columns=["Nt", "E_Nt1_given_Nt", "n_obs"])

    counts = np.asarray(counts, dtype=np.int64)
    if counts.size < 2:
        return np.nan, np.nan, empty_df

    Nt  = counts[:-1]
    Nt1 = counts[1:]

    # keep only bins where Nt > 0 (avoid conditioning on silence)
    mask = Nt > 0
    Nt   = Nt[mask]
    Nt1  = Nt1[mask]

    if Nt.size < 2:
        return np.nan, np.nan, empty_df

    dfp = pd.DataFrame({"Nt": Nt, "Nt1": Nt1})
    grp = dfp.groupby("Nt")["Nt1"]
    pts = pd.DataFrame({
        "Nt"              : grp.mean().index.to_numpy(dtype=np.float64),
        "E_Nt1_given_Nt"  : grp.mean().to_numpy(dtype=np.float64),
        "n_obs"           : grp.size().to_numpy(dtype=np.int64),
    })

    if pts.shape[0] < 2:
        return np.nan, np.nan, pts

    x    = pts["Nt"].to_numpy(dtype=np.float64)
    y    = pts["E_Nt1_given_Nt"].to_numpy(dtype=np.float64)
    xbar = x.mean()
    ybar = y.mean()
    denom = float(np.sum((x - xbar) ** 2))
    if denom <= 0:
        return np.nan, np.nan, pts

    BR        = float(np.sum((x - xbar) * (y - ybar)) / denom)
    intercept = float(ybar - BR * xbar)
    return BR, intercept, pts


# ── Shape collapse ─────────────────────────────────────────────────────────────
def shape_collapse(counts: np.ndarray,
                   beta_pred: float,
                   Dmin: int = 4,
                   Dmax: int = 20,
                   x_step: float = 0.05,
                   min_size: int = 1):
    """
    Find the β that best collapses avalanche temporal profiles.

    The avalanche shape (mean profile vs scaled time t/D) should be
    universal under the scaling F(t/D) = D^(−γ) · <s(t|D)>, γ = β−1.
    We find β by minimising the variance of scaled profiles across durations.

    Based on Sethna (2001) / Friedman et al. (2012).

    Parameters
    ----------
    counts   : 1-D array of population counts per bin (from counts_per_dt)
    beta_pred: predicted scaling exponent (used to compute SC error)
    Dmin/Dmax: duration range (bins) to include in the collapse
    x_step   : resolution of the common scaled-time axis
    min_size : minimum total spike count for a profile to be included

    Returns
    -------
    beta_opt  : float  – β that minimises the collapse metric
    sce       : float  – |beta_opt − beta_pred|
    sc_minerr : float  – value of the collapse metric at beta_opt
    details   : dict   – full data for plotting
    """
    counts = np.asarray(counts, dtype=np.float64)
    segs   = avalanche_segments(counts)
    beta_grid = np.linspace(0.5, 2.5, 81)

    # group profiles by duration
    byD = {}
    for a, b in segs:
        D    = int(b - a)
        if D < Dmin or D > Dmax:
            continue
        prof = counts[a:b]
        if prof.sum() < float(min_size):
            continue
        byD.setdefault(D, []).append(prof)

    durations = sorted(byD.keys())
    if len(durations) < 2:
        return np.nan, np.nan, np.nan, dict(durations=durations)

    mean_profiles = {D: np.vstack([p for p in byD[D] if len(p) == D]).mean(axis=0)
                     for D in durations}

    x_common = np.arange(0.0, 1.0 + 1e-9, x_step)
    best_beta = np.nan
    best_err  = np.inf

    for beta in beta_grid:
        gamma  = beta - 1.0
        F_list = []
        for D in durations:
            s  = mean_profiles[D]
            F  = s / (float(D) ** gamma)
            x  = np.arange(D, dtype=np.float64) / float(D)
            Fi = np.interp(x_common, x, F)
            F_list.append(Fi)

        F_stack = np.vstack(F_list)
        varF    = np.var(F_stack, axis=0).mean()
        span    = float(np.max(F_stack) - np.min(F_stack))
        err     = float(varF / span ** 2) if span > 0 else np.inf

        if err < best_err:
            best_err  = err
            best_beta = float(beta)

    sce = float(abs(best_beta - float(beta_pred))) \
        if np.isfinite(best_beta) and np.isfinite(beta_pred) else np.nan

    details = dict(
        durations     = durations,
        x_common      = x_common,
        mean_profiles = mean_profiles,
        beta_grid     = beta_grid,
        beta_opt      = best_beta,
        err_min       = best_err,
    )
    return float(best_beta), sce, float(best_err), details
