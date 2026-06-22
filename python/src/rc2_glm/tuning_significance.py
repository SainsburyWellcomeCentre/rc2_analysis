"""Tuning-curve model selection + shuffle-null significance (pure numpy/scipy).

A Python re-implementation of the MATLAB ``lib/classes/tuning`` family
(`ShuffleTuning` / `ModelSelectionTuning` / `AsymmetricGaussianFit`), adapted
for the goggles trial-structure figure. It tests whether a cluster's observed
FR-vs-stimulus-value tuning curve (TF / SF / OR) is real, and picks the model
family that best describes its shape so the figure can overlay it.

No GLM coupling — this is a standalone helper (like ``precomputed_bins`` /
``cross_validation``). The figure and any cohort table both call
``per_trial_bin_matrix`` so there is exactly **one** binning path.

Decisions (2026-06-19, Laura) — and where they diverge from the MATLAB:

* **Test statistic = R² on the MEAN curve** (the per-bin trial-mean, ``n_bins``
  points), not the MATLAB flattened per-trial-per-bin R². Per-trial-per-bin R²
  is ~0.03 on Poisson V1 (trial noise dominates) and the test is then nearly
  powerless; the mean curve asks the question we actually care about — *is the
  mean tuning SHAPE beyond chance*. (`AsymmetricGaussianFit` already fits the
  mean curve, so this matches one of the two MATLAB conventions.)
* **Model selection = BIC** (lowest), not the MATLAB R²-on-mean (which always
  favours the most flexible 4-param curve). BIC penalises the Gaussian/sigmoid
  families so we don't over-select them.
* **Null = the MATLAB bootstrap** — resample the whole ``(n_trials × n_bins)``
  matrix with replacement (``rng.integers(0, size, size)``), recompute the mean
  curve, refit the selected family, ``p = Σ(rsq_shuff ≥ rsq_obs) / n_reps``.
  This destroys bin↔FR association and re-weights values, exactly as
  ``rng(1); I = randi(numel(tuning), size(tuning))`` did.
* **OR is circular** (0–180°): fit only a 180°-period von Mises
  ``A·exp(κ(cos 2(θ−μ) − 1)) + b``; SF / TF keep linear/quadratic/cubic/
  Gaussian/asymmetric-Gaussian/sigmoid.

Determinism: the bootstrap draws from a seeded ``numpy.random.Generator``
(default seed 1, matching the MATLAB ``rng(1)``). The per-curve least-squares
fits run on at most ``n_bins`` (=20) points, so they are not the threaded-IRLS
path that wanders across BLAS thread counts
(cf. ``reference_rc2_glm_cvbps_blas_nondeterminism``) — no thread pinning needed.

Degenerate cells return ``p = NaN`` (never a fabricated curve): too few trials,
too few finite bins, no value variation, or a flat mean curve (TSS = 0).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit

# Families fit on a linear value axis (TF, SF). OR uses the circular family.
LINEAR_FAMILIES = (
    "linear",
    "quadratic",
    "cubic",
    "gaussian",
    "asym_gaussian",
    "sigmoid",
)
CIRCULAR_FAMILIES = ("vonmises_180",)

# Guards (don't fabricate a curve from noise) — see module docstring.
MIN_TRIALS = 3
MIN_FINITE_BINS = 4


# --------------------------------------------------------------------------- #
# Shared binning — the ONE per-trial-per-bin matrix builder.
# --------------------------------------------------------------------------- #
def per_trial_bin_matrix(df, value_col, bw, *, condition, n_bins=20, edges=None):
    """``(n_trials × n_bins)`` per-trial mean-FR matrix + bin centres.

    Per trial, the mean FR (``spike_count / bw``) of the bins falling in each
    value bin; empty bins for a trial stay ``NaN``. Bins are 20 equal-count
    (5%-quantile) ``_pooled_quantile_edges``. **Pass ``edges`` to bin into a
    SHARED set computed elsewhere** (e.g. pooled across V+VT — the project
    convention, so bin k is the same interval in every condition); if ``edges``
    is None they are computed from this condition's values alone (only correct
    when the caller has already restricted to a single comparable pool).

    Returns ``(matrix, centres)`` with ``matrix`` shape ``(n_trials, n_bins)``,
    or ``(None, None)`` if there is nothing to bin (no rows, missing column,
    degenerate value range).
    """
    # Local import to avoid a heavy plots.py import at module load.
    from rc2_glm.plots import _pooled_quantile_edges

    if value_col not in df.columns:
        return None, None
    sub = df[df["condition"] == condition]
    if sub.empty:
        return None, None

    vals = sub[value_col].to_numpy(np.float64)
    finite = np.isfinite(vals)
    if finite.sum() < 2:
        return None, None
    if edges is None:
        if finite.sum() < n_bins:
            return None, None
        edges, centres = _pooled_quantile_edges(vals[finite], n_bins=n_bins)
        if edges is None:
            return None, None
    else:
        edges = np.asarray(edges, np.float64)
        centres = 0.5 * (edges[:-1] + edges[1:])
    n_b = edges.size - 1

    fr = sub["spike_count"].to_numpy(np.float64) / bw
    tids = sub["trial_id"].to_numpy(np.int64)

    rows = []
    for tid in np.unique(tids):
        tm = (tids == tid) & finite
        if tm.sum() < 2:
            continue
        bidx = np.clip(np.digitize(vals[tm], edges) - 1, 0, n_b - 1)
        ftm = fr[tm]
        row = np.full(n_b, np.nan, dtype=np.float64)
        for b in range(n_b):
            inb = bidx == b
            if inb.any():
                row[b] = float(np.nanmean(ftm[inb]))
        rows.append(row)
    if not rows:
        return None, None
    return np.asarray(rows, dtype=np.float64), centres


# --------------------------------------------------------------------------- #
# Model families (fit on the MEAN curve).
# --------------------------------------------------------------------------- #
def _vonmises_180(x_deg, amp, kappa, mu_deg, base):
    """180°-period (axial) von Mises: peak ``amp + base`` at θ = mu."""
    th = np.deg2rad(x_deg)
    mu = np.deg2rad(mu_deg)
    return amp * np.exp(kappa * (np.cos(2.0 * (th - mu)) - 1.0)) + base


def _gaussian(x, amp, mu, sigma, base):
    return amp * np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2)) + base


def _asym_gaussian(x, r_max, x_max, sigma_minus, sigma_plus):
    x = np.asarray(x, float)
    sig = np.where(x < x_max, sigma_minus, sigma_plus)
    return r_max * np.exp(-((x - x_max) ** 2) / sig)


def _sigmoid(x, amp, k, x0, base):
    # Clip the logistic argument to the float64-safe range so curve_fit probing
    # large k doesn't raise overflow warnings (the saturated value is unchanged).
    z = np.clip(-k * (np.asarray(x, float) - x0), -700.0, 700.0)
    return amp / (1.0 + np.exp(z)) + base


def _metrics(y, y_pred, k):
    """(rsq, bic) on the mean curve; matches the MATLAB Gaussian-LL BIC.

    ``rsq = 1 − RSS/TSS``; ``BIC = −2·LL + k·ln(n)`` with the Gaussian
    log-likelihood ``LL = −0.5·n·(ln(2π·RSS/n) + 1)``. Returns ``(-inf, inf)``
    on a numerically dead fit so it loses both selection and the test.
    """
    n = y.size
    rss = float(np.sum((y - y_pred) ** 2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    if tss <= 0:
        return np.nan, np.inf
    rsq = 1.0 - rss / tss
    if rss <= 0:
        return rsq, -np.inf
    sigma2 = rss / n
    ll = -0.5 * n * (np.log(2.0 * np.pi * sigma2) + 1.0)
    bic = -2.0 * ll + k * np.log(n)
    return rsq, bic


def _fit_one(x, y, family):
    """Fit one family to finite mean-curve points ``(x, y)``.

    Returns ``dict(name, params, n_params, rsq, bic, ok)``; ``ok=False`` (with
    ``rsq=-inf, bic=inf``) if the curve_fit raised or could not be evaluated.
    """
    xr = float(x.max() - x.min())
    try:
        if family in ("linear", "quadratic", "cubic"):
            deg = {"linear": 1, "quadratic": 2, "cubic": 3}[family]
            beta = np.polyfit(x, y, deg)
            y_pred = np.polyval(beta, x)
            k = deg + 1
            params = beta
        elif family == "gaussian":
            amax = int(np.argmax(y))
            p0 = [y.max() - y.min(), x[amax], max(xr / 4.0, 1e-6), y.min()]
            lb = [0.0, x.min(), 1e-3, -np.inf]
            ub = [np.inf, x.max(), np.inf, np.inf]
            params, _ = curve_fit(_gaussian, x, y, p0=p0, bounds=(lb, ub), maxfev=10000)
            y_pred = _gaussian(x, *params)
            k = 4
        elif family == "asym_gaussian":
            amax = int(np.argmax(y))
            s0 = max(abs(xr) / 4.0, 1e-6)
            p0 = [y.max(), x[amax], s0, s0]
            lb = [0.0, x.min(), 1e-3, 1e-3]
            ub = [y.max() * 2.0 + 1e-9, x.max(), abs(xr) * 2.0 + 1e-9, abs(xr) * 2.0 + 1e-9]
            params, _ = curve_fit(_asym_gaussian, x, y, p0=p0, bounds=(lb, ub), maxfev=10000)
            y_pred = _asym_gaussian(x, *params)
            k = 4
        elif family == "sigmoid":
            sx = float(np.std(x)) or 1.0
            p0 = [y.max() - y.min(), 1.0 / sx, float(np.mean(x)), y.min()]
            lb = [0.0, 0.0, x.min(), -np.inf]
            ub = [np.inf, np.inf, x.max(), np.inf]
            params, _ = curve_fit(_sigmoid, x, y, p0=p0, bounds=(lb, ub), maxfev=10000)
            y_pred = _sigmoid(x, *params)
            k = 4
        elif family == "vonmises_180":
            amax = int(np.argmax(y))
            p0 = [max(y.max() - y.min(), 1e-6), 2.0, float(x[amax]) % 180.0, y.min()]
            lb = [0.0, 1e-3, 0.0, -np.inf]
            ub = [np.inf, 50.0, 180.0, np.inf]
            params, _ = curve_fit(_vonmises_180, x, y, p0=p0, bounds=(lb, ub), maxfev=10000)
            y_pred = _vonmises_180(x, *params)
            k = 4
        else:
            raise ValueError(f"unknown family {family!r}")
        if not np.all(np.isfinite(y_pred)):
            raise ValueError("non-finite prediction")
    except Exception:
        return dict(name=family, params=None, n_params=0, rsq=-np.inf, bic=np.inf, ok=False)

    rsq, bic = _metrics(y, y_pred, k)
    return dict(name=family, params=np.asarray(params, float), n_params=k,
                rsq=rsq, bic=bic, ok=True)


def evaluate(fit, x_grid):
    """Evaluate a fit dict on ``x_grid`` (for the figure overlay)."""
    name, p = fit["name"], fit["params"]
    if p is None:
        return np.full_like(np.asarray(x_grid, float), np.nan)
    if name in ("linear", "quadratic", "cubic"):
        return np.polyval(p, x_grid)
    if name == "gaussian":
        return _gaussian(x_grid, *p)
    if name == "asym_gaussian":
        return _asym_gaussian(x_grid, *p)
    if name == "sigmoid":
        return _sigmoid(x_grid, *p)
    if name == "vonmises_180":
        return _vonmises_180(x_grid, *p)
    raise ValueError(f"unknown family {name!r}")


def _mean_curve(matrix):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(matrix, axis=0)


def _fit_data_xy(matrix, centres, aggregate):
    """The (x, y) the model is fit to / scored on.

    - ``"mean"``: the per-bin mean curve (``n_bins`` points) — high R², asks "is
      the mean tuning SHAPE real".
    - ``"flat"``: **every** per-trial-per-bin point, ``x`` = bin centres tiled
      across trials (the MATLAB ``ShuffleTuning`` convention — R² computed across
      all trials). Low R² (per-trial Poisson scatter dominates), but it's the
      lab's prior convention.

    NaNs dropped.
    """
    centres = np.asarray(centres, float)
    if aggregate == "flat":
        x = np.tile(centres, matrix.shape[0])
        y = np.asarray(matrix, float).ravel()
    else:
        x = centres
        y = _mean_curve(matrix)
    finite = np.isfinite(x) & np.isfinite(y)
    return x[finite], y[finite]


def rsq_against_mean(matrix, centres, family, params):
    """R² of a fixed fitted curve against the per-bin MEAN tuning curve.

    How well the selected model (its ``params``, whatever data they were fit on)
    describes the mean curve — i.e. the fit-vs-mean agreement, computed on the
    ``n_bins`` mean points. NaN if too few finite bins or a flat mean (TSS=0)."""
    if params is None:
        return np.nan
    x, y = _fit_data_xy(matrix, centres, "mean")
    if x.size < MIN_FINITE_BINS:
        return np.nan
    y_pred = evaluate({"name": family, "params": params}, x)
    if not np.all(np.isfinite(y_pred)):
        return np.nan
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return np.nan if ss_tot <= 0 else 1.0 - ss_res / ss_tot


def fit_tuning(centres, mean_curve, *, kind, linear_families=LINEAR_FAMILIES):
    """Fit the candidate families for ``kind`` ('linear' or 'circular').

    ``linear_families`` restricts the set tried on the linear axis (TF/SF) — e.g.
    ``("asym_gaussian",)`` to fit only the asymmetric Gaussian. The circular axis
    (OR) always uses :data:`CIRCULAR_FAMILIES`. Returns a list of fit dicts; use
    :func:`select_best`.
    """
    finite = np.isfinite(mean_curve)
    x, y = np.asarray(centres, float)[finite], np.asarray(mean_curve, float)[finite]
    families = CIRCULAR_FAMILIES if kind == "circular" else tuple(linear_families)
    if x.size < MIN_FINITE_BINS:
        return []
    return [_fit_one(x, y, fam) for fam in families]


def select_best(fits):
    """Lowest-BIC fit (ties broken by higher rsq). ``None`` if none usable."""
    usable = [f for f in fits if f["ok"] and np.isfinite(f["bic"])]
    if not usable:
        return None
    return min(usable, key=lambda f: (f["bic"], -f["rsq"]))


def select_best_rsq_mean(fits, matrix, centres):
    """Highest R²-on-the-mean-curve fit (the MATLAB ``ModelSelectionTuning`` rule).

    Each candidate is fit on the (flat) data, then scored by how well its curve
    matches the per-bin MEAN — exactly ``select_best_model``'s criterion. ``None``
    if none usable. Attaches the mean-curve R² to each fit as ``rsq_mean_sel``."""
    usable = []
    for f in fits:
        if not f["ok"]:
            continue
        f["rsq_mean_sel"] = rsq_against_mean(matrix, centres, f["name"], f["params"])
        if np.isfinite(f["rsq_mean_sel"]):
            usable.append(f)
    if not usable:
        return None
    return max(usable, key=lambda f: f["rsq_mean_sel"])


# --------------------------------------------------------------------------- #
# Shuffle null + orchestration.
# --------------------------------------------------------------------------- #
@dataclass
class TuningSignificance:
    value: str                  # "tf" | "sf" | "or"
    condition: str              # "V" | "VT"
    kind: str                   # "linear" | "circular"
    n_trials: int = 0
    n_bins: int = 0
    best_model: str | None = None
    params: np.ndarray | None = None
    rsq: float = np.nan          # R² of the selected model (see `aggregate`)
    rsq_mean: float = np.nan     # R² of that fitted curve vs the per-bin MEAN curve
    bic: float = np.nan
    p: float = np.nan
    aggregate: str = "mean"      # "mean" (mean-curve R²) | "flat" (across-trials R²)
    select_criterion: str = "bic"  # how the family was chosen: "bic" | "rsq_mean"
    data_source: str = ""        # where the tuning matrix came from (matlab_cache | recomputed)
    n_reps: int = 0
    seed: int = 1
    null_scheme: str = "bootstrap_resample_with_replacement"
    centres: np.ndarray | None = None
    mean_curve: np.ndarray | None = None
    null_rsq: np.ndarray | None = field(default=None, repr=False)


def shuffle_test(matrix, centres, family, *, kind, aggregate="mean",
                 n_reps=1000, seed=1):
    """Bootstrap-null p for ``family`` on ``matrix``.

    Resample the whole matrix with replacement ``n_reps`` times, rebuild the
    fit data for ``aggregate`` ('mean' curve | 'flat' across-trials points),
    refit ``family``, ``p = Σ(rsq_shuff ≥ rsq_obs)/n_reps``. Returns
    ``(p, rsq_obs, null_rsq)``.
    """
    x_obs, y_obs = _fit_data_xy(matrix, centres, aggregate)
    if x_obs.size < MIN_FINITE_BINS:
        return np.nan, np.nan, np.array([])
    rsq_obs = _fit_one(x_obs, y_obs, family)["rsq"]
    if not np.isfinite(rsq_obs):
        return np.nan, rsq_obs, np.array([])

    rng = np.random.default_rng(seed)
    flat = matrix.ravel()
    size = flat.size
    null = np.full(n_reps, np.nan)
    for i in range(n_reps):
        samp = flat[rng.integers(0, size, size)].reshape(matrix.shape)
        xs, ys = _fit_data_xy(samp, centres, aggregate)
        if xs.size < MIN_FINITE_BINS:
            continue
        null[i] = _fit_one(xs, ys, family)["rsq"]
    valid = null[np.isfinite(null)]
    if valid.size == 0:
        return np.nan, rsq_obs, valid
    p = float(np.sum(valid >= rsq_obs) / valid.size)
    return p, rsq_obs, valid


def tuning_significance(matrix, centres, *, value, condition, kind,
                        aggregate="mean", select_criterion="bic",
                        n_reps=1000, seed=1, linear_families=LINEAR_FAMILIES,
                        min_trials=MIN_TRIALS, min_finite_bins=MIN_FINITE_BINS):
    """End-to-end: build fit data → fit families → select → shuffle-null p.

    ``matrix`` is ``(n_trials × n_bins)`` from :func:`per_trial_bin_matrix`.
    ``aggregate`` chooses the R²/fit target: ``"mean"`` (per-bin mean curve) or
    ``"flat"`` (every per-trial-per-bin point — the MATLAB ``ShuffleTuning`` R²
    across all trials). ``select_criterion`` chooses the family: ``"bic"`` (lowest
    BIC) or ``"rsq_mean"`` (highest R²-on-the-mean-curve, the MATLAB
    ``ModelSelectionTuning`` rule). The selected family is then tested the same
    way regardless (shuffle null on the ``aggregate`` R²). ``linear_families``
    restricts the TF/SF family set; OR always uses the circular family.
    Degenerate inputs yield ``p=NaN`` with ``best_model=None``.
    """
    res = TuningSignificance(value=value, condition=condition, kind=kind,
                             aggregate=aggregate, select_criterion=select_criterion,
                             n_reps=n_reps, seed=seed)
    if matrix is None or matrix.size == 0:
        return res
    res.n_trials, res.n_bins = matrix.shape
    res.centres = np.asarray(centres, float)
    mc = _mean_curve(matrix)
    res.mean_curve = mc

    finite = np.isfinite(mc)
    if res.n_trials < min_trials or finite.sum() < min_finite_bins:
        return res  # degenerate → p stays NaN

    x_fit, y_fit = _fit_data_xy(matrix, centres, aggregate)
    families = CIRCULAR_FAMILIES if kind == "circular" else tuple(linear_families)
    fits = [_fit_one(x_fit, y_fit, fam) for fam in families]
    best = (select_best_rsq_mean(fits, matrix, centres)
            if select_criterion == "rsq_mean" else select_best(fits))
    if best is None:
        return res
    res.best_model = best["name"]
    res.params = best["params"]
    res.rsq = best["rsq"]
    res.rsq_mean = rsq_against_mean(matrix, centres, best["name"], best["params"])
    res.bic = best["bic"]

    p, _rsq, null = shuffle_test(matrix, centres, best["name"], kind=kind,
                                 aggregate=aggregate, n_reps=n_reps, seed=seed)
    res.p = p
    res.null_rsq = null
    return res
