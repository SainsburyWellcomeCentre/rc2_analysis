"""Default parameters matching MATLAB v10 of glm_single_cluster_analysis.m."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class GLMConfig:
    # --- Time binning ---
    time_bin_width: float = 0.1                # seconds (100 ms)
    motion_fraction_threshold: float = 0.5      # min fraction of motion samples per bin

    # --- Motion mask (matches MATLAB Trial.treadmill_motion_mask defaults) ---
    velocity_threshold: float = 1.0             # cm/s
    acceleration_threshold: float = 0.5         # m/s^2
    min_stationary_duration: float = 0.2        # seconds

    # --- Velocity filter (matches lib/fcn/general/filter_trace.m) ---
    apply_velocity_filter: bool = True
    filter_cutoff_hz: float = 50.0
    filter_order: int = 3

    # --- Basis functions ---
    n_speed_bases: int = 5
    speed_range: tuple[float, float] = (0.0, 50.0)        # cm/s
    n_tf_bases: int = 5
    tf_range: tuple[float, float] = (0.0, 7.3)            # Hz
    n_onset_bases: int = 6
    onset_range: tuple[float, float] = (0.0, 2.0)         # seconds

    # --- GLM fitting ---
    # Ridge on all non-intercept columns. Non-zero default picks a
    # specific rotation in the flat direction of the Poisson likelihood
    # created by correlated raised-cosine bases (Speed/TF bases are
    # correlated by construction — Park et al. 2014), which is what
    # makes tuning-curve parity with MATLAB glmnet's lambda_1se
    # achievable. 1e-3 was tuned empirically on CAA-1123243_rec1 as
    # the largest value that still improves tuning-curve Pearson r vs
    # MATLAB; 1e-2 over-shrinks (hurts parity), 0 leaves β rotation
    # unconstrained (also hurts parity). FullInteraction keeps its own
    # (larger) lambda because p approaches n there.
    lambda_ridge: float = 1e-3
    full_interaction_lambda: float = 1.0
    lambda_ridge_min: float = 1e-6
    irls_max_iter: int = 100
    irls_tol: float = 1e-8
    eta_clip: float = 20.0
    mu_floor: float = 1e-10

    # --- Cross-validation ---
    n_folds: int = 5
    cv_seed: int = 0
    # "condition-stratified" (default, MATLAB sp_fold parity) stratifies
    # trial-level k-fold over unique (trial_id, condition) pairs.
    # "speed-profile" uses 2 folds keyed on ``TrialData.profile_id`` —
    # train on one reproduced velocity trajectory, test on the other.
    # Mirrors MATLAB glm_single_cluster_analysis.m:2291-2293 and lets
    # us quantify how well the GLM generalises across speed profiles.
    cv_strategy: str = "condition-stratified"
    # When True, run a post-hoc speed-profile CV diagnostic on top of the
    # normal forward selection (which keeps condition-stratified folds):
    # for each cluster re-compute CV-bps on Null, Selected, and
    # Selected-without-Speed under profile folds, emit comparison CSV
    # columns + the MATLAB-parity PDF figure. Does NOT change the
    # primary fit — only adds the extra diagnostic pass. Mirrors MATLAB
    # glm_single_cluster_analysis.m:2246-2367.
    profile_cv_diagnostic: bool = False

    # --- Forward selection ---
    main_effects: tuple[str, ...] = ("Speed", "TF", "SF", "OR")
    interactions: tuple[str, ...] = (
        "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
        "TF_x_SF", "TF_x_OR", "SF_x_OR",
    )
    delta_bps_threshold: float = 0.005

    # --- Reference levels (matches MATLAB sf_levels / or_levels) ---
    sf_levels: tuple[float, ...] = (0.003, 0.006, 0.012)
    or_levels: tuple[float, ...] = field(
        default_factory=lambda: (-0.7853981633974483, 0.0, 0.7853981633974483, 1.5707963267948966)
    )  # -π/4, 0, π/4, π/2 (sorted)

    # --- Prefilter ---
    apply_prefilter: bool = True
    prefilter_seed: int = 0

    # --- Compute backend ---
    device: str = "auto"  # "auto" | "cpu" | "gpu"

    # --- Tuning-curve rendering ---
    # "trial-averaged" (default): predicted tuning at each grid point is
    # averaged over the cluster's observed motion-bin distribution per
    # condition (marginalises out time_since_onset and the other covariates
    # not fixed by the MATLAB convention). This gives the trial-averaged
    # expected rate the neuron would actually show, not a hypothetical
    # steady-state evaluation. "steady-state": evaluate the onset kernel
    # at t=1.5s and the other continuous variable at its per-cluster
    # mean (the pre-2026-04-23 behaviour, kept for back-compat).
    tuning_curve_mode: str = "trial-averaged"

    # Per-trial uncertainty band drawn around each model-row line on the
    # cluster_<id>_tuning.pdf panels. Computed only when tuning_curve_mode
    # is "trial-averaged" (steady-state evaluates one point per condition;
    # there is no per-trial spread to summarise). One of:
    #   "none"           — no band; reproduces the pre-prompt-12 line-only output.
    #   "iqr" (default)  — fill between q25 and q75 of per-trial predicted rates.
    #   "wide-quantile"  — fill between q05 and q95.
    #   "std"            — fill mean ± std of per-trial predicted rates.
    # Sparse-bin guard: grid points with < MIN_TRIALS_FOR_BAND per-trial
    # contributions are skipped (no fake bands at single-trial bins).
    tuning_curve_uncertainty: str = "iqr"


# Minimum number of trials at a tuning-curve grid point for the
# uncertainty band to be drawn there. Below this, the band is skipped at
# that grid point and a one-line warning is logged.
MIN_TRIALS_FOR_BAND: int = 3


INTERACTION_PARENTS: dict[str, tuple[str, str]] = {
    "Speed_x_TF": ("Speed", "TF"),
    "Speed_x_SF": ("Speed", "SF"),
    "Speed_x_OR": ("Speed", "OR"),
    "TF_x_SF": ("TF", "SF"),
    "TF_x_OR": ("TF", "OR"),
    "SF_x_OR": ("SF", "OR"),
}
