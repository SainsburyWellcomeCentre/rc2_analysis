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

    # --- Spike history (prompt 03, 2026-04-28; default ON since 2026-04-29) ---
    # 10 log-spaced raised-cosine bases over a 200 ms post-spike window,
    # added as a Phase-1 forward-selection candidate. Each cluster's
    # history features are convolved trial-aware (zero-padded at trial
    # starts; lag 0 excluded for causality).
    #
    # Default flipped from False to True on 2026-04-29 after the prompt-03
    # ablation showed history dominates Phase-1 selection in 32/33 clusters
    # of the smoke probe and adds median Δ +0.071 bps on the 88-cluster
    # filtered set (76/88 positive). Turn off via the `--no-history` CLI
    # flag for one-off comparisons (or for legacy reruns of the
    # MATLAB-parity model).
    include_history: bool = True
    n_history_bases: int = 10
    history_window_s: float = 0.2
    # When False (default), History interacts with nothing in Phase 2 —
    # interaction interpretations are rarely useful for spike history.
    allow_history_interactions: bool = False
    # Onset kernel inclusion (default OFF since 2026-04-29).
    #
    # Default flipped from True to False on 2026-04-29 after the prompt-03
    # ablation showed the onset kernel adds ~0 CV-bps once history is
    # included (median Δ(C-A) = +0.066 ≈ Δ(B-A) = +0.071, i.e. removing
    # onset costs nothing when history is present). The kernel basis +
    # design-matrix wiring stay in place for occasional ablation reruns
    # via the `--with-onset-kernel` opt-in flag.
    #
    # Consequence: MATLAB parity claim retired on 2026-04-29 — the Null
    # model is no longer "intercept + onset" so cv-bps comparisons against
    # the MATLAB reference no longer hold gate-by-gate. New trust signal:
    # python/tests/test_pipeline_regression.py.
    include_onset_kernel: bool = False

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
    # is "trial-averaged". One of:
    #   "none"             — no band; reproduces pre-prompt-12 line-only output.
    #   "covariate-spread" — IQR across trials of the predicted rate at
    #                        fixed sweep-x given each trial's actual non-
    #                        target covariates. Reflects model-structural
    #                        sensitivity to TF/SF/OR heterogeneity.
    #                        (Was called "iqr" pre-2026-04-28.)
    #   "simulated" (default) — parametric bootstrap: predict λ at the
    #                        training granularity (100 ms bins) using each
    #                        trial's actual covariates; draw y_sim ~
    #                        Poisson(λ · Δt); collapse simulated rates to
    #                        the cache's display bins; band = IQR across
    #                        trials per display bin (mean across MC iterations).
    #                        Directly comparable to the Observed row's
    #                        whiskers — same reconstruction procedure,
    #                        simulated vs real spike counts.
    # Aliases:
    #   "iqr"  → maps to "covariate-spread" (deprecated alias for back-compat).
    # Sparse-bin guard: bins with < MIN_TRIALS_FOR_BAND contributions are
    # skipped (no fake bands at single-trial bins).
    tuning_curve_uncertainty: str = "simulated"

    # Number of Monte Carlo iterations for the parametric-bootstrap
    # ("simulated") band. 100 is the standard floor for IQR stability;
    # drop to 50 if compute is tight.
    n_bootstrap_iterations: int = 100


# Minimum number of trials at a display-bin for the uncertainty band to
# be drawn there. Below this, the band is skipped and a one-line warning
# is logged.
MIN_TRIALS_FOR_BAND: int = 3


INTERACTION_PARENTS: dict[str, tuple[str, str]] = {
    "Speed_x_TF": ("Speed", "TF"),
    "Speed_x_SF": ("Speed", "SF"),
    "Speed_x_OR": ("Speed", "OR"),
    "TF_x_SF": ("TF", "SF"),
    "TF_x_OR": ("TF", "OR"),
    "SF_x_OR": ("SF", "OR"),
}
