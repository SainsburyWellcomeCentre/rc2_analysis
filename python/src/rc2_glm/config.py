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

    # Spacing for the Speed/TF value-axis raised-cosine bases. "log"
    # (default, Weber-law via raised_cosine_basis — dense at low values,
    # MATLAB parity) or "linear" (even tiling via raised_cosine_basis_linear
    # — better mid-range resolution for a central tuning peak, at the cost
    # of low-value resolution). Applies to the fit bases, the kernel
    # reconstruction (_kernel_for_var) and the tuning-curve grid. NOTE:
    # secondary plot panels (model overview, Speed×TF interaction heatmap)
    # still assume log spacing; this switch is for value-axis fit/kernel/
    # tuning experiments (scripts/smoke_basis_count_cluster.py, 2026-06-08).
    speed_tf_basis_spacing: str = "log"

    # --- Spike history (prompt 03, 2026-04-28; default ON since 2026-04-29; default OFF since 2026-04-30) ---
    # 10 log-spaced raised-cosine bases over a 200 ms post-spike window,
    # added as a Phase-1 forward-selection candidate. Each cluster's
    # history features are convolved trial-aware (zero-padded at trial
    # starts; lag 0 excluded for causality).
    #
    # Default flipped from False to True on 2026-04-29 after the prompt-03
    # ablation showed history dominates Phase-1 selection in 32/33 clusters
    # of the smoke probe and adds median Δ +0.071 bps on the 88-cluster
    # filtered set (76/88 positive).
    #
    # Default flipped back from True to False on 2026-04-30 (prompt 06).
    # The 2026-04-29 history-overweighted analysis showed 66/79 clusters
    # peak in trial-level Pearson r at α < 1 — a sign that the History
    # term was absorbing structure beyond its autoregressive role
    # (NLL-Pearson disagreement consistent with misspecification
    # absorption). The component reads as a catch-all rather than an
    # interpretable scientific quantity, defeating the explanatory
    # purpose of the GLM. Standing by until a more interpretable
    # parameterisation lands (per-cluster instead of session-shared,
    # decoupled history-bin from GLM-bin, or hierarchical priors).
    # CLI flag --include-history still works for ad-hoc experiments.
    include_history: bool = False
    # Reduced from 10 to 5 on 2026-04-29 after the basis-count sweep.
    # At 100 ms only 2 lag bins are distinct (n_lag_bins = window/bin_width
    # = 0.2/0.1 = 2), so 10 bases over 2 lag bins is purely a basis-rotation
    # ambiguity — predictions are identical. 5 bases gives the same fits
    # with half the History coefficient columns. At 20 ms with 10 lag
    # bins the 5-basis fit also matches the 10-basis fit (kernels in
    # lag space overlap). Verified: median Δ cv_bps = 0.000000 across
    # 88 clusters at 100 ms; 5/10 selected_vars identical 88/88.
    n_history_bases: int = 5
    history_window_s: float = 0.2
    # History basis kind: "raised_cosine" (default) or "identity" (one dummy
    # per lag bin). For short windows with few lag bins (e.g. 20 ms bins,
    # 80–100 ms window → 4–5 lags) raised cosines are just a rotation of the
    # dummies — use "identity" for interpretability (coef = effect of a spike
    # N bins ago). n_history_bases is ignored when "identity".
    history_basis_kind: str = "raised_cosine"
    # When False (default), History interacts with nothing in Phase 2 —
    # interaction interpretations are rarely useful for spike history.
    allow_history_interactions: bool = False
    # Onset kernel inclusion (default OFF since 2026-04-29; default ON since 2026-04-30).
    #
    # Default flipped from True to False on 2026-04-29 after the prompt-03
    # ablation showed the onset kernel adds ~0 CV-bps once history is
    # included (median Δ(C-A) = +0.066 ≈ Δ(B-A) = +0.071, i.e. removing
    # onset costs nothing when history is present). The kernel basis +
    # design-matrix wiring stay in place for occasional ablation reruns
    # via the `--with-onset-kernel` opt-in flag.
    #
    # Default flipped back from False to True on 2026-04-30 (prompt 06).
    # When history is on stand-by, the onset kernel reclaims its role of
    # capturing the trial-onset response shape — without it, the Null
    # model loses the "intercept + onset" baseline that anchors the
    # cv_bps scale. Restoring onset returns the GLM to the pre-2026-04-29
    # configuration, with the addition of ME_face as the new Phase-1
    # candidate. Consequence: cv_bps numbers on the new run ARE comparable
    # in absolute scale to legacy_with_onset/, but the new run also has
    # ME_face columns so model-comparison is one-step-removed from
    # legacy parity.
    include_onset_kernel: bool = True

    # --- Face motion energy (prompt 06, 2026-04-30) ---
    # Pre-computed pixel-variance trace from session.camera0 (face camera),
    # parameterised exactly like Speed: bin to 100 ms, z-score per session
    # on motion bins, then evaluate 5 raised-cosine bases tiling the
    # z-scored value range. Per-bin behavioural covariate, NOT a history-
    # shaped lag basis (the lag-basis design from prompt 05 was rejected).
    # Captures non-linear ME-tuning curves the same way Speed/TF tuning
    # curves are captured.
    #
    # The earlier prompt-05 spec proposed convolving ME with a temporal
    # lag basis (5 raised cosines × 300 ms window). That parameterisation
    # was rejected on 2026-04-30: ME is a per-bin behavioural input, not
    # a history-shaped signal. The Speed-style value-axis basis lets the
    # marginal ME tuning curve be inspectable without smearing it across
    # past bins.
    #
    # Range default (-2.0, 3.0) z-score units (slightly right-skewed to
    # cover the high-ME tail of whisking/grooming bouts). Widen if the
    # empirical per-session ME distribution clips at the boundaries.
    # Switch to camera1 (body cam) via motion_energy_camera="camera1"
    # without re-plumbing.
    #
    # CLI flag: --no-me-face turns off the candidate; --motion-energy-camera
    # picks camera0 vs camera1.
    n_me_face_bases: int = 5
    me_face_range: tuple[float, float] = (-2.0, 3.0)
    motion_energy_camera: str = "camera0"
    # Master gate for the ME_face Phase-1 candidate. True (default):
    # build B_me_face whenever camera0 is present. False (--no-me-face):
    # skip ME_face construction entirely; "ME_face" is dropped from
    # remaining candidates in forward_select. ME_face_x_Speed Phase-2
    # eligibility falls through automatically.
    include_me_face: bool = True

    # --- Acceleration (signed translation acceleration; binned in time_binning
    # as the 'acceleration' column, 0 in V/stationary). Value-axis Phase-1
    # candidate like Speed/ME when include_acceleration is True AND "Acceleration"
    # is in main_effects. Built per cluster from df['acceleration'] (z-scored,
    # raised_cosine_basis_linear over accel_range). Default OFF — existing runs
    # are unaffected. ---
    include_acceleration: bool = False
    n_accel_bases: int = 5
    accel_range: tuple[float, float] = (-3.0, 3.0)

    # --- GLM fitting ---
    # Ridge on all non-intercept columns. Tuned by held-out cv_bps on
    # probe 243 (no-ME no-history config, 33 retained clusters,
    # 2026-05-07 sweep at λ ∈ {1e-8, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10}):
    # mean Selected cv_bps maximises at λ=1.0 (mean = -3.7962, vs the
    # previous default 1e-3 = -3.8460, +0.05 cv_bps mean improvement).
    # The improvement is tail-driven (median is flat at -3.337 across
    # all λ); a few bad-fit clusters benefit from the stronger shrinkage,
    # the median cluster is unaffected. Selection sets are essentially
    # unchanged across λ. Plot:
    # ~/local_data/motion_clouds/figures/glm/exploration/ridge_cv_sweep.{pdf,png}.
    #
    # Earlier history: λ=1e-3 was originally chosen 2026-04-22 (commit
    # 11226a8) for MATLAB-parity calibration of per-cluster tuning-curve
    # Pearson — pinned Python's IRLS rotation close to MATLAB glmnet's
    # lambda_1se rotation. MATLAB parity was retired as a production
    # gate 2026-04-29; the value was inherited until the cv-bps re-tune
    # on 2026-05-07. FullInteraction keeps its own (larger) lambda
    # because p approaches n there.
    lambda_ridge: float = 1.0
    full_interaction_lambda: float = 1.0
    lambda_ridge_min: float = 1e-6
    # Spike-history smoothness prior (ASD-style, 2026-06-05). None → plain
    # isotropic ridge on the History block (historical behaviour). When set,
    # the History coefficients are penalised by the squared second difference
    # of the reconstructed lag-space filter (lambda * B.T D2.T D2 B + ridge
    # floor), so fine-bin history filters stay smooth instead of oscillating.
    # Only bites at >=3 lag bins (e.g. 20 ms); a no-op at the 100 ms / 2-lag
    # default. See rc2_glm.penalty.build_penalty_matrix.
    history_smooth_lambda: float | None = None
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
    # Re-label profile_id from the recorded velocity trajectory (KMeans-2 on
    # onset-aligned |velocity|) instead of the stimulus trial-order halving.
    # The halving (StimulusLookup.trial_profile_id = trial_id <= midpoint)
    # does NOT track the two reproduced velocity trajectories — they are
    # interleaved across trial_id (verified 2026-06-08: dip/flat split
    # 18/18/18/18 against profile_id). The recorded velocity is the ground
    # truth of the trajectory each trial actually had. Default True (the
    # fix); set False to reproduce the legacy trial-halving profile_id.
    # See io._assign_profile_ids_by_velocity.
    profile_from_velocity: bool = True

    # --- Split-by-condition fitting (2026-06-05) ---
    # When set to one of {"V", "T_Vstatic", "VT"}, restrict the fit to
    # the trials of that single condition: that condition's motion bins
    # PLUS their stationary prelude (the stationary rows carry the trial's
    # trial_id but are tagged condition='stationary'; see time_binning
    # bin_trial). Keeping the prelude preserves the intercept+onset
    # "baseline vs motion" Null model while dropping the other two
    # conditions' trials entirely. The caller is responsible for also
    # restricting ``main_effects`` / ``interactions`` to the condition's
    # non-degenerate regressors (V: TF/SF/OR — Speed is identically 0;
    # T_Vstatic: Speed only — TF is 0 and SF/OR are NaN; VT: all four).
    # None (default) = the standard pooled-conditions fit. Driver:
    # scripts/run_glm_split_by_condition.py
    # (figures/glm/current_splitted_by_condition/).
    fit_condition: str | None = None

    # --- Forward selection ---
    # ME_face joined main_effects 2026-04-30 (prompt 06) as a Speed-style
    # value-axis candidate. ME_face_x_Speed joined interactions same day
    # — face-motion vs locomotion-speed are correlated 0.3–0.7 (Musall
    # 2019, Stringer 2019), so the interaction is the obvious one.
    main_effects: tuple[str, ...] = ("Speed", "TF", "SF", "OR", "ME_face")
    interactions: tuple[str, ...] = (
        "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
        "TF_x_SF", "TF_x_OR", "SF_x_OR",
        "ME_face_x_Speed",
    )
    delta_bps_threshold: float = 0.005

    # --- Forward-selection robustness (multi-seed admission) ---
    # Admit a candidate at a forward-selection round iff its Δ cv_bps
    # clears delta_bps_threshold in at least selection_threshold_count of
    # n_selection_seeds independent cv-fold partitions (built by
    # make_trial_folds with seeds 0..n_selection_seeds-1). Defaults below
    # reduce to today's single-seed behaviour for back-compat. Production
    # reruns 2026-05-08 set n_selection_seeds=10, selection_threshold_count=7
    # — the 7/10 admission rule was chosen against the 2026-05-07 seed
    # sweep where cluster 376 (TF on the knife edge) was 9/10 above
    # threshold and cluster 377 (SF/OR straddle) was on the boundary.
    # The "best" candidate of a passing round is the admitted candidate
    # with the highest mean Δ across the n_selection_seeds partitions.
    n_selection_seeds: int = 1
    selection_threshold_count: int = 1

    # --- Reference levels (matches MATLAB sf_levels / or_levels) ---
    sf_levels: tuple[float, ...] = (0.003, 0.006, 0.012)
    or_levels: tuple[float, ...] = field(
        default_factory=lambda: (-0.7853981633974483, 0.0, 0.7853981633974483, 1.5707963267948966)
    )  # -π/4, 0, π/4, π/2 (sorted)

    # --- RF-local SF/OR (goggles, opt-in; 2026-06-12) ---
    # How the SF and OR regressors are built:
    #   "tokens" (default): the categorical stimulus generation tokens
    #     (sf_levels / or_levels), reference-coded dummies — the MATLAB-
    #     parity behaviour. ALL screens and existing goggles token runs use
    #     this; leaving the default here keeps them byte-identical.
    #   "rf_local": the per-cluster, per-bin *local* SF (cpd) and OR (deg)
    #     read from each cluster's receptive field on the motion cloud
    #     (Gabor best-match extraction, gabor_extract_gpu.py → the cohort
    #     parquet). SF becomes a continuous log-Weber raised-cosine basis and
    #     OR a circular (π-periodic) basis — like Speed/TF, NOT dummies. Only
    #     the goggles cohort has RFs + cloud frames, so this is goggles-only.
    #     Requires rf_sf_or_parquet_dir; clusters without an RF fall out of
    #     the cohort (see rc2_glm.rf_sf_or + pipeline cohort intersection).
    # The continuous-basis code path (B_sf / B_or) only activates under
    # "rf_local"; "tokens" never touches it.
    sf_or_source: str = "tokens"
    # Directory of per-cloud parquet files written by gabor_extract_gpu.py
    # (cols: probe, cluster, rf_type, cx, cy, frame, sf_cpd, or_deg,
    # concentration, edge, cloud). Default points at the local cohort mirror;
    # only read when sf_or_source == "rf_local".
    rf_sf_or_parquet_dir: str = (
        "~/local_data/motion_clouds/saved_goggles/_extract/cohort"
    )
    # SF continuous basis: log-Weber raised cosines over the observed local-SF
    # range (extraction spans ~0.021–0.072 cpd; widen for headroom). OR
    # continuous basis: circular bumps over [0, 180) deg (orientation is
    # π-periodic — 0° ≡ 180°). Modest knot counts: per cluster the local SF
    # clusters near ~3 tokens and OR near ~4, so few knots avoid over-
    # parameterising a near-discrete variable.
    n_sf_bases: int = 4
    # Covers the pooled per-bin local-SF span (observed ~0.023–0.139 cpd across
    # the 36 clouds once the velocity-locked frame fluctuation is included), so
    # the high-SF tail isn't clipped onto the top knot.
    sf_cpd_range: tuple[float, float] = (0.02, 0.15)
    n_or_bases: int = 6
    # Drop an RF whose median Gabor concentration (peak/mean energy, a
    # reliability proxy) is below this. 0.0 = keep all; the goggle cohort
    # sits ~3–5 so this is non-binding, but it guards future noisier RFs.
    rf_min_concentration: float = 0.0
    # rf_local "_all" mode: instead of EXCLUDING clusters without a clean RF,
    # keep them in the cohort and give their SF/OR the per-cloud NOMINAL value
    # (cohort-mean SF cpd + circular-mean OR deg over all RF clusters/frames for
    # that cloud) — constant per trial, the continuous "imitate the dummy
    # values" stand-in. RF clusters are unchanged. Only read when rf_local.
    # Default False = the RF-only cohort (excludes no-RF clusters), so existing
    # rf_local runs are byte-identical.
    rf_sf_or_nominal_fallback: bool = False

    # --- Prefilter ---
    # The stationary-vs-motion Wilcoxon is now a DIAGNOSTIC, not the default
    # selection gate: the whole selected cohort is fit, and the prefilter table
    # is still computed/written so the motion-responsive funnel stays on record
    # (see pipeline: compute is decoupled from gating). A spike-count floor will
    # become the quality gate in its place. Set True to restore gating (the
    # cohort = should_run_glm rows only).
    apply_prefilter: bool = False
    prefilter_seed: int = 0

    # --- Spike-count quality floor (the gate that replaces the prefilter; 2026-06-15) ---
    # A cluster is dropped from the cohort BEFORE fitting unless, in the data the
    # run actually fits (per-condition for the split runs), it has at least
    # ``min_spikes_floor`` total spikes AND fires in at least
    # ``min_trial_occupancy`` of its trials. Below this the per-spike cv-bps is
    # too noisy to interpret (sim: SD ~0.17 bps at 50 spikes, ~0.35 at 15; the
    # cluster-105 goggles outlier had 30 spikes — see
    # project_motion_clouds_cvbps_stability_defaults). Defaults ON for ALL runs;
    # set both to 0 to reproduce the unfiltered legacy cohort. This is a
    # PRINCIPLED EXCLUSION (too few spikes to estimate the quantity), not a clip
    # on an anomalous value. Gate: prefilter.passes_spike_floor.
    min_spikes_floor: int = 50
    min_trial_occupancy: float = 0.5

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
    "ME_face_x_Speed": ("ME_face", "Speed"),
    # Acceleration / ME_face full-pairwise set (2026-06-15): every pairwise
    # interaction on the value/visual predictors, so each condition can carry
    # all interactions on its non-degenerate candidates (goggles rf_sfor _all
    # + the V / T_Vstatic splits). Acceleration is a value-axis basis like
    # Speed/ME_face; SF/OR are categorical (tokens) or continuous (rf_local) —
    # the design_matrix branches handle both representations.
    "Speed_x_Acceleration": ("Speed", "Acceleration"),
    "TF_x_ME_face": ("TF", "ME_face"),
    "TF_x_Acceleration": ("TF", "Acceleration"),
    "SF_x_ME_face": ("SF", "ME_face"),
    "SF_x_Acceleration": ("SF", "Acceleration"),
    "OR_x_ME_face": ("OR", "ME_face"),
    "OR_x_Acceleration": ("OR", "Acceleration"),
    "ME_face_x_Acceleration": ("ME_face", "Acceleration"),
}
