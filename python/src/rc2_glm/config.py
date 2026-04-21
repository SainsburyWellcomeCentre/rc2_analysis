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
    lambda_ridge: float = 0.0
    full_interaction_lambda: float = 1.0
    lambda_ridge_min: float = 1e-6
    irls_max_iter: int = 100
    irls_tol: float = 1e-8
    eta_clip: float = 20.0
    mu_floor: float = 1e-10

    # --- Cross-validation ---
    n_folds: int = 5
    cv_seed: int = 0

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


INTERACTION_PARENTS: dict[str, tuple[str, str]] = {
    "Speed_x_TF": ("Speed", "TF"),
    "Speed_x_SF": ("Speed", "SF"),
    "Speed_x_OR": ("Speed", "OR"),
    "TF_x_SF": ("TF", "SF"),
    "TF_x_OR": ("TF", "OR"),
    "SF_x_OR": ("SF", "OR"),
}
