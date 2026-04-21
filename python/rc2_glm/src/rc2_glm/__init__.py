"""Python replication of MATLAB GLM (v10 Hierarchical Forward Selection).

Pipeline:
    1. io.load_trial_data — wraps rc2_formatted_data_reader to return per-trial
       structured records.
    2. time_binning.bin_cluster — 100 ms time-binned table per cluster.
    3. basis.raised_cosine_basis / onset_kernel_basis — basis matrices.
    4. design_matrix.assemble_design_matrix_selected — assemble X for a
       chosen list of selected variables.
    5. fitting.fit_poisson_glm — IRLS Poisson fit (matches MATLAB).
    6. cross_validation.cross_validate_glm — trial-level k-fold + bps.
    7. forward_selection.forward_select — Hardcastle-style 2-phase selection.
"""

from rc2_glm import basis, config, design_matrix, fitting, io, masks_helpers
from rc2_glm import (
    cross_validation,
    forward_selection,
    pipeline,
    prefilter,
    time_binning,
)
from rc2_glm.config import GLMConfig

__all__ = [
    "GLMConfig",
    "basis",
    "config",
    "cross_validation",
    "design_matrix",
    "fitting",
    "forward_selection",
    "io",
    "masks_helpers",
    "pipeline",
    "prefilter",
    "time_binning",
]
