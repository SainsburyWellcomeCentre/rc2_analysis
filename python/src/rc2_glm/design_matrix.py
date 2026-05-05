"""Assemble the GLM design matrix.

Two builders, mirroring the two MATLAB functions:

- `assemble_design_matrix(model_label, ...)` — fixed model labels
  ('Null', 'M0', 'Additive', 'FullInteraction', 'Additive_no_*'). See
  MATLAB `assemble_design_matrix` (line 5572).
- `assemble_design_matrix_selected(selected_vars, ...)` — driven by a
  list of selected variable / interaction names; mirrors MATLAB
  `assemble_design_matrix_selected` (line 6147). This is the version
  used by the forward-selection loop.

Both always include intercept + onset bases. SF / OR are reference-coded
dummies (first level dropped). Zero-variance columns (except the
intercept) are removed during training; pass non-empty `sf_ref_levels` /
`or_ref_levels` to disable that pruning (prediction mode).
"""

from __future__ import annotations

import numpy as np


def assemble_design_matrix_selected(
    B_speed: np.ndarray,
    B_tf: np.ndarray,
    B_onset: np.ndarray,
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    selected_vars: list[str],
    sf_ref_levels: np.ndarray | None = None,
    or_ref_levels: np.ndarray | None = None,
    *,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    include_onset_kernel: bool = True,
) -> tuple[np.ndarray, list[str]]:
    """Assemble the design matrix from selected variable names.

    Keyword-only arguments:

    - ``B_history``: per-cluster spike-history features ``(n_bins, n_history_bases)``.
      When provided AND ``"History"`` is in ``selected_vars``, append the
      history columns. ``None`` → no history term regardless of selection.
      (Added 2026-04-28 for prompt 03; default flipped off 2026-04-30.)
    - ``B_me_face``: per-bin face motion energy basis ``(n_bins, n_me_face_bases)``.
      Speed-style raised-cosine basis evaluated at each bin's z-scored ME
      value. When provided AND ``"ME_face"`` is in ``selected_vars``,
      append the ME_face columns; ``"ME_face_x_Speed"`` builds the
      row-wise outer product with ``B_speed`` for the interaction.
      ``None`` → no ME term. (Added 2026-04-30 for prompt 06.)
    - ``include_onset_kernel``: when False, the Onset basis is OMITTED
      from the design matrix entirely. Used by the prompt-03 ablation
      experiment. Default True for parity-preserving baseline behaviour.
    """
    n = B_speed.shape[0]
    n_speed_b = B_speed.shape[1]
    n_tf_b = B_tf.shape[1]

    sf_levels, or_levels = _resolve_reference_levels(
        sf_vals, or_vals, sf_ref_levels, or_ref_levels
    )

    cols: list[np.ndarray] = [np.ones((n, 1))]
    names: list[str] = ["Intercept"]

    if include_onset_kernel and B_onset is not None and B_onset.shape[1] > 0:
        cols.append(B_onset)
        names += [f"Onset_{i + 1}" for i in range(B_onset.shape[1])]

    selected = set(selected_vars)

    if "History" in selected and B_history is not None and B_history.shape[1] > 0:
        cols.append(B_history)
        names += [f"History_{i + 1}" for i in range(B_history.shape[1])]

    if "Speed" in selected:
        cols.append(B_speed)
        names += [f"Speed_{i + 1}" for i in range(n_speed_b)]
    if "TF" in selected:
        cols.append(B_tf)
        names += [f"TF_{i + 1}" for i in range(n_tf_b)]
    if "SF" in selected:
        for level in sf_levels:
            cols.append(_dummy(sf_vals, level).reshape(-1, 1))
            names.append(f"SF_{level:.4f}")
    if "OR" in selected:
        for level in or_levels:
            cols.append(_dummy(or_vals, level).reshape(-1, 1))
            names.append(f"OR_{level:.3f}")
    if "ME_face" in selected and B_me_face is not None and B_me_face.shape[1] > 0:
        cols.append(B_me_face)
        names += [f"ME_face_{i + 1}" for i in range(B_me_face.shape[1])]

    if "Speed_x_TF" in selected:
        for si in range(n_speed_b):
            for ti in range(n_tf_b):
                cols.append((B_speed[:, si] * B_tf[:, ti]).reshape(-1, 1))
                names.append(f"Spd{si + 1}_x_TF{ti + 1}")
    if "Speed_x_SF" in selected:
        for si in range(n_speed_b):
            for level in sf_levels:
                d = _dummy(sf_vals, level)
                cols.append((B_speed[:, si] * d).reshape(-1, 1))
                names.append(f"Spd{si + 1}_x_SF{level:.4f}")
    if "Speed_x_OR" in selected:
        for si in range(n_speed_b):
            for level in or_levels:
                d = _dummy(or_vals, level)
                cols.append((B_speed[:, si] * d).reshape(-1, 1))
                names.append(f"Spd{si + 1}_x_OR{level:.3f}")
    if "TF_x_SF" in selected:
        for ti in range(n_tf_b):
            for level in sf_levels:
                d = _dummy(sf_vals, level)
                cols.append((B_tf[:, ti] * d).reshape(-1, 1))
                names.append(f"TF{ti + 1}_x_SF{level:.4f}")
    if "TF_x_OR" in selected:
        for ti in range(n_tf_b):
            for level in or_levels:
                d = _dummy(or_vals, level)
                cols.append((B_tf[:, ti] * d).reshape(-1, 1))
                names.append(f"TF{ti + 1}_x_OR{level:.3f}")
    if "SF_x_OR" in selected:
        for sf_level in sf_levels:
            d_sf = _dummy(sf_vals, sf_level)
            for or_level in or_levels:
                d_or = _dummy(or_vals, or_level)
                cols.append((d_sf * d_or).reshape(-1, 1))
                names.append(f"SF{sf_level:.4f}_x_OR{or_level:.3f}")
    if "ME_face_x_Speed" in selected and B_me_face is not None and B_me_face.shape[1] > 0:
        # 5 ME bases × 5 Speed bases = 25 interaction columns.
        # Same row-wise product convention as Speed_x_TF above.
        for mi in range(B_me_face.shape[1]):
            for si in range(n_speed_b):
                cols.append((B_me_face[:, mi] * B_speed[:, si]).reshape(-1, 1))
                names.append(f"MEf{mi + 1}_x_Spd{si + 1}")

    X = np.hstack(cols).astype(np.float64, copy=False)

    is_prediction = sf_ref_levels is not None or or_ref_levels is not None
    if not is_prediction:
        X, names = _drop_zero_variance(X, names)
    return X, names


def assemble_design_matrix(
    B_speed: np.ndarray,
    B_tf: np.ndarray,
    B_onset: np.ndarray,
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    model_label: str,
    sf_ref_levels: np.ndarray | None = None,
    or_ref_levels: np.ndarray | None = None,
    *,
    B_history: np.ndarray | None = None,
    B_me_face: np.ndarray | None = None,
    include_onset_kernel: bool = True,
) -> tuple[np.ndarray, list[str]]:
    """Build a fixed model matrix labelled by `model_label`.

    See ``assemble_design_matrix_selected`` for the History / ME_face /
    onset-toggle keyword args. ``M0`` and ``Null`` ignore History and
    ME_face (both are Phase-1 candidates, not in the always-on baseline).
    """
    common_kwargs = dict(
        B_history=B_history,
        B_me_face=B_me_face,
        include_onset_kernel=include_onset_kernel,
    )
    if model_label == "M0":
        n = B_speed.shape[0]
        return np.ones((n, 1)), ["Intercept"]
    if model_label == "Null":
        return assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals,
            [], sf_ref_levels, or_ref_levels,
            **common_kwargs,
        )
    if model_label == "Additive":
        return assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals,
            ["Speed", "TF", "SF", "OR"], sf_ref_levels, or_ref_levels,
            **common_kwargs,
        )
    if model_label == "FullInteraction":
        # ME_face + ME_face_x_Speed added 2026-04-30 (prompt 06): keep
        # FullInteraction as the true ceiling so a Selected model that
        # picks ME_face_x_Speed can never out-cv-bps FullInteraction.
        # Branches in assemble_design_matrix_selected guard on
        # B_me_face is not None, so legacy callers without a B_me_face
        # arg still produce the pre-2026-04-30 design.
        return assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals,
            ["Speed", "TF", "SF", "OR", "ME_face",
             "Speed_x_TF", "Speed_x_SF", "Speed_x_OR",
             "TF_x_SF", "TF_x_OR", "SF_x_OR", "ME_face_x_Speed"],
            sf_ref_levels, or_ref_levels,
            **common_kwargs,
        )
    if model_label.startswith("Additive_no_"):
        drop = model_label.removeprefix("Additive_no_")
        keep = [v for v in ("Speed", "TF", "SF", "OR") if v != drop]
        return assemble_design_matrix_selected(
            B_speed, B_tf, B_onset, sf_vals, or_vals,
            keep, sf_ref_levels, or_ref_levels,
            **common_kwargs,
        )
    raise ValueError(f"Unknown model label: {model_label}")


def _dummy(vals: np.ndarray, level: float) -> np.ndarray:
    out = (vals == level).astype(np.float64)
    out[np.isnan(vals)] = 0.0
    return out


def _resolve_reference_levels(
    sf_vals: np.ndarray,
    or_vals: np.ndarray,
    sf_ref: np.ndarray | None,
    or_ref: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    if sf_ref is not None and len(sf_ref) > 0:
        sf_unique = np.asarray(sf_ref, dtype=np.float64).ravel()
    else:
        sf_unique = np.sort(np.unique(sf_vals[~np.isnan(sf_vals)]))
    sf_levels = sf_unique[1:] if sf_unique.size >= 1 else sf_unique

    if or_ref is not None and len(or_ref) > 0:
        or_unique = np.asarray(or_ref, dtype=np.float64).ravel()
    else:
        or_unique = np.sort(np.unique(or_vals[~np.isnan(or_vals)]))
    or_levels = or_unique[1:] if or_unique.size >= 1 else or_unique
    return sf_levels, or_levels


def _drop_zero_variance(
    X: np.ndarray, names: list[str], tol: float = 1e-10
) -> tuple[np.ndarray, list[str]]:
    if X.shape[1] <= 1:
        return X, names
    keep = np.ones(X.shape[1], dtype=bool)
    keep[1:] = X[:, 1:].var(axis=0) > tol
    if keep.all():
        return X, names
    return X[:, keep], [n for n, k in zip(names, keep) if k]
