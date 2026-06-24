"""Tests for the V→VT tuning-significance transition classifier.

``tuning_transition_category`` is the pure logic behind the cohort transition
figure (``plot_tuning_transitions``): given the tuned/not-tuned call in each
condition plus the selected best_model, it returns the 1–4 transition class
(or 0 for cells tuned in neither condition, which the figure drops). Only this
branching is tested; the rendering follows the script's untested convention.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from make_fens_poster_figures import tuning_transition_category  # noqa: E402


def test_tuned_v_only_is_cat1() -> None:
    assert tuning_transition_category(True, False, "cubic", "sigmoid") == 1
    # model identity is irrelevant once VT is not tuned
    assert tuning_transition_category(True, False, "cubic", "cubic") == 1


def test_tuned_vt_only_is_cat2() -> None:
    assert tuning_transition_category(False, True, "cubic", "sigmoid") == 2
    assert tuning_transition_category(False, True, "cubic", "cubic") == 2


def test_both_tuned_different_model_is_cat3() -> None:
    assert tuning_transition_category(True, True, "cubic", "sigmoid") == 3


def test_both_tuned_same_model_is_cat4() -> None:
    assert tuning_transition_category(True, True, "cubic", "cubic") == 4
    assert tuning_transition_category(True, True, "vonmises_180", "vonmises_180") == 4


def test_neither_tuned_is_cat0() -> None:
    assert tuning_transition_category(False, False, "cubic", "cubic") == 0
    assert tuning_transition_category(False, False, "cubic", "sigmoid") == 0


def test_truthy_inputs_are_coerced() -> None:
    # pandas hands in numpy bools / 0-1 ints; the classifier must coerce.
    assert tuning_transition_category(1, 0, "a", "b") == 1
    assert tuning_transition_category(0, 1, "a", "b") == 2


def test_model_identity_is_exact_string_match() -> None:
    # different float param strings that are textually distinct → diff model
    assert tuning_transition_category(True, True, "gaussian", "asym_gaussian") == 3
    # numeric-equal but string-equal stays cat4
    assert tuning_transition_category(True, True, 1.0, 1.0) == 4
