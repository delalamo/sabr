import warnings

import numpy as np
import pytest

from sabr.corrections import (
    apply_corrections,
    cdr_columns,
    correct_cdr_loop,
)


@pytest.mark.parametrize(
    ("residues", "positions"),
    [
        (5, [0, 1, 2, 11, 12]),
        (6, [0, 1, 2, 10, 11, 12]),
    ],
)
def test_cdr_gaps_follow_the_imgt_alternating_pattern(residues, positions):
    assert cdr_columns(residues, 13) == positions


def test_cdr_overflow_leaves_central_residues_as_insertions():
    assert cdr_columns(14, 13) == [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        None,
        7,
        8,
        9,
        10,
        11,
        12,
    ]


@pytest.mark.parametrize(
    ("name", "start", "end", "anchors"),
    [
        ("CDR1", 27, 38, (23, 40)),
        ("CDR2", 56, 65, (54, 67)),
        ("CDR3", 105, 117, (104, 118)),
    ],
)
def test_structural_gap_skips_each_affected_cdr(name, start, end, anchors):
    alignment = np.zeros((20, 128), dtype=int)
    alignment[2, anchors[0] - 1] = 1
    alignment[12, anchors[1] - 1] = 1
    original = alignment.copy()
    with pytest.warns(UserWarning, match=f"Skipping {name}.*structural gap"):
        corrected = correct_cdr_loop(
            alignment,
            name,
            start,
            end,
            gap_indices=frozenset({7}),
        )
    np.testing.assert_array_equal(corrected, original)


def test_multiple_structural_gaps_warn_for_only_the_affected_regions():
    alignment = np.eye(128, dtype=int)
    with pytest.warns(UserWarning) as captured:
        corrected = apply_corrections(
            alignment.copy(), gap_indices=frozenset({30, 58})
        )
    assert [str(item.message) for item in captured] == [
        "Skipping CDR1 deterministic correction: structural gap detected "
        "between rows 22 and 39; using the raw embedding alignment for "
        "this region.",
        "Skipping CDR2 deterministic correction: structural gap detected "
        "between rows 53 and 66; using the raw embedding alignment for "
        "this region.",
    ]
    np.testing.assert_array_equal(corrected, alignment)


def test_gap_outside_corrected_regions_emits_no_warning():
    alignment = np.eye(128, dtype=int)
    with warnings.catch_warnings(record=True) as captured:
        apply_corrections(alignment.copy(), gap_indices=frozenset({10}))
    assert captured == []


def test_de_loop_and_c_terminus_keep_the_learned_alignment():
    alignment = np.zeros((130, 128), dtype=int)
    for row in range(125):
        alignment[row, row] = 1
    original = alignment.copy()
    corrected = apply_corrections(alignment)
    np.testing.assert_array_equal(corrected, original)
