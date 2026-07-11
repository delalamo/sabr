import warnings

import numpy as np
import pytest

from sabr.corrections import (
    apply_corrections,
    cdr_columns,
    correct_c_terminus,
    correct_cdr_loop,
    correct_light_chain_de_loop,
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


def test_light_chain_de_loop_moves_81_and_82_to_83_and_84():
    alignment = np.zeros((4, 128), dtype=int)
    alignment[0, 80] = 1
    alignment[1, 81] = 1
    corrected = correct_light_chain_de_loop(alignment)
    assert corrected[0, 82] == 1
    assert corrected[1, 83] == 1
    assert corrected[:, 80:82].sum() == 0


def test_light_chain_de_loop_clears_sources_when_targets_are_occupied():
    alignment = np.zeros((4, 128), dtype=int)
    alignment[0, 80] = 1
    alignment[1, 81] = 1
    alignment[2, 82] = 1
    alignment[3, 83] = 1
    corrected = correct_light_chain_de_loop(alignment)
    assert corrected[:, 80:82].sum() == 0
    assert corrected[2, 82] == 1
    assert corrected[3, 83] == 1


def test_gap_outside_de_loop_does_not_skip_its_correction():
    alignment = np.zeros((90, 128), dtype=int)
    alignment[70, 78] = 1
    alignment[71, 79] = 1
    alignment[72, 80] = 1
    alignment[73, 81] = 1
    alignment[75, 83] = 1
    corrected = correct_light_chain_de_loop(
        alignment, gap_indices=frozenset({10})
    )
    assert corrected[72, 82] == 1
    assert corrected[73, 81] == 0
    assert corrected[75, 83] == 1


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
            alignment.copy(), "H", gap_indices=frozenset({30, 58})
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
        apply_corrections(alignment.copy(), "H", gap_indices=frozenset({10}))
    assert captured == []


def test_de_loop_gap_does_not_prevent_c_terminal_correction():
    alignment = np.zeros((130, 128), dtype=int)
    alignment[70, 78] = 1
    alignment[71, 79] = 1
    alignment[72, 80] = 1
    alignment[73, 81] = 1
    alignment[74, 82] = 1
    alignment[75, 83] = 1
    alignment[124, 124] = 1
    with pytest.warns(UserWarning, match="FR3/DE loop"):
        corrected = apply_corrections(
            alignment, "K", gap_indices=frozenset({72})
        )
    assert corrected[72, 80] == 1
    assert corrected[73, 81] == 1
    assert corrected[125, 125] == 1
    assert corrected[126, 126] == 1
    assert corrected[127, 127] == 1


def test_c_terminal_assigns_available_positions_and_stops_at_128():
    alignment = np.zeros((130, 128), dtype=int)
    for row in range(125):
        alignment[row, row] = 1
    corrected = correct_c_terminus(alignment)
    assert corrected[125, 125] == 1
    assert corrected[126, 126] == 1
    assert corrected[127, 127] == 1
    assert corrected[128:].sum() == 0


def test_c_terminal_does_nothing_before_imgt_125():
    alignment = np.zeros((100, 128), dtype=int)
    for row in range(90):
        alignment[row, row] = 1
    original = alignment.copy()
    np.testing.assert_array_equal(correct_c_terminus(alignment), original)


def test_c_terminal_does_nothing_when_every_row_is_assigned():
    alignment = np.zeros((120, 128), dtype=int)
    for row in range(120):
        alignment[row, row] = 1
    original = alignment.copy()
    np.testing.assert_array_equal(correct_c_terminus(alignment), original)


def test_c_terminal_assigns_one_trailing_residue():
    alignment = np.zeros((118, 128), dtype=int)
    for row in range(117):
        alignment[row, row] = 1
    alignment[116] = 0
    alignment[116, 125] = 1
    corrected = correct_c_terminus(alignment)
    assert corrected[117, 126] == 1


def test_c_terminal_handles_an_empty_alignment():
    alignment = np.zeros((100, 128), dtype=int)
    np.testing.assert_array_equal(correct_c_terminus(alignment), alignment)
