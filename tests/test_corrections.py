import numpy as np
import pytest

import sabr.corrections as corrections
from sabr.corrections import (
    apply_corrections,
    correct_c_terminus,
    correct_cdr_loop,
    correct_fr3_alignment,
    correct_gap_numbering,
)


@pytest.mark.parametrize(
    ("residues", "positions"),
    [
        (5, [0, 1, 2, 11, 12]),
        (6, [0, 1, 2, 10, 11, 12]),
    ],
)
def test_cdr_gaps_follow_the_imgt_alternating_pattern(residues, positions):
    corrected = correct_gap_numbering(np.zeros((residues, 13), dtype=int))
    assert corrected.sum() == residues
    assert list(np.where(corrected == 1)[1]) == positions
    assert corrected[0, 0] == 1
    assert corrected[-1, -1] == 1


def test_light_chain_de_loop_moves_81_and_82_to_83_and_84():
    alignment = np.zeros((4, 128), dtype=int)
    alignment[0, 80] = 1
    alignment[1, 81] = 1
    corrected = correct_fr3_alignment(alignment)
    assert corrected[0, 82] == 1
    assert corrected[1, 83] == 1
    assert corrected[:, 80:82].sum() == 0


def test_light_chain_de_loop_clears_sources_when_targets_are_occupied():
    alignment = np.zeros((4, 128), dtype=int)
    alignment[0, 80] = 1
    alignment[1, 81] = 1
    alignment[2, 82] = 1
    alignment[3, 83] = 1
    corrected = correct_fr3_alignment(alignment)
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
    corrected = correct_fr3_alignment(alignment, gap_indices=frozenset({10}))
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


def test_fr1_position_10_correction_is_completely_removed():
    assert not hasattr(corrections, "correct_fr1_alignment")
    assert "FR1" not in corrections.__doc__
