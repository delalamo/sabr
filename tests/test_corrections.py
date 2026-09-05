import warnings

import numpy as np
import pytest

from sabr.corrections import (
    apply_corrections,
    cdr_columns,
    correct_cdr_loop,
    correct_de_loop,
    de_loop_positions,
)
from sabr.numbering import alignment_to_states, number_alignment


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


@pytest.mark.parametrize(
    ("residue_count", "expected"),
    [
        (0, []),
        (1, [80]),
        (2, [80, 84]),
        (3, [80, 83, 84]),
        (4, [80, 82, 83, 84]),
        (5, [80, 81, 82, 83, 84]),
        (6, [80, 81, 82, 82, 83, 84]),
    ],
)
def test_de_loop_positions_follow_imgt_79_85_convention(
    residue_count, expected
):
    assert de_loop_positions(residue_count) == expected


@pytest.mark.parametrize(
    ("expected_columns", "expected_82_states"),
    [
        ([], [(82, "d")]),
        ([79], [(82, "d")]),
        ([79, 83], [(82, "d")]),
        ([79, 82, 83], [(82, "d")]),
        ([79, 81, 82, 83], [(82, "m")]),
        ([79, 80, 81, 82, 83], [(82, "m")]),
        ([79, 80, 81, None, 82, 83], [(82, "m"), (82, "i")]),
    ],
)
def test_de_loop_is_rebuilt_between_79_and_85(
    expected_columns, expected_82_states
):
    residue_count = len(expected_columns)
    alignment = np.zeros((residue_count + 2, 128), dtype=int)
    alignment[0, 78] = 1
    alignment[-1, 84] = 1
    for row in range(1, residue_count + 1):
        alignment[row, 78] = 1

    corrected = correct_de_loop(alignment)
    expected = np.zeros_like(alignment)
    expected[0, 78] = 1
    expected[-1, 84] = 1
    for row, column in enumerate(expected_columns, start=1):
        if column is not None:
            expected[row, column] = 1
    np.testing.assert_array_equal(corrected, expected)

    states, _, _ = alignment_to_states(corrected)
    assert [state for state, _ in states if state[0] == 82] == (
        expected_82_states
    )


def test_six_residue_de_loop_numbers_the_insertion_as_82a():
    alignment = np.zeros((8, 128), dtype=int)
    alignment[0, 78] = 1
    alignment[-1, 84] = 1
    alignment[1:7, 78] = 1

    records = number_alignment(correct_de_loop(alignment), "A" * 8, "imgt", "H")

    assert [(number, code.strip()) for _, number, code, _ in records] == [
        (79, ""),
        (80, ""),
        (81, ""),
        (82, ""),
        (82, "A"),
        (83, ""),
        (84, ""),
        (85, ""),
    ]


def test_structural_gap_skips_de_loop_correction():
    alignment = np.zeros((7, 128), dtype=int)
    for row, position in enumerate(range(79, 86)):
        alignment[row, position - 1] = 1
    original = alignment.copy()

    with pytest.warns(UserWarning, match=r"Skipping DE loop \(79-85\)"):
        corrected = correct_de_loop(alignment, gap_indices=frozenset({1}))
    np.testing.assert_array_equal(corrected, original)


def test_standard_de_loop_keeps_the_learned_alignment():
    alignment = np.zeros((130, 128), dtype=int)
    for row in range(125):
        alignment[row, row] = 1
    original = alignment.copy()
    corrected = apply_corrections(alignment)
    np.testing.assert_array_equal(corrected, original)


@pytest.mark.parametrize(
    "region,start,end,anchors",
    [
        ("CDR1", 27, 38, (23, 40)),
        ("CDR2", 56, 65, (54, 67)),
        ("CDR3", 105, 117, (104, 118)),
    ],
)
@pytest.mark.parametrize("scenario", ("missing", "reversed", "same_row"))
def test_unusable_cdr_anchors_preserve_alignment(
    caplog, region, start, end, anchors, scenario
):
    alignment = np.zeros((30, 128), dtype=int)
    if scenario != "missing":
        alignment[20, anchors[0] - 1] = 1
        alignment[20 if scenario == "same_row" else 2, anchors[1] - 1] = 1
    original = alignment.copy()
    np.testing.assert_array_equal(
        correct_cdr_loop(alignment, region, start, end), original
    )
    assert f"Skipping {region}" in caplog.text


@pytest.mark.parametrize(
    "region,start,end,anchors",
    [("CDR1", 27, 38, (23, 40)), ("CDR2", 56, 65, (54, 67))],
)
def test_cdr_with_insufficient_framework_residues_is_unchanged(
    caplog, region, start, end, anchors
):
    alignment = np.zeros((2, 128), dtype=int)
    alignment[0, anchors[0] - 1] = 1
    alignment[1, anchors[1] - 1] = 1
    original = alignment.copy()
    np.testing.assert_array_equal(
        correct_cdr_loop(alignment, region, start, end), original
    )
    assert "not enough residues" in caplog.text


@pytest.mark.parametrize("offset", (-2, -1, 1, 2))
def test_nearby_anchors_still_rebuild_cdr(offset):
    alignment = np.zeros((18, 128), dtype=int)
    alignment[0, 22 + offset] = 1
    alignment[17, 39 + offset] = 1
    corrected = correct_cdr_loop(alignment, "CDR1", 27, 38)
    assert np.flatnonzero(corrected[1]).tolist() == [23]
    assert np.flatnonzero(corrected[16]).tolist() == [38]
    assert corrected[4:16, 26:38].trace() == 12


@pytest.mark.parametrize(
    "count,expected",
    [
        (0, []),
        (3, [26, 27, 37]),
        (12, list(range(26, 38))),
        (14, [*range(26, 32), None, None, *range(32, 38)]),
    ],
)
def test_cdr1_rebuilds_short_standard_and_overflowing_loops(count, expected):
    alignment = np.zeros((count + 6, 128), dtype=int)
    alignment[0, 22] = 1
    alignment[-1, 39] = 1
    corrected = correct_cdr_loop(alignment, "CDR1", 27, 38)
    for row, column in enumerate(expected, start=4):
        assert np.flatnonzero(corrected[row]).tolist() == (
            [] if column is None else [column]
        )
    assert [
        np.flatnonzero(corrected[i]).item() for i in (1, 2, 3, count + 4)
    ] == [23, 24, 25, 38]


@pytest.mark.parametrize(
    "region,start,end",
    [("CDR1", 22, 39), ("CDR2", 53, 66), ("CDR3", 103, 117), ("DE", 78, 84)],
)
@pytest.mark.parametrize(
    "boundary,skipped",
    [("before", False), ("first", True), ("last", True), ("after", False)],
)
def test_structural_gap_boundaries_skip_only_crossed_bonds(
    region, start, end, boundary, skipped
):
    alignment = np.eye(128, dtype=int)
    alignment[start + 1] = 0
    original = alignment.copy()
    gap = {"before": start - 1, "first": start, "last": end - 1, "after": end}[
        boundary
    ]
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        if region == "DE":
            result = correct_de_loop(alignment, frozenset({gap}))
        else:
            bounds = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}[
                region
            ]
            result = correct_cdr_loop(
                alignment, region, *bounds, frozenset({gap})
            )
    assert bool(captured) is skipped
    if skipped:
        np.testing.assert_array_equal(result, original)
    else:
        np.testing.assert_array_equal(result, np.eye(128, dtype=int))


def test_gap_in_cdr1_does_not_suppress_cdr2_repair():
    alignment = np.eye(128, dtype=int)
    alignment[30] = 0
    alignment[58] = 0
    with pytest.warns(UserWarning, match="Skipping CDR1") as captured:
        corrected = apply_corrections(alignment, frozenset({30}))
    assert len(captured) == 1
    assert not corrected[30].any()
    assert corrected[58, 58] == 1


@pytest.mark.parametrize("scenario", ("missing", "reversed"))
def test_unusable_de_anchors_preserve_alignment(caplog, scenario):
    alignment = np.zeros((3, 128), dtype=int)
    if scenario == "reversed":
        alignment[2, 78] = 1
        alignment[0, 84] = 1
    original = alignment.copy()
    np.testing.assert_array_equal(correct_de_loop(alignment), original)
    assert "Skipping DE loop" in caplog.text
