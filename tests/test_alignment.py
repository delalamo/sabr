from pathlib import Path

import numpy as np
import pytest

from sabr import constants
from sabr.alignment import (
    _affine_gap_penalty,
    _affine_score,
    _align_reference,
    _alignment_path,
    _concatenate_reference,
    _terminal_gap_penalty,
    _validate_multidomain_alignment,
    align,
    load_gap_penalties,
    load_references,
)
from sabr.api import _normalize_chain_type
from sabr.numbering import alignment_to_states

DATA = Path(__file__).parent / "data"
ASSETS = Path(__file__).parents[1] / "src" / "sabr" / "assets"


@pytest.mark.parametrize(
    ("mode", "noise_level"),
    (("sabr", 0.0), ("softalign", 2.0)),
)
def test_scfv_references_are_exact_ordered_concatenations(mode, noise_level):
    references = load_references(noise_level, mode, scfv=True)
    assert tuple(references) == (
        *constants.CHAIN_TYPES,
        *constants.SCFV_CHAIN_TYPES,
    )
    for representation in constants.SCFV_CHAIN_TYPES:
        first_type, second_type = representation
        embeddings, positions = references[representation]
        first_embeddings, first_positions = references[first_type]
        second_embeddings, second_positions = references[second_type]
        np.testing.assert_array_equal(
            embeddings,
            np.concatenate((first_embeddings, second_embeddings), axis=0),
        )
        assert positions == (
            *first_positions,
            *(
                position + constants.IMGT_MAX_POSITION
                for position in second_positions
            ),
        )
        assert not embeddings.flags.writeable


def test_scfv_alignment_allows_an_unassigned_linker_between_domains():
    alignment = np.zeros((5, 256), dtype=int)
    alignment[0, 0] = 1
    alignment[1, 1] = 1
    alignment[3, 128] = 1
    alignment[4, 129] = 1
    _validate_multidomain_alignment(alignment, "HK")


def test_higher_order_reference_is_an_ordered_concatenation():
    references = load_references(0.0)
    embeddings, positions = _concatenate_reference(references, "HHK")
    expected_embeddings = np.concatenate(
        (
            references["H"][0],
            references["H"][0],
            references["K"][0],
        ),
        axis=0,
    )
    np.testing.assert_array_equal(embeddings, expected_embeddings)
    assert positions == tuple(
        position + domain_index * constants.IMGT_MAX_POSITION
        for domain_index, domain_type in enumerate("HHK")
        for position in references[domain_type][1]
    )
    assert not embeddings.flags.writeable


@pytest.mark.parametrize("linker_length", (1, 3))
def test_scfv_reference_boundary_has_no_affine_linker_penalty(linker_length):
    query_length = linker_length + 2
    similarities = np.full((query_length, 4), -100.0, dtype=np.float32)
    similarities[:, (0, 3)] = 0.0
    similarities[0, 1] = 10.0
    similarities[-1, 2] = 10.0
    lengths = np.asarray((query_length, similarities.shape[1]))
    gap_extend = -1.0
    gap_open = -3.0

    penalized = _affine_score(
        similarities,
        lengths,
        1e-4,
        gap_extend,
        gap_open,
    )
    free = _affine_score(
        similarities,
        lengths,
        1e-4,
        gap_extend,
        gap_open,
        free_gap_boundary=1,
    )

    expected_penalty = gap_open + (linker_length - 1) * gap_extend
    assert float(penalized) == pytest.approx(20.0 + expected_penalty)
    assert float(free) == pytest.approx(20.0)


def test_higher_order_reference_makes_every_linker_boundary_free():
    similarities = np.full((5, 5), -100.0, dtype=np.float32)
    similarities[:, (0, 4)] = 0.0
    similarities[0, 1] = 10.0
    similarities[2, 2] = 10.0
    similarities[4, 3] = 10.0
    lengths = np.asarray(similarities.shape)

    penalized = _affine_score(
        similarities,
        lengths,
        1e-4,
        gap_extend=-1.0,
        gap_open=-3.0,
    )
    free = _affine_score(
        similarities,
        lengths,
        1e-4,
        gap_extend=-1.0,
        gap_open=-3.0,
        free_gap_boundary=(1, 2),
    )

    assert float(penalized) == pytest.approx(24.0)
    assert float(free) == pytest.approx(30.0)


def test_softalign_gap_penalties_reach_affine_alignment(monkeypatch):
    captured = {}

    def fake_affine(
        similarity,
        lengths,
        temperature,
        gap_extend,
        gap_open,
        free_gap_boundary,
    ):
        captured["gap_extend"] = gap_extend
        captured["gap_open"] = gap_open
        captured["free_gap_boundary"] = free_gap_boundary
        return np.asarray([0.0]), np.zeros(np.asarray(similarity).shape)

    monkeypatch.setattr("sabr.alignment._AFFINE_ALIGNMENT", fake_affine)
    _align_reference(
        np.zeros((1, constants.EMBED_DIM)),
        np.zeros((1, constants.EMBED_DIM)),
        "softalign",
    )
    assert captured == {
        "gap_extend": pytest.approx(0.1942468136548996),
        "gap_open": pytest.approx(-2.5441808700561523),
        "free_gap_boundary": -1,
    }


def test_auto_reference_ties_resolve_in_h_k_l_order(monkeypatch):
    references = {
        chain_type: (np.zeros((1, 64)), [1])
        for chain_type in constants.CHAIN_TYPES
    }
    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: references,
    )
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, mode: (
            np.ones((len(query), 1)),
            np.zeros((len(query), 3)),
            1.0,
        ),
    )
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )
    _, selected, _ = align(np.zeros((1, 64)), frozenset(), "auto", 0.0)
    assert selected == "H"


@pytest.mark.parametrize("chain_type", constants.CHAIN_TYPES)
def test_explicit_chain_type_aligns_only_its_reference(monkeypatch, chain_type):
    references = {
        candidate: (np.full((1, 64), index), [1])
        for index, candidate in enumerate(constants.CHAIN_TYPES)
    }
    seen = []

    def fake_align(query, reference, mode):
        seen.append(int(reference[0, 0]))
        return np.ones((len(query), 1)), np.zeros((len(query), 3)), 1.0

    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: references,
    )
    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )
    _, selected, _ = align(np.zeros((1, 64)), frozenset(), chain_type, 0.0)
    assert selected == chain_type
    assert seen == [constants.CHAIN_TYPES.index(chain_type)]


@pytest.mark.parametrize("chain_type", constants.TCR_CHAIN_TYPES)
def test_tcr_chain_type_aligns_only_kappa_reference(monkeypatch, chain_type):
    references = {
        candidate: (np.full((1, 64), index), [1])
        for index, candidate in enumerate(constants.CHAIN_TYPES)
    }
    seen = []

    def fake_align(query, reference, mode):
        seen.append(int(reference[0, 0]))
        return np.ones((len(query), 1)), np.zeros((len(query), 3)), 1.0

    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: references,
    )
    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )

    _, selected, _ = align(np.zeros((1, 64)), frozenset(), chain_type, 0.0)

    assert selected == "K"
    assert seen == [constants.CHAIN_TYPES.index("K")]


def test_scfv_mode_selects_and_corrects_a_composite_reference(monkeypatch):
    representations = (
        *constants.CHAIN_TYPES,
        *constants.SCFV_CHAIN_TYPES,
    )
    references = {}
    for index, representation in enumerate(representations):
        domain_count = len(representation)
        references[representation] = (
            np.full((domain_count, 64), index),
            [1, 129] if domain_count == 2 else [1],
        )

    seen_boundaries = {}

    def fake_align(query, reference, mode, free_gap_boundary=-1):
        reduced = np.zeros((len(query), len(reference)))
        np.fill_diagonal(reduced, 1)
        marker = int(reference[0, 0])
        representation = representations[marker]
        seen_boundaries[representation] = free_gap_boundary
        if representation == "H":
            score = 100.0
        else:
            score = 10.0 if representation == "HK" else 1.0
        return reduced, np.zeros((len(query), len(reference) + 2)), score

    corrected_shapes = []
    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: references if scfv else None,
    )
    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: (
            corrected_shapes.append(alignment.shape) or alignment
        ),
    )

    alignment, selected, score = align(
        np.zeros((2, 64)), frozenset(), "auto", 0.0, scfv=True
    )
    assert selected == "HK"
    assert score == 10.0
    assert alignment.shape == (2, 256)
    assert corrected_shapes == [(2, 128), (2, 128)]
    assert seen_boundaries == {
        "HK": (1,),
        "HL": (1,),
        "KH": (1,),
        "LH": (1,),
    }


def test_requested_candidates_and_higher_order_boundaries_are_exact(
    monkeypatch,
):
    references = {
        "H": (np.zeros((1, 64)), (1,)),
        "K": (np.ones((1, 64)), (1,)),
        "L": (np.full((1, 64), 2), (1,)),
    }
    seen = []

    def fake_align(query, reference, mode, free_gap_boundary=-1):
        seen.append((int(reference[0, 0]), free_gap_boundary))
        reduced = np.zeros((len(query), len(reference)))
        np.fill_diagonal(reduced, 1)
        return reduced, np.zeros((len(query), len(reference) + 2)), 1.0

    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: references,
    )
    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )

    alignment, selected, _ = align(
        np.zeros((3, 64)), frozenset(), "HHK,HHL", 0.0
    )
    assert selected == "HHK"
    assert alignment.shape == (3, 3 * constants.IMGT_MAX_POSITION)
    assert seen == [(0, (1, 2)), (0, (1, 2))]


def test_terminal_gap_penalty_uses_normal_affine_gap_costs():
    alignment = np.zeros((8, 10))
    alignment[2, 3] = 1
    alignment[5, 7] = 1
    expected = 4 * constants.SW_GAP_OPEN + 5 * constants.SW_GAP_EXTEND
    assert _affine_gap_penalty(0) == 0.0
    assert _affine_gap_penalty(1) == constants.SW_GAP_OPEN
    assert _terminal_gap_penalty(alignment) == pytest.approx(expected)
    soft_extend, soft_open = load_gap_penalties("softalign")
    assert _terminal_gap_penalty(alignment, "softalign") == pytest.approx(
        4 * soft_open + 5 * soft_extend
    )


def test_multidomain_terminal_penalty_is_used_only_for_candidate_selection(
    monkeypatch,
):
    references = {
        "H": (np.zeros((1, 64)), [1]),
        "HK": (np.ones((4, 64)), [1, 2, 129, 130]),
    }

    def fake_align(query, reference, mode, free_gap_boundary=-1):
        if len(reference) == 1:
            return np.ones((1, 1)), np.zeros((1, 3)), 8.0
        return (
            np.asarray([[1.0, 0.0, 0.0, 0.0]]),
            np.zeros((1, 6)),
            10.0,
        )

    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: references,
    )
    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )

    alignment, selected, score = align(
        np.zeros((1, 64)), frozenset(), "H,HK", 0.0
    )
    assert selected == "H"
    assert score == 8.0
    assert alignment.shape == (1, constants.IMGT_MAX_POSITION)


def test_huge_internal_run_is_valid_inside_one_cdr():
    alignment = np.zeros((102, 128), dtype=int)
    alignment[0, 26] = 1
    alignment[101, 37] = 1
    _alignment_path(alignment)


def test_de_loop_insertions_form_a_valid_alignment():
    alignment = np.zeros((8, 128), dtype=int)
    alignment[0, 78] = 1
    alignment[1, 79] = 1
    alignment[2, 80] = 1
    alignment[3, 81] = 1
    alignment[5, 82] = 1
    alignment[6, 83] = 1
    alignment[7, 84] = 1
    _alignment_path(alignment)


@pytest.mark.parametrize(
    ("alignment", "message"),
    [
        (np.asarray([[0.0, np.nan]]), "finite"),
        (np.asarray([[0, 2]]), "zeroes and ones"),
        (np.asarray([[1, 1], [0, 0]]), "query row"),
        (np.asarray([[1, 0], [1, 0]]), "IMGT column"),
        (np.asarray([[0, 1], [1, 0]]), "monotonic"),
    ],
)
def test_malformed_alignments_are_rejected(alignment, message):
    with pytest.raises(ValueError, match=message):
        _alignment_path(alignment)


def test_internal_framework_run_is_allowed():
    alignment = np.zeros((4, 128), dtype=int)
    alignment[0, 19] = 1
    alignment[3, 20] = 1
    _alignment_path(alignment)
    states, _, _ = alignment_to_states(alignment)
    assert [state for state, _ in states] == [
        (20, "m"),
        (20, "i"),
        (20, "i"),
        (21, "m"),
    ]


def test_leading_and_trailing_alignment_boundaries_are_allowed():
    for first_row in (2, 3, 4):
        leading = np.zeros((first_row + 1, 128), dtype=int)
        leading[first_row, 2] = 1
        _alignment_path(leading)

    for last_position in (123, 124, 125):
        trailing = np.zeros((3, 128), dtype=int)
        trailing[0, last_position - 1] = 1
        _alignment_path(trailing)


def test_non_finite_reference_score_is_rejected(monkeypatch):
    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: {"H": (np.zeros((1, 64)), (1,))},
    )
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, mode: (
            np.ones((len(query), 1)),
            np.zeros((len(query), 1)),
            float("nan"),
        ),
    )
    with pytest.raises(ValueError, match="non-finite"):
        align(np.zeros((1, 64)), frozenset(), "H", 0.0)


def test_non_finite_similarity_matrix_is_rejected(monkeypatch):
    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise, mode, scfv=False: {"H": (np.zeros((1, 64)), (1,))},
    )
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, mode: (
            np.ones((len(query), 1)),
            np.full((len(query), 1), np.inf),
            1.0,
        ),
    )
    with pytest.raises(ValueError, match="non-finite"):
        align(np.zeros((1, 64)), frozenset(), "H", 0.0)


@pytest.mark.parametrize("candidate_order", ("K,H", "h,k,H"))
def test_public_normalization_and_alignment_ties_use_canonical_order(
    monkeypatch, candidate_order
):
    normalized = _normalize_chain_type(candidate_order, 0.0, "sabr")
    assert normalized == "H,K"
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, mode: (
            np.eye(len(query), len(reference)),
            np.zeros((len(query), len(reference))),
            1.0,
        ),
    )
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )
    _, selected, _ = align(np.zeros((2, 64)), frozenset(), normalized, 0.0)
    assert selected == "H"


def test_winning_multidomain_returns_raw_score_after_terminal_penalty(
    monkeypatch,
):
    references = {
        "H": (np.zeros((1, 64)), (1,)),
        "HK": (np.ones((4, 64)), (1, 2, 129, 130)),
    }
    monkeypatch.setattr(
        "sabr.alignment.load_references", lambda *args, **kwargs: references
    )

    def fake_align(query, reference, mode, free_gap_boundary=-1):
        if len(reference) == 1:
            return np.array([[1.0], [0.0]]), np.zeros((2, 1)), 1.0
        return (
            np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
            np.zeros((2, 4)),
            20.0,
        )

    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )
    alignment, selected, score = align(
        np.zeros((2, 64)), frozenset(), "H,HK", 0.0
    )
    assert selected == "HK"
    assert score == 20.0
    assert np.argwhere(alignment).tolist() == [[0, 1], [1, 128]]


@pytest.mark.parametrize(
    "representation,shape,assignments,message",
    [
        ("HK", (2, 128), [(0, 0)], "256 columns"),
        ("HK", (2, 256), [(0, 0)], "K domain 2"),
        ("HHK", (2, 384), [(0, 0), (1, 256)], "H domain 2"),
        ("HK", (2, 256), [(0, 128), (1, 0)], "monotonic"),
    ],
)
def test_multidomain_rejects_wrong_shape_missing_domains_and_reversed_order(
    representation, shape, assignments, message
):
    alignment = np.zeros(shape, dtype=int)
    for row, column in assignments:
        alignment[row, column] = 1
    with pytest.raises(ValueError, match=message):
        _validate_multidomain_alignment(alignment, representation)


@pytest.mark.parametrize("alignment", (np.zeros((2, 128)), np.ones(2)))
def test_empty_or_nonmatrix_alignment_is_rejected(alignment):
    with pytest.raises(ValueError, match="no assigned|two-dimensional"):
        _alignment_path(alignment)
