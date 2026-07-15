import hashlib
from pathlib import Path

import gemmi
import numpy as np
import pytest

from sabr import constants
from sabr.alignment import (
    _align_reference,
    _validate_alignment,
    align,
    create_gap_penalties,
    load_references,
)
from sabr.model import encode, load_parameters
from sabr.structure import extract_chain

DATA = Path(__file__).parent / "data"
ASSETS = Path(__file__).parents[1] / "src" / "sabr" / "assets"


def test_scientific_assets_are_unchanged():
    expected = {
        "mpnn_encoder.npz": (
            "ee70ee43b4c9fdcf7528ad1c7ccfcba20"
            "f47d99de7568939db2f034f21b37ca7"
        ),
        "embeddings_noise_0.0.npz": (
            "1a738d26d0decb75f2fe0b83da1c71ac"
            "91433d80cc3d8270e669e99cf00e24bb"
        ),
        "embeddings_noise_0.2.npz": (
            "619286656b05d0921a4bf7a587bfc767"
            "d93f9e193fd2d2ff6787d7a4b7fa9ef3"
        ),
        "embeddings_noise_0.5.npz": (
            "3a4f51604dd0b2d27352ad269ed9ceb62"
            "42d08dc0a837bf5d3846af145a4444c"
        ),
        "embeddings_noise_1.0.npz": (
            "70495a54ae02a4db70ade0be2252c86b"
            "b50b4b02197b9ac33b042493dad88deb"
        ),
        "embeddings_noise_2.0.npz": (
            "b2f8e886e73655d25c68d23488ad17a"
            "2e441af9fe016da50d28b3e05928861e2"
        ),
    }
    for filename, digest in expected.items():
        assert (
            hashlib.sha256((ASSETS / filename).read_bytes()).hexdigest()
            == digest
        )
    assert not (ASSETS / "embeddings.npz").exists()


def test_encoder_and_affine_alignment_match_captured_main_baseline():
    baseline = np.load(DATA / "math_baseline.npz")
    structure = gemmi.read_structure(str(DATA / "test_heavy_chain.pdb"))
    data = extract_chain(structure, "F", None)
    embeddings = encode(data.coords)
    np.testing.assert_allclose(
        embeddings, baseline["embeddings"], rtol=1e-5, atol=1e-6
    )

    reference, positions = load_references(0.0)["H"]
    reduced, similarity, score = _align_reference(
        embeddings, reference, positions
    )
    np.testing.assert_allclose(
        reduced, baseline["reduced_alignment"], rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(
        similarity, baseline["similarity"], rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(score, baseline["score"], rtol=1e-5, atol=1e-6)
    np.testing.assert_array_equal(
        np.round(reduced), np.round(baseline["reduced_alignment"])
    )
    np.testing.assert_array_equal(positions, baseline["positions"])


def test_gap_penalties_are_fixed_and_uniform():
    positions = [0, 10, 27, 38, 39, 56, 65, 66, 105, 117, 118, 129]
    gap_extend, gap_open = create_gap_penalties(7, positions)
    assert gap_extend.dtype == np.float32
    assert gap_open.dtype == np.float32
    assert np.all(gap_extend == constants.SW_GAP_EXTEND)
    assert np.all(gap_open == constants.SW_GAP_OPEN)


def test_every_noise_asset_has_all_chain_references():
    for level in constants.NOISE_LEVELS:
        references = load_references(level)
        assert tuple(references) == constants.CHAIN_TYPES
        for embeddings, positions in references.values():
            assert embeddings.shape[1] == constants.EMBED_DIM
            assert embeddings.shape[0] == len(positions)


def test_model_assets_are_cached_and_immutable():
    load_parameters.cache_clear()
    first_parameters = load_parameters()
    assert load_parameters() is first_parameters

    load_references.cache_clear()
    first_references = load_references(0.0)
    assert load_references(0.0) is first_references
    for embeddings, positions in first_references.values():
        assert not embeddings.flags.writeable
        assert isinstance(positions, tuple)


def test_auto_reference_ties_resolve_in_h_k_l_order(monkeypatch):
    references = {
        chain_type: (np.zeros((1, 64)), [1])
        for chain_type in constants.CHAIN_TYPES
    }
    monkeypatch.setattr(
        "sabr.alignment.load_references", lambda noise: references
    )
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, positions: (
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

    def fake_align(query, reference, positions):
        seen.append(int(reference[0, 0]))
        return np.ones((len(query), 1)), np.zeros((len(query), 3)), 1.0

    monkeypatch.setattr(
        "sabr.alignment.load_references", lambda noise: references
    )
    monkeypatch.setattr("sabr.alignment._align_reference", fake_align)
    monkeypatch.setattr(
        "sabr.alignment.apply_corrections",
        lambda alignment, gap_indices: alignment,
    )
    _, selected, _ = align(np.zeros((1, 64)), frozenset(), chain_type, 0.0)
    assert selected == chain_type
    assert seen == [constants.CHAIN_TYPES.index(chain_type)]


def test_huge_internal_run_is_valid_inside_one_cdr():
    alignment = np.zeros((102, 128), dtype=int)
    alignment[0, 26] = 1
    alignment[101, 37] = 1
    _validate_alignment(alignment, "K")


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
        _validate_alignment(alignment, "H")


def test_internal_framework_run_is_rejected():
    alignment = np.zeros((4, 128), dtype=int)
    alignment[0, 19] = 1
    alignment[3, 20] = 1
    with pytest.raises(ValueError, match="residue_range"):
        _validate_alignment(alignment, "H")


def test_leading_and_trailing_alignment_boundaries():
    valid_leading = np.zeros((4, 128), dtype=int)
    valid_leading[2, 2] = 1
    valid_leading[3, 3] = 1
    _validate_alignment(valid_leading, "H")

    invalid_leading = np.zeros((4, 128), dtype=int)
    invalid_leading[3, 2] = 1
    with pytest.raises(ValueError, match="non-positive"):
        _validate_alignment(invalid_leading, "H")

    valid_trailing = np.zeros((3, 128), dtype=int)
    valid_trailing[0, 124] = 1
    _validate_alignment(valid_trailing, "H")

    invalid_trailing = np.zeros((3, 128), dtype=int)
    invalid_trailing[0, 123] = 1
    with pytest.raises(ValueError, match="trailing"):
        _validate_alignment(invalid_trailing, "H")


def test_non_finite_reference_score_is_rejected(monkeypatch):
    monkeypatch.setattr(
        "sabr.alignment.load_references",
        lambda noise: {"H": (np.zeros((1, 64)), (1,))},
    )
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, positions: (
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
        lambda noise: {"H": (np.zeros((1, 64)), (1,))},
    )
    monkeypatch.setattr(
        "sabr.alignment._align_reference",
        lambda query, reference, positions: (
            np.ones((len(query), 1)),
            np.full((len(query), 1), np.inf),
            1.0,
        ),
    )
    with pytest.raises(ValueError, match="non-finite"):
        align(np.zeros((1, 64)), frozenset(), "H", 0.0)
