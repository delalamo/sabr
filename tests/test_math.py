import hashlib
from pathlib import Path

import gemmi
import numpy as np
import pytest

from sabr import constants
from sabr.alignment import (
    _align_reference,
    align,
    create_gap_penalties,
    load_references,
)
from sabr.model import encode
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


def test_gap_penalties_are_fixed_and_only_cdr_openings_are_free():
    positions = [0, 10, 27, 38, 39, 56, 65, 66, 105, 117, 118, 129]
    gap_extend, gap_open = create_gap_penalties(7, positions)
    assert gap_extend.dtype == np.float32
    assert gap_open.dtype == np.float32
    assert np.all(gap_extend == constants.SW_GAP_EXTEND)
    for column, position in enumerate(positions):
        expected = (
            0.0
            if position in {27, 38, 56, 65, 105, 117}
            else constants.SW_GAP_OPEN
        )
        assert np.all(gap_open[:, column] == expected)


def test_every_noise_asset_has_all_chain_references():
    for level in constants.NOISE_LEVELS:
        references = load_references(level)
        assert tuple(references) == constants.CHAIN_TYPES
        for embeddings, positions in references.values():
            assert embeddings.shape[1] == constants.EMBED_DIM
            assert embeddings.shape[0] == len(positions)


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
        lambda alignment, chain_type, gap_indices: alignment,
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
        lambda alignment, selected, gap_indices: alignment,
    )
    _, selected, _ = align(np.zeros((1, 64)), frozenset(), chain_type, 0.0)
    assert selected == chain_type
    assert seen == [constants.CHAIN_TYPES.index(chain_type)]
