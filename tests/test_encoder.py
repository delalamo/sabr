from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

from sabr.alignment import _align_reference, load_references
from sabr.model import encode
from sabr.structure import extract_chain

DATA = Path(__file__).parent / "data"
ASSETS = Path(__file__).parents[1] / "src" / "sabr" / "assets"


def test_encoder_and_affine_alignment_match_captured_main_baseline():
    baseline = np.load(DATA / "math_baseline.npz")
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    data = extract_chain(structure, "F", None)
    embeddings = encode(data.coords)
    np.testing.assert_allclose(
        embeddings, baseline["embeddings"], rtol=1e-5, atol=2e-6
    )
    softalign_embeddings = encode(data.coords, "softalign")
    assert softalign_embeddings.shape == embeddings.shape
    assert not np.allclose(softalign_embeddings, embeddings)

    reference, positions = load_references(0.0)["H"]
    reduced, similarity, score = _align_reference(embeddings, reference)
    np.testing.assert_allclose(
        reduced, baseline["reduced_alignment"], rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(
        similarity, baseline["similarity"], rtol=1e-5, atol=1e-5
    )
    np.testing.assert_allclose(score, baseline["score"], rtol=1e-5, atol=1e-6)
    np.testing.assert_array_equal(
        np.round(reduced), np.round(baseline["reduced_alignment"])
    )
    np.testing.assert_array_equal(positions, baseline["positions"])
