import copy
import json
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser

from sabr import constants, renumber_structure
from sabr.alignment import _align_reference, align, load_references
from sabr.model import encode
from sabr.numbering import number_alignment
from sabr.structure import apply_numbering, extract_chain
from tests.helpers import assert_atomic_content_equal, residue_ids

DATA = Path(__file__).parent / "data"
ASSETS = Path(__file__).parents[1] / "src" / "sabr" / "assets"


@pytest.mark.integration
@pytest.mark.parametrize("mode", constants.MODES)
def test_encoder_and_affine_alignment_match_captured_baseline(mode):
    filename = (
        "math_baseline.npz" if mode == "sabr" else "softalign_math_baseline.npz"
    )
    baseline = np.load(DATA / filename)
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    data = extract_chain(structure, "F", None)
    embeddings = encode(data.coords, mode)
    np.testing.assert_allclose(
        embeddings, baseline["embeddings"], rtol=1e-5, atol=2e-6
    )
    reference, positions = load_references(0.0, mode)["H"]
    reduced, similarity, score = _align_reference(embeddings, reference, mode)
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


@pytest.mark.integration
@pytest.mark.parametrize("mode", constants.MODES)
def test_rigid_transform_and_repeat_preserve_embeddings_and_numbering(
    aligned_case, mode
):
    case = aligned_case("H", mode)
    original = PDBParser(QUIET=True).get_structure("original", case["path"])
    transformed = copy.deepcopy(original)
    # A proper 90-degree rotation, exactly representable before translation.
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert np.linalg.det(rotation) == 1
    for atom in transformed.get_atoms():
        atom.coord = atom.coord.astype(np.float64) @ rotation + np.array(
            [12.0, -7.0, 4.0]
        )
    data = extract_chain(transformed, case["chain"], None)
    embeddings = encode(data.coords, mode)
    np.testing.assert_allclose(
        embeddings, case["embeddings"], rtol=2e-4, atol=2e-5
    )
    np.testing.assert_array_equal(
        encode(case["data"].coords, mode), case["embeddings"]
    )
    baseline = renumber_structure(original, case["chain"], mode=mode)
    numbered = renumber_structure(transformed, case["chain"], mode=mode)
    assert residue_ids(numbered, case["chain"]) == residue_ids(
        baseline, case["chain"]
    )
    repeated = renumber_structure(baseline, case["chain"], mode=mode)
    assert residue_ids(repeated, case["chain"]) == residue_ids(
        baseline, case["chain"]
    )
    assert_atomic_content_equal(original, baseline)
    assert_atomic_content_equal(baseline, repeated)


@pytest.mark.integration
@pytest.mark.parametrize("kind", constants.CHAIN_TYPES)
@pytest.mark.parametrize("noise", constants.NOISE_LEVELS)
def test_every_noise_reference_runs_real_numbering(aligned_case, kind, noise):
    case = aligned_case(kind)
    alignment, selected, score = align(
        case["embeddings"], case["data"].gap_indices, "auto", noise
    )
    assert selected == kind
    assert np.isfinite(score)
    records = number_alignment(
        alignment, case["data"].sequence, "imgt", selected
    )
    golden = json.loads((DATA / "numbering_baseline.json").read_text())[kind][
        "schemes"
    ]["imgt"]
    first = golden["first_row"]
    last = first + len(golden["numbered"])
    assert [row for row, _, _, _ in records] == list(range(first, last))
    assert [aa for _, _, _, aa in records] == list(
        case["data"].sequence[first:last]
    )
    assert len(
        {(number, code.strip()) for _, number, code, _ in records}
    ) == len(records)
    structure = PDBParser(QUIET=True).get_structure("input", case["path"])
    result = apply_numbering(structure, case["chain"], case["data"], records)
    assert_atomic_content_equal(structure, result)
    assert len(residue_ids(result, case["chain"])) == len(case["data"].sequence)
    expected_tail = [
        (records[-1][1] + j, "")
        for j in range(1, len(case["data"].sequence) - last + 1)
    ]
    assert residue_ids(result, case["chain"])[last:] == expected_tail


@pytest.mark.integration
@pytest.mark.parametrize("noise", constants.NOISE_LEVELS)
def test_softalign_noise_does_not_change_alignment_or_numbering(
    aligned_case, noise
):
    case = aligned_case("H", "softalign")
    alignment, selected, score = align(
        case["embeddings"],
        case["data"].gap_indices,
        "auto",
        noise,
        mode="softalign",
    )
    np.testing.assert_array_equal(alignment, case["alignment"])
    assert (selected, score) == (case["selected"], case["score"])
    structure = PDBParser(QUIET=True).get_structure("input", case["path"])
    result = renumber_structure(
        structure, case["chain"], mode="softalign", noise_level=noise
    )
    golden = json.loads(
        (DATA / "softalign_numbering_baseline.json").read_text()
    )["H"]["schemes"]["imgt"]
    assert residue_ids(result, case["chain"]) == [
        (n, code.strip()) for n, code, _ in golden["numbered"]
    ]
