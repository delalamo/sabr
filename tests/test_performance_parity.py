"""Numerical guards for changes to inference execution."""

from pathlib import Path

import haiku as hk
import jax
import numpy as np
import pytest
from Bio.PDB import PDBParser

from sabr.alignment import _align_reference, _align_references, load_references
from sabr.model import ProteinFeatures, encode
from sabr.structure import extract_chain


@pytest.mark.parametrize("length", (5, 65, 250))
def test_neighbor_rbf_matches_dense_distances(length):
    """Preserve directed atom-pair distances and neighbor ordering."""
    rng = np.random.default_rng(42)
    first = rng.normal(size=(2, length, 3)).astype(np.float32)
    second = rng.normal(size=(2, length, 3)).astype(np.float32)
    indices = rng.integers(0, length, size=(2, length, min(length, 64)))

    def features(a, b, neighbors):
        return ProteinFeatures(64)._get_rbf(a, b, neighbors)

    transformed = hk.without_apply_rng(hk.transform(features))
    params = transformed.init(jax.random.PRNGKey(0), first, second, indices)
    actual = jax.jit(transformed.apply)(params, first, second, indices)
    dense = np.sqrt(
        np.sum((first[:, :, None, :] - second[:, None, :, :]) ** 2, axis=-1)
        + 1e-6
    )
    distances = np.take_along_axis(dense, indices, axis=2)
    centers = np.linspace(2.0, 22.0, 16, dtype=np.float32)
    expected = np.exp(-(((distances[..., None] - centers) / 1.25) ** 2))
    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-6)


@pytest.mark.parametrize("mode", ("sabr", "softalign"))
@pytest.mark.parametrize(
    "filename,chain",
    (
        ("test_heavy_chain.pdb", "F"),
        ("12e8_L.pdb", "L"),
        ("1bjm_A.pdb", "A"),
        ("8sve_L.pdb", "M"),
    ),
)
def test_batched_references_match_scalar_alignments(filename, chain, mode):
    path = Path(__file__).parent / "data" / filename
    structure = PDBParser(QUIET=True).get_structure("case", path)
    data = extract_chain(structure, chain, None)
    query = encode(data.coords, mode)
    references = [entry[0] for entry in load_references(0.0, mode).values()]
    batched = _align_references(query, references, mode)
    for reference, actual in zip(references, batched, strict=True):
        expected = _align_reference(query, reference, mode)
        np.testing.assert_allclose(actual[0], expected[0], rtol=1e-5, atol=1e-6)
        np.testing.assert_array_equal(
            np.round(actual[0]), np.round(expected[0])
        )
        np.testing.assert_allclose(actual[1], expected[1], rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(actual[2], expected[2], rtol=1e-5, atol=1e-6)
