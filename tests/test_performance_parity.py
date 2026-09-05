"""Numerical guards for changes to inference execution."""

from pathlib import Path

import haiku as hk
import jax
import numpy as np
import pytest
from Bio.PDB import PDBParser

from sabr.alignment import _align_reference, _align_references, load_references
from sabr.model import _APPLY_ENCODER, ProteinFeatures, encode, load_parameters
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


@pytest.mark.parametrize("mode", ("sabr", "softalign"))
@pytest.mark.parametrize("length", (63, 64, 65, 95, 96, 97, 122, 250))
def test_encoder_buckets_preserve_unpadded_embeddings(length, mode):
    coords = (
        np.random.default_rng(length)
        .normal(size=(length, 4, 3))
        .astype(np.float32)
    )
    expected = _APPLY_ENCODER(
        load_parameters(mode),
        jax.random.PRNGKey(0),
        coords[None],
        np.ones((1, length)),
        np.ones((1, length)),
        np.arange(length)[None],
    )[0]
    actual = encode(coords, mode)
    assert actual.shape == (length, 64)
    assert isinstance(actual, np.ndarray)
    np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=2e-6)


def test_padding_never_enters_valid_neighbor_sets():
    """Even tied coordinates cannot select a padded neighbor."""
    coordinates = np.zeros((1, 96, 3), dtype=np.float32)
    mask = np.zeros((1, 96), dtype=np.float32)
    mask[:, :65] = 1

    def neighbors(x, valid):
        return ProteinFeatures(64, top_k=64)._dist(x, valid)

    transformed = hk.without_apply_rng(hk.transform(neighbors))
    params = transformed.init(jax.random.PRNGKey(0), coordinates, mask)
    distances, indices = jax.jit(transformed.apply)(params, coordinates, mask)
    assert np.all(np.asarray(indices)[:, :65] < 65)
    assert np.isfinite(np.asarray(distances)[:, :65]).all()
