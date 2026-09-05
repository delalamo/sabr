"""Numerical guards for changes to inference execution."""

import haiku as hk
import jax
import numpy as np
import pytest

from sabr.model import ProteinFeatures


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
    expected = np.exp(-((distances[..., None] - centers) / 1.25) ** 2)
    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-6)
