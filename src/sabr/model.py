"""NumPy inference for the saved ProteinMPNN/SoftAlign encoders.

SciPy supplies the exact erf-based GELU. Training, stochastic augmentation,
JAX tracing, and Haiku module construction are unnecessary for inference.
"""

import functools
from importlib.resources import files

import numpy as np
from scipy.special import erfc

from sabr import constants


def _linear(x, parameters):
    result = x @ parameters["w"]
    if "b" in parameters:
        result += parameters["b"]
    return result


def _norm(x, parameters):
    mean = np.mean(x, axis=-1, keepdims=True)
    variance = np.var(x, axis=-1, keepdims=True)
    inv = parameters["scale"] * np.reciprocal(np.sqrt(variance + 1e-5))
    return inv * (x - mean) + parameters["offset"]


def _gelu(x):
    return (0.5 * x) * erfc(-x * np.float32(2**-0.5))


def _rbf(distances):
    centers = np.linspace(2, 22, 16, dtype=np.float32)
    return np.exp(-(((distances[..., None] - centers) / 1.25) ** 2))


def _features(coords, parameters):
    n, ca, c, oxygen = np.moveaxis(coords, 1, 0)
    b = ca - n
    d = c - ca
    cb = -0.58273431 * np.cross(b, d) + 0.56802827 * b - 0.54067466 * d + ca
    distances = np.sqrt(
        np.sum((ca[:, None] - ca[None, :]) ** 2, axis=-1) + 1e-6
    )
    # Stable sorting preserves lower-index precedence when distances tie.
    neighbors = np.argsort(distances, axis=-1, kind="stable")[:, :64]
    radial = [_rbf(np.take_along_axis(distances, neighbors, axis=1))]
    pairs = (
        (n, n),
        (c, c),
        (oxygen, oxygen),
        (cb, cb),
        (ca, n),
        (ca, c),
        (ca, oxygen),
        (ca, cb),
        (n, c),
        (n, oxygen),
        (n, cb),
        (cb, c),
        (cb, oxygen),
        (oxygen, c),
        (n, ca),
        (c, ca),
        (oxygen, ca),
        (cb, ca),
        (c, n),
        (oxygen, n),
        (cb, n),
        (c, cb),
        (oxygen, cb),
        (c, oxygen),
    )
    # Compute only selected neighbor distances instead of 24 dense N x N
    # distance matrices. Atom-pair order is part of the trained model.
    for first, second in pairs:
        distance = np.sqrt(
            np.sum((first[:, None] - second[neighbors]) ** 2, axis=-1) + 1e-6
        )
        radial.append(_rbf(distance))
    # Preserve the historical argument inversion: residue indices were
    # chain labels and all actual positional offsets were zero.
    positions = np.where(neighbors == np.arange(len(coords))[:, None], 32, 65)
    positional = parameters[
        "protein_features/~/positional_encodings/~/embedding_linear"
    ]
    positional = positional["w"][positions] + positional["b"]
    edges = np.concatenate((positional, *radial), axis=-1)
    edges = _linear(edges, parameters["protein_features/~/edge_embedding"])
    return _norm(edges, parameters["protein_features/~/norm_edges"]), neighbors


def _messages(nodes, edges, neighbors, parameters, prefix, names):
    center = np.broadcast_to(nodes[:, None], edges.shape)
    joined = np.concatenate((center, edges, nodes[neighbors]), axis=-1)
    first, second, third = (parameters[prefix + name] for name in names)
    return _linear(_gelu(_linear(_gelu(_linear(joined, first)), second)), third)


def _unflatten_parameters(flat_parameters: dict) -> dict:
    parameters = {}
    for key, value in flat_parameters.items():
        current = parameters
        parts = key.split(".")
        for part in parts[:-1]:
            current = current.setdefault(part, {})
        value = np.array(value)
        value.flags.writeable = False
        current[parts[-1]] = value
    return parameters


@functools.cache
def load_parameters(mode: str = "sabr") -> dict:
    """Load the immutable encoder parameters for one alignment mode."""
    filenames = {
        "sabr": "mpnn_encoder.npz",
        "softalign": "softalign_encoder.npz",
    }
    try:
        filename = filenames[mode]
    except KeyError as error:
        raise ValueError(f"mode must be one of {constants.MODES}.") from error
    path = files("sabr.assets") / filename
    with path.open("rb") as handle:
        flat_parameters = dict(np.load(handle, allow_pickle=False))
    return _unflatten_parameters(flat_parameters)


def encode(coords: np.ndarray, mode: str = "sabr") -> np.ndarray:
    """Return 64-dimensional residue embeddings using CPU inference."""
    parameters = load_parameters(mode)
    coords = np.asarray(coords, dtype=np.float32)
    edges, neighbors = _features(coords, parameters)
    nodes = np.zeros((len(coords), constants.EMBED_DIM), dtype=np.float32)
    edges = _linear(edges, parameters["W_e"])
    for layer in range(constants.N_MPNN_LAYERS):
        module = "enc_layer" + (f"_{layer}" if layer else "") + "/~/"
        prefix = module + f"enc{layer}_"
        messages = _messages(
            nodes, edges, neighbors, parameters, prefix, ("W1", "W2", "W3")
        )
        nodes = _norm(
            nodes + np.sum(messages, axis=-2) / 30, parameters[prefix + "norm1"]
        )
        dense = module + f"position_wise_feed_forward/~/enc{layer}_dense_"
        hidden = _gelu(_linear(nodes, parameters[dense + "W_in"]))
        nodes = _norm(
            nodes + _linear(hidden, parameters[dense + "W_out"]),
            parameters[prefix + "norm2"],
        )
        # The last edge update cannot affect the returned node embeddings.
        if layer + 1 < constants.N_MPNN_LAYERS:
            messages = _messages(
                nodes,
                edges,
                neighbors,
                parameters,
                prefix,
                ("W11", "W12", "W13"),
            )
            edges = _norm(edges + messages, parameters[prefix + "norm3"])
    return nodes
