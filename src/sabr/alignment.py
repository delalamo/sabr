#!/usr/bin/env python3
"""Reference alignment using the validated affine Smith-Waterman method."""

import functools
import logging
from importlib.resources import files

import jax
import jax.numpy as jnp
import numpy as np

from sabr import constants
from sabr.corrections import apply_corrections

LOGGER = logging.getLogger(__name__)


def _rotate_for_dp(x, NINF):
    """Rotate a matrix for striped dynamic programming."""
    a, b = x.shape
    ar = jnp.arange(a)[::-1, None]
    br = jnp.arange(b)[None, :]
    i, j = (br - ar) + (a - 1), (ar + br) // 2
    n, m = (a + b - 1), (a + b) // 2
    output = {
        "x": jnp.full([n, m], NINF).at[i, j].set(x),
        "o": (jnp.arange(n) + a % 2) % 2,
    }
    prev = (jnp.full((m, 3), NINF), jnp.full((m, 3), NINF))
    return output, prev, (i, j)


def _soft_maximum(x, temp, NINF, axis=None, mask=None):
    """Compute soft maximum using log-sum-exp."""

    def _logsumexp(y):
        y = jnp.maximum(y, NINF)
        if mask is None:
            return jax.nn.logsumexp(y, axis=axis)
        return y.max(axis) + jnp.log(
            jnp.sum(
                mask * jnp.exp(y - y.max(axis, keepdims=True)),
                axis=axis,
            )
        )

    return temp * _logsumexp(x / temp)


def _cond(cond, true_val, false_val):
    """Conditional selection."""
    return cond * true_val + (1 - cond) * false_val


def _pad(x, shape, NINF):
    """Pad array with NINF values."""
    return jnp.pad(x, shape, constant_values=(NINF, NINF))


def _apply_length_mask(x, lengths, NINF):
    """Apply length mask to similarity matrix."""
    a, b = x.shape
    real_a, real_b = lengths
    mask = (jnp.arange(a) < real_a)[:, None] * (jnp.arange(b) < real_b)[None, :]
    return x + NINF * (1 - mask), mask


def _rotate_gap_matrix(x):
    """Rotate a gap penalty matrix for striped DP, using 0 as fill value."""
    a, b = x.shape
    ar = jnp.arange(a)[::-1, None]
    br = jnp.arange(b)[None, :]
    i, j = (br - ar) + (a - 1), (ar + br) // 2
    n, m = (a + b - 1), (a + b) // 2
    return jnp.full([n, m], 0.0).at[i, j].set(x)


def _affine_score(
    similarities,
    lengths,
    temperature,
    gap_extend,
    gap_open,
    NINF=-1e30,
):
    """Return the smooth affine Smith-Waterman score for one matrix."""

    def step(previous, stripe):
        h2, h1 = previous
        aligned = jnp.pad(h2, [[0, 0], [0, 1]]) + stripe["x"][:, None]
        right = _cond(
            stripe["o"],
            _pad(h1[:-1], ([1, 0], [0, 0]), NINF),
            h1,
        )
        down = _cond(
            stripe["o"],
            h1,
            _pad(h1[1:], ([0, 1], [0, 0]), NINF),
        )
        extension = stripe["gap"]
        opening = stripe["open"]
        right += jnp.stack([opening, extension, opening], axis=-1)
        down += jnp.stack([opening, opening, extension], axis=-1)
        right = right[:, :2]

        h0 = jnp.stack(
            [
                _soft_maximum(aligned, temperature, NINF, axis=-1),
                _soft_maximum(right, temperature, NINF, axis=-1),
                _soft_maximum(down, temperature, NINF, axis=-1),
            ],
            axis=-1,
        )
        return (h1, h0), h0

    similarities, mask = _apply_length_mask(similarities, lengths, NINF)
    stripes, previous, indices = _rotate_for_dp(similarities[:-1, :-1], NINF)
    stripes["gap"] = _rotate_gap_matrix(gap_extend[:-1, :-1])
    stripes["open"] = _rotate_gap_matrix(gap_open[:-1, :-1])
    scores = jax.lax.scan(step, previous, stripes, unroll=2)[-1][indices]
    return _soft_maximum(
        scores + similarities[1:, 1:, None],
        temperature,
        NINF,
        mask=mask[1:, 1:, None],
    )


_AFFINE_ALIGNMENT = jax.jit(
    jax.vmap(jax.value_and_grad(_affine_score), (0, 0, None, 0, 0))
)


def create_gap_penalties(query_length: int, positions: list) -> tuple:
    """Create the fixed position-dependent affine gap penalties."""
    shape = (query_length, len(positions))
    gap_extend = np.full(shape, constants.SW_GAP_EXTEND, dtype=np.float32)
    gap_open = np.full(shape, constants.SW_GAP_OPEN, dtype=np.float32)
    cdr_positions = set()
    for start, end in constants.IMGT_LOOPS.values():
        cdr_positions.update(range(start, end + 1))
    for column, position in enumerate(positions):
        if position in cdr_positions:
            gap_open[:, column] = 0.0
    return gap_extend, gap_open


@functools.cache
def load_references(noise_level: float) -> dict:
    """Load the H, K, and L reference embeddings for one noise level."""
    filename = f"embeddings_noise_{noise_level:.1f}.npz"
    path = files("sabr.assets") / filename
    with path.open("rb") as handle:
        data = np.load(handle, allow_pickle=True)["arr_0"].item()
    references = {}
    for chain_type in constants.CHAIN_TYPES:
        embeddings = np.asarray(data[chain_type]["array"])
        embeddings.flags.writeable = False
        references[chain_type] = (
            embeddings,
            tuple(int(position) for position in data[chain_type]["idxs"]),
        )
    return references


def _validate_alignment(alignment: np.ndarray, chain_type: str) -> None:
    """Reject ambiguous paths before they are converted to numbering states."""
    if alignment.ndim != 2 or not np.isfinite(alignment).all():
        raise ValueError("Alignment must be a finite two-dimensional matrix.")
    if not np.isin(alignment, (0, 1)).all():
        raise ValueError("Alignment matrix must contain only zeroes and ones.")
    if np.any(alignment.sum(axis=1) > 1):
        raise ValueError("Alignment assigns a query row more than once.")
    if np.any(alignment.sum(axis=0) > 1):
        raise ValueError("Alignment assigns an IMGT column more than once.")

    path = np.argwhere(alignment == 1)
    if not len(path):
        raise ValueError("Alignment contains no assigned residues.")
    rows = path[:, 0]
    columns = path[:, 1]
    if np.any(np.diff(rows) <= 0) or np.any(np.diff(columns) <= 0):
        raise ValueError(
            "Alignment assignments are not strictly monotonic; use "
            "residue_range to select one antibody domain."
        )

    first_row = int(rows[0])
    first_position = int(columns[0]) + 1
    if first_row and first_position - first_row < 1:
        raise ValueError(
            "N-terminal residues would require non-positive numbering; "
            "use residue_range to select the antibody domain."
        )

    regions = list(constants.IMGT_LOOPS.values())
    if chain_type in ("K", "L"):
        regions.append((79, 84))
    for index, (left_row, right_row) in enumerate(zip(rows, rows[1:])):
        if right_row == left_row + 1:
            continue
        left_position = int(columns[index]) + 1
        right_position = int(columns[index + 1]) + 1
        if not any(
            start <= left_position <= end and start <= right_position <= end
            for start, end in regions
        ):
            raise ValueError(
                f"Unassigned query rows {left_row + 1}-{right_row - 1} "
                f"are bracketed by IMGT {left_position} and "
                f"{right_position}; use residue_range to select one "
                "antibody domain."
            )

    last_row = int(rows[-1])
    last_position = int(columns[-1]) + 1
    if last_row < alignment.shape[0] - 1 and last_position < 125:
        raise ValueError(
            f"Unassigned trailing query rows follow IMGT {last_position}; "
            "use residue_range to select the antibody domain."
        )


def _align_reference(query: np.ndarray, reference: np.ndarray, positions: list):
    anchor = np.zeros((1, reference.shape[1]), dtype=reference.dtype)
    augmented_reference = np.concatenate((anchor, reference, anchor), axis=0)
    augmented_positions = [0, *positions, 129]
    gap_extend, gap_open = create_gap_penalties(
        query.shape[0], augmented_positions
    )
    query_batch = jnp.asarray(query[None, :])
    reference_batch = jnp.asarray(augmented_reference[None, :])
    lengths = jnp.asarray([[query.shape[0], augmented_reference.shape[0]]])
    similarity = jnp.einsum("nia,nja->nij", query_batch, reference_batch)
    scores, soft_alignment = _AFFINE_ALIGNMENT(
        similarity,
        lengths,
        constants.DEFAULT_TEMPERATURE,
        jnp.asarray(gap_extend[None, :]),
        jnp.asarray(gap_open[None, :]),
    )
    return (
        np.asarray(soft_alignment[0])[:, 1:-1],
        np.asarray(similarity[0]),
        float(scores[0]),
    )


def align(
    query: np.ndarray,
    gap_indices: frozenset,
    chain_type: str,
    noise_level: float,
) -> tuple:
    """Align query embeddings and return corrected IMGT alignment metadata."""
    references = load_references(noise_level)
    candidates = (
        constants.CHAIN_TYPES if chain_type == "auto" else (chain_type,)
    )
    best = None
    for candidate in candidates:
        reference, positions = references[candidate]
        reduced, similarity, score = _align_reference(
            query, reference, positions
        )
        if (
            not np.isfinite(score)
            or not np.isfinite(reduced).all()
            or not np.isfinite(similarity).all()
        ):
            raise ValueError(
                f"{candidate} reference produced a non-finite alignment."
            )
        LOGGER.info("%s reference score: %.4f", candidate, score)
        if best is None or score > best[0]:
            best = (score, candidate, reduced, positions)

    score, selected_type, reduced, positions = best
    full_alignment = np.zeros(
        (query.shape[0], constants.IMGT_MAX_POSITION), dtype=reduced.dtype
    )
    full_alignment[:, np.asarray(positions) - 1] = reduced
    full_alignment = np.round(full_alignment).astype(int)
    corrected = apply_corrections(full_alignment, gap_indices=gap_indices)
    _validate_alignment(corrected, selected_type)
    return corrected, selected_type, score
