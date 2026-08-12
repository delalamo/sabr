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


def _rotate_for_dp(x, NINF, free_gap_boundary=-1):
    """Rotate a matrix for striped dynamic programming."""
    a, b = x.shape
    ar = jnp.arange(a)[::-1, None]
    br = jnp.arange(b)[None, :]
    i, j = (br - ar) + (a - 1), (ar + br) // 2
    n, m = (a + b - 1), (a + b) // 2
    free_gap = jnp.broadcast_to(br == free_gap_boundary, (a, b))
    output = {
        "x": jnp.full([n, m], NINF).at[i, j].set(x),
        "o": (jnp.arange(n) + a % 2) % 2,
        "free_gap": jnp.full([n, m], False).at[i, j].set(free_gap),
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


def _affine_score(
    similarities,
    lengths,
    temperature,
    gap_extend,
    gap_open,
    free_gap_boundary=-1,
    NINF=-1e30,
):
    """Return a score, optionally making one reference gap boundary free."""
    right_penalties = jnp.asarray(
        [
            gap_open,
            gap_extend,
            gap_open,
        ],
        dtype=similarities.dtype,
    )
    down_penalties = jnp.asarray(
        [
            gap_open,
            gap_open,
            gap_extend,
        ],
        dtype=similarities.dtype,
    )

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
        right += right_penalties
        down += jnp.where(
            stripe["free_gap"][:, None],
            0,
            down_penalties,
        )
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
    stripes, previous, indices = _rotate_for_dp(
        similarities[:-1, :-1], NINF, free_gap_boundary
    )
    scores = jax.lax.scan(step, previous, stripes, unroll=2)[-1][indices]
    return _soft_maximum(
        scores + similarities[1:, 1:, None],
        temperature,
        NINF,
        mask=mask[1:, 1:, None],
    )


_AFFINE_ALIGNMENT = jax.jit(
    jax.vmap(
        jax.value_and_grad(_affine_score),
        (0, 0, None, None, None, None),
    )
)


@functools.cache
def load_references(
    noise_level: float,
    mode: str = "sabr",
    scfv: bool = False,
) -> dict:
    """Load single-domain and optional scFv references for one mode."""
    if mode == "sabr":
        filename = f"embeddings_noise_{noise_level:.1f}.npz"
    elif mode == "softalign":
        filename = "softalign_embeddings.npz"
    else:
        raise ValueError(f"mode must be one of {constants.MODES}.")
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

    if scfv:
        for representation in constants.SCFV_CHAIN_TYPES:
            first_type, second_type = representation.split(":")
            first_embeddings, first_positions = references[first_type]
            second_embeddings, second_positions = references[second_type]
            embeddings = np.concatenate(
                (first_embeddings, second_embeddings), axis=0
            )
            embeddings.flags.writeable = False
            references[representation] = (
                embeddings,
                (
                    *first_positions,
                    *(
                        position + constants.IMGT_MAX_POSITION
                        for position in second_positions
                    ),
                ),
            )
    return references


@functools.cache
def load_gap_penalties(mode: str = "sabr") -> tuple[float, float]:
    """Return the gap extension and opening values for one mode."""
    if mode == "sabr":
        return constants.SW_GAP_EXTEND, constants.SW_GAP_OPEN
    if mode != "softalign":
        raise ValueError(f"mode must be one of {constants.MODES}.")
    path = files("sabr.assets") / "softalign_gap.npz"
    with path.open("rb") as handle:
        data = np.load(handle, allow_pickle=False)
        return float(data["gap_extend"]), float(data["gap_open"])


def _alignment_path(alignment: np.ndarray) -> np.ndarray:
    """Validate common matrix invariants and return its assigned path."""
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
    return path


def _validate_alignment(alignment: np.ndarray, chain_type: str) -> None:
    """Reject ambiguous paths before they are converted to numbering states."""
    path = _alignment_path(alignment)
    rows = path[:, 0]
    columns = path[:, 1]

    first_row = int(rows[0])
    first_position = int(columns[0]) + 1
    if first_row and first_position - first_row < 1:
        raise ValueError(
            "N-terminal residues would require non-positive numbering; "
            "use residue_range to select the antibody domain."
        )

    regions = [*constants.IMGT_LOOPS.values(), (79, 85)]
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


def _validate_scfv_alignment(
    alignment: np.ndarray, representation: str
) -> None:
    """Validate two ordered IMGT-domain blocks with an optional linker."""
    expected_columns = 2 * constants.IMGT_MAX_POSITION
    if alignment.ndim != 2 or alignment.shape[1] != expected_columns:
        raise ValueError(
            f"scFv alignment must have {expected_columns} columns."
        )
    _alignment_path(alignment)

    domain_rows = []
    for domain_index, chain_type in enumerate(representation.split(":")):
        start = domain_index * constants.IMGT_MAX_POSITION
        end = start + constants.IMGT_MAX_POSITION
        domain = alignment[:, start:end]
        assigned_rows = np.flatnonzero(domain.sum(axis=1))
        if not len(assigned_rows):
            raise ValueError(
                f"scFv {chain_type} domain contains no assigned residues."
            )
        domain_rows.append((int(assigned_rows[0]), int(assigned_rows[-1])))
        if domain_index == 0:
            _validate_alignment(domain[: assigned_rows[-1] + 1], chain_type)
        else:
            _validate_alignment(domain[assigned_rows[0] :], chain_type)

    if domain_rows[0][1] >= domain_rows[1][0]:
        raise ValueError("scFv domain assignments overlap or are out of order.")


def _align_reference(
    query: np.ndarray,
    reference: np.ndarray,
    mode: str = "sabr",
    free_gap_boundary: int = -1,
):
    """Align one reference, optionally allowing a free query insertion."""
    anchor = np.zeros((1, reference.shape[1]), dtype=reference.dtype)
    augmented_reference = np.concatenate((anchor, reference, anchor), axis=0)
    query_batch = jnp.asarray(query[None, :])
    reference_batch = jnp.asarray(augmented_reference[None, :])
    lengths = jnp.asarray([[query.shape[0], augmented_reference.shape[0]]])
    similarity = jnp.einsum("nia,nja->nij", query_batch, reference_batch)
    gap_extend, gap_open = load_gap_penalties(mode)
    scores, soft_alignment = _AFFINE_ALIGNMENT(
        similarity,
        lengths,
        constants.DEFAULT_TEMPERATURE,
        gap_extend,
        gap_open,
        free_gap_boundary,
    )
    return (
        np.asarray(soft_alignment[0])[:, 1:-1],
        np.asarray(similarity[0]),
        float(scores[0]),
    )


def _affine_gap_penalty(
    length: int,
    gap_extend: float = constants.SW_GAP_EXTEND,
    gap_open: float = constants.SW_GAP_OPEN,
) -> float:
    """Return the normal affine score for one gap of the given length."""
    if length <= 0:
        return 0.0
    return gap_open + (length - 1) * gap_extend


def _terminal_gap_penalty(
    alignment: np.ndarray,
    mode: str = "sabr",
) -> float:
    """Score unaligned query and reference termini as affine gaps."""
    discrete = np.round(alignment).astype(int)
    path = np.argwhere(discrete == 1)
    if not len(path):
        raise ValueError("Alignment contains no assigned residues.")

    first_query, first_reference = path[0]
    last_query, last_reference = path[-1]
    terminal_lengths = (
        int(first_query),
        int(first_reference),
        alignment.shape[0] - int(last_query) - 1,
        alignment.shape[1] - int(last_reference) - 1,
    )
    gap_extend, gap_open = load_gap_penalties(mode)
    return sum(
        _affine_gap_penalty(length, gap_extend, gap_open)
        for length in terminal_lengths
    )


def align(
    query: np.ndarray,
    gap_indices: frozenset,
    chain_type: str,
    noise_level: float,
    mode: str = "sabr",
    scfv: bool = False,
) -> tuple:
    """Align query embeddings and return corrected IMGT alignment metadata."""
    references = load_references(noise_level, mode, scfv=scfv)
    candidates = tuple(references) if chain_type == "auto" else (chain_type,)
    best = None
    for candidate in candidates:
        reference, positions = references[candidate]
        if ":" in candidate:
            free_gap_boundary = sum(
                position <= constants.IMGT_MAX_POSITION
                for position in positions
            )
            reduced, similarity, score = _align_reference(
                query,
                reference,
                mode,
                free_gap_boundary=free_gap_boundary,
            )
        else:
            reduced, similarity, score = _align_reference(
                query, reference, mode
            )
        if (
            not np.isfinite(score)
            or not np.isfinite(reduced).all()
            or not np.isfinite(similarity).all()
        ):
            raise ValueError(
                f"{candidate} reference produced a non-finite alignment."
            )
        selection_score = score
        if ":" in candidate:
            terminal_penalty = _terminal_gap_penalty(reduced, mode)
            selection_score += terminal_penalty
            LOGGER.info(
                "%s reference score: %.4f; terminal gap penalty: %.4f; "
                "selection score: %.4f",
                candidate,
                score,
                terminal_penalty,
                selection_score,
            )
        else:
            LOGGER.info("%s reference score: %.4f", candidate, score)
        if best is None or selection_score > best[0]:
            best = (
                selection_score,
                score,
                candidate,
                reduced,
                positions,
            )

    _, score, selected_type, reduced, positions = best
    domain_count = (max(positions) - 1) // constants.IMGT_MAX_POSITION + 1
    full_alignment = np.zeros(
        (query.shape[0], domain_count * constants.IMGT_MAX_POSITION),
        dtype=reduced.dtype,
    )
    full_alignment[:, np.asarray(positions) - 1] = reduced
    full_alignment = np.round(full_alignment).astype(int)
    if ":" in selected_type:
        corrected = full_alignment.copy()
        for domain_index, _ in enumerate(selected_type.split(":")):
            start = domain_index * constants.IMGT_MAX_POSITION
            end = start + constants.IMGT_MAX_POSITION
            corrected[:, start:end] = apply_corrections(
                corrected[:, start:end],
                gap_indices=gap_indices,
            )
        _validate_scfv_alignment(corrected, selected_type)
    else:
        corrected = apply_corrections(full_alignment, gap_indices=gap_indices)
        _validate_alignment(corrected, selected_type)
    return corrected, selected_type, score
