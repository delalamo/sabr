"""The public, in-memory SAbR API."""

import logging

from sabr import constants
from sabr.alignment import align
from sabr.model import encode
from sabr.numbering import number_alignment
from sabr.structure import apply_numbering, extract_chain

LOGGER = logging.getLogger(__name__)


def _validate_options(
    scheme: str,
    chain_type: str,
    noise_level: float,
    residue_range: tuple | None,
) -> tuple:
    if (
        not isinstance(scheme, str)
        or scheme.lower() not in constants.NUMBERING_SCHEMES
    ):
        raise ValueError(
            f"scheme must be one of {', '.join(constants.NUMBERING_SCHEMES)}."
        )
    if not isinstance(chain_type, str):
        raise ValueError("chain_type must be 'auto', 'H', 'K', or 'L'.")
    normalized_chain_type = (
        "auto" if chain_type.lower() == "auto" else chain_type.upper()
    )
    if normalized_chain_type not in ("auto", *constants.CHAIN_TYPES):
        raise ValueError("chain_type must be 'auto', 'H', 'K', or 'L'.")
    try:
        normalized_noise = float(noise_level)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"noise_level must be one of {constants.NOISE_LEVELS}."
        ) from error
    if normalized_noise not in constants.NOISE_LEVELS:
        raise ValueError(
            f"noise_level must be one of {constants.NOISE_LEVELS}."
        )
    if residue_range is not None:
        if (
            not isinstance(residue_range, tuple)
            or len(residue_range) != 2
            or not all(isinstance(value, int) for value in residue_range)
        ):
            raise ValueError(
                "residue_range must be an inclusive (start, end) tuple."
            )
        if residue_range[0] > residue_range[1]:
            raise ValueError("residue_range start must not exceed its end.")
    return scheme.lower(), normalized_chain_type, normalized_noise


def renumber_structure(
    structure,
    chain: str,
    scheme: str = "imgt",
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
):
    """Return a non-mutating, same-type renumbered structure."""
    scheme, chain_type, noise_level = _validate_options(
        scheme, chain_type, noise_level, residue_range
    )
    data = extract_chain(structure, chain, residue_range)
    embeddings = encode(data.coords)
    imgt_alignment, selected_type, score = align(
        embeddings,
        data.gap_indices,
        chain_type,
        noise_level,
    )
    LOGGER.info(
        "Selected %s reference with alignment score %.4f.",
        selected_type,
        score,
    )
    numbered, first_row = number_alignment(
        imgt_alignment,
        data.sequence,
        scheme,
        selected_type,
    )
    return apply_numbering(structure, chain, data, numbered, first_row)
