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
    mode: str,
    scfv: bool,
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
    if isinstance(noise_level, bool):
        raise ValueError(
            f"noise_level must be one of {constants.NOISE_LEVELS}."
        )
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
            or not all(
                isinstance(value, int) and not isinstance(value, bool)
                for value in residue_range
            )
        ):
            raise ValueError(
                "residue_range must be an inclusive (start, end) tuple."
            )
        if residue_range[0] > residue_range[1]:
            raise ValueError("residue_range start must not exceed its end.")
    if not isinstance(mode, str) or mode.lower() not in constants.MODES:
        raise ValueError(f"mode must be one of {', '.join(constants.MODES)}.")
    if not isinstance(scfv, bool):
        raise ValueError("scfv must be a boolean.")
    if scfv and normalized_chain_type != "auto":
        raise ValueError("scfv requires chain_type='auto'.")
    return (
        scheme.lower(),
        normalized_chain_type,
        normalized_noise,
        mode.lower(),
        scfv,
    )


def renumber_structure(
    structure,
    chain: str,
    scheme: str = "imgt",
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
    mode: str = "sabr",
    scfv: bool = False,
):
    """Return a non-mutating, same-type renumbered structure.

    ``mode="softalign"`` selects the SoftAlign encoder, references, and gap
    penalties as one scientifically consistent parameter set. ``scfv=True``
    adds the four supported two-domain composite references.
    """
    scheme, chain_type, noise_level, mode, scfv = _validate_options(
        scheme, chain_type, noise_level, residue_range, mode, scfv
    )
    data = extract_chain(structure, chain, residue_range)
    embeddings = encode(data.coords, mode)
    imgt_alignment, selected_type, score = align(
        embeddings,
        data.gap_indices,
        chain_type,
        noise_level,
        mode=mode,
        scfv=scfv,
    )
    LOGGER.info(
        "Selected %s reference with alignment score %.4f.",
        selected_type,
        score,
    )
    numbered = number_alignment(
        imgt_alignment,
        data.sequence,
        scheme,
        selected_type,
    )
    return apply_numbering(structure, chain, data, numbered)
