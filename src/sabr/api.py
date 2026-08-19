"""The public, in-memory SAbR API."""

import logging
from dataclasses import dataclass

from sabr import constants
from sabr.alignment import align
from sabr.model import encode
from sabr.numbering import number_alignment
from sabr.structure import apply_numbering, extract_chain

LOGGER = logging.getLogger(__name__)


def _normalize_choice(value: str, name: str, choices: tuple[str, ...]) -> str:
    """Normalize a case-insensitive string constrained to known choices."""
    if not isinstance(value, str) or value.lower() not in choices:
        raise ValueError(f"{name} must be one of {', '.join(choices)}.")
    return value.lower()


def _normalize_chain_type(chain_type: str) -> str:
    if not isinstance(chain_type, str):
        raise ValueError("chain_type must be 'auto', 'H', 'K', or 'L'.")
    normalized = "auto" if chain_type.lower() == "auto" else chain_type.upper()
    if normalized not in ("auto", *constants.CHAIN_TYPES):
        raise ValueError("chain_type must be 'auto', 'H', 'K', or 'L'.")
    return normalized


def _normalize_noise_level(noise_level: float) -> float:
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
    return normalized_noise


def _validate_residue_range(
    residue_range: tuple[int, int] | None,
) -> None:
    if residue_range is None:
        return
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


@dataclass(frozen=True, slots=True)
class _RenumberOptions:
    """Validated and normalized arguments for one renumbering operation."""

    scheme: str
    chain_type: str
    noise_level: float
    residue_range: tuple[int, int] | None
    mode: str
    scfv: bool

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "scheme",
            _normalize_choice(
                self.scheme, "scheme", constants.NUMBERING_SCHEMES
            ),
        )
        object.__setattr__(
            self, "chain_type", _normalize_chain_type(self.chain_type)
        )
        object.__setattr__(
            self, "noise_level", _normalize_noise_level(self.noise_level)
        )
        object.__setattr__(
            self, "mode", _normalize_choice(self.mode, "mode", constants.MODES)
        )
        _validate_residue_range(self.residue_range)
        if not isinstance(self.scfv, bool):
            raise ValueError("scfv must be a boolean.")
        if self.scfv and self.chain_type != "auto":
            raise ValueError("scfv requires chain_type='auto'.")


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
    options = _RenumberOptions(
        scheme=scheme,
        chain_type=chain_type,
        noise_level=noise_level,
        residue_range=residue_range,
        mode=mode,
        scfv=scfv,
    )
    data = extract_chain(structure, chain, options.residue_range)
    embeddings = encode(data.coords, options.mode)
    imgt_alignment, selected_type, score = align(
        embeddings,
        data.gap_indices,
        options.chain_type,
        options.noise_level,
        mode=options.mode,
        scfv=options.scfv,
    )
    LOGGER.info(
        "Selected %s reference with alignment score %.4f.",
        selected_type,
        score,
    )
    numbered = number_alignment(
        imgt_alignment,
        data.sequence,
        options.scheme,
        selected_type,
    )
    return apply_numbering(structure, chain, data, numbered)
