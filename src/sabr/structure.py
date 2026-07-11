"""Type-preserving BioPython and Gemmi structure handling."""

import copy

import gemmi
import numpy as np
from Bio.PDB.Structure import Structure as BioStructure

from sabr import constants


class _ChainData:
    """Private, normalized view of the selected polymer residues."""

    def __init__(
        self,
        coords: np.ndarray,
        sequence: str,
        residue_ids: list,
        residue_indices: list,
        gap_indices: frozenset,
    ):
        self.coords = coords
        self.sequence = sequence
        self.residue_ids = residue_ids
        self.residue_indices = residue_indices
        self.gap_indices = gap_indices


def _compute_cb(n_coord: np.ndarray, ca_coord: np.ndarray, c_coord: np.ndarray):
    epsilon = 1e-8

    def normalize(vector):
        length = np.sqrt(
            np.square(vector).sum(axis=-1, keepdims=True) + epsilon
        )
        return vector / length

    direction = normalize(n_coord - ca_coord)
    normal = normalize(np.cross(n_coord - c_coord, direction))
    return ca_coord + (
        constants.CB_BOND_LENGTH * np.cos(constants.CB_BOND_ANGLE) * direction
        + constants.CB_BOND_LENGTH
        * np.sin(constants.CB_BOND_ANGLE)
        * np.cos(constants.CB_DIHEDRAL)
        * np.cross(normal, direction)
        + constants.CB_BOND_LENGTH
        * np.sin(constants.CB_BOND_ANGLE)
        * np.sin(constants.CB_DIHEDRAL)
        * (-normal)
    )


def _validate_structure(structure, chain: str):
    if not isinstance(structure, (BioStructure, gemmi.Structure)):
        raise TypeError(
            "structure must be a Bio.PDB.Structure.Structure or "
            "gemmi.Structure object."
        )
    if not isinstance(chain, str) or not chain:
        raise ValueError("chain must be a non-empty string.")
    if len(structure) != 1:
        raise ValueError(
            "SAbR requires a structure containing exactly one model."
        )


def _find_chain(structure, chain: str):
    _validate_structure(structure, chain)
    model = next(iter(structure))
    for candidate in model:
        candidate_name = (
            candidate.id
            if isinstance(structure, BioStructure)
            else candidate.name
        )
        if candidate_name == chain:
            return candidate
    names = [
        candidate.id if isinstance(structure, BioStructure) else candidate.name
        for candidate in model
    ]
    raise ValueError(f"Chain '{chain}' not found. Available chains: {names}.")


def _bio_residue_data(residue) -> tuple:
    if residue.id[0].strip():
        return None
    missing = [atom for atom in ("N", "CA", "C") if atom not in residue]
    if missing:
        label = f"{residue.resname} {residue.id[1]}{residue.id[2].strip()}"
        raise ValueError(
            f"Residue {label} is missing required backbone atom(s): "
            f"{', '.join(missing)}."
        )
    coords = tuple(
        np.asarray(residue[atom].coord, dtype=np.float64)
        for atom in ("N", "CA", "C")
    )
    return residue.id[1], residue.id[2].strip(), residue.resname, coords


def _gemmi_residue_data(residue) -> tuple:
    if residue.het_flag != "A":
        return None
    atoms = [residue.find_atom(name, "*") for name in ("N", "CA", "C")]
    missing = [name for name, atom in zip(("N", "CA", "C"), atoms) if not atom]
    if missing:
        label = f"{residue.name} {residue.seqid}"
        raise ValueError(
            f"Residue {label} is missing required backbone atom(s): "
            f"{', '.join(missing)}."
        )
    coords = tuple(
        np.asarray((atom.pos.x, atom.pos.y, atom.pos.z), dtype=np.float64)
        for atom in atoms
    )
    return residue.seqid.num, residue.seqid.icode.strip(), residue.name, coords


def _detect_gaps(coords: np.ndarray) -> frozenset:
    if len(coords) < 2:
        return frozenset()
    distances = np.linalg.norm(coords[:-1, 2] - coords[1:, 0], axis=1)
    return frozenset(
        int(index)
        for index in np.where(distances > constants.PEPTIDE_BOND_MAX_DISTANCE)[
            0
        ]
    )


def extract_chain(
    structure, chain: str, residue_range: tuple | None
) -> _ChainData:
    """Normalize the selected polymer residues without modifying the input."""
    target = _find_chain(structure, chain)
    is_biopython = isinstance(structure, BioStructure)
    coords = []
    sequence = []
    residue_ids = []
    residue_indices = []

    for index, residue in enumerate(target):
        if is_biopython:
            is_polymer = not residue.id[0].strip()
            number = residue.id[1]
        else:
            is_polymer = residue.het_flag == "A"
            number = residue.seqid.num
        if not is_polymer:
            continue
        if residue_range is not None and not (
            residue_range[0] <= number <= residue_range[1]
        ):
            continue

        data = (
            _bio_residue_data(residue)
            if is_biopython
            else _gemmi_residue_data(residue)
        )
        number, insertion_code, residue_name, backbone = data
        n_coord, ca_coord, c_coord = backbone
        residue_coords = np.stack(
            (
                n_coord,
                ca_coord,
                c_coord,
                _compute_cb(n_coord, ca_coord, c_coord),
            )
        )
        if not np.isfinite(residue_coords).all():
            raise ValueError(
                f"Residue {number}{insertion_code} has non-finite coordinates."
            )
        coords.append(residue_coords)
        sequence.append(constants.AA_3TO1.get(residue_name, "X"))
        residue_ids.append((number, insertion_code))
        residue_indices.append(index)

    if not coords:
        description = "the requested range" if residue_range else "the chain"
        raise ValueError(
            "No polymer residues with complete backbone found in "
            f"{description}."
        )
    stacked_coords = np.stack(coords)
    return _ChainData(
        stacked_coords,
        "".join(sequence),
        residue_ids,
        residue_indices,
        _detect_gaps(stacked_coords),
    )


def _new_residue_ids(data: _ChainData, numbered: list, first_row: int) -> dict:
    if not numbered:
        raise ValueError("ANARCI returned no numbered residues.")
    mapping = {}
    first_number = numbered[0][0][0]
    for query_index in range(first_row):
        mapping[data.residue_indices[query_index]] = (
            first_number - (first_row - query_index),
            "",
        )

    numbered_index = 0
    last_number = first_number - 1
    for query_index in range(first_row, len(data.sequence)):
        if numbered_index < len(numbered):
            (number, insertion_code), expected = numbered[numbered_index]
            actual = data.sequence[query_index]
            if expected != actual:
                original = data.residue_ids[query_index]
                raise ValueError(
                    f"Residue mismatch at {original[0]}{original[1]}: "
                    f"ANARCI expected {expected}, structure contains {actual}."
                )
            last_number = number
            numbered_index += 1
        else:
            last_number += 1
            number, insertion_code = last_number, ""
        mapping[data.residue_indices[query_index]] = (number, insertion_code)
    return mapping


def _check_collisions(structure, chain: str, mapping: dict) -> None:
    target = _find_chain(structure, chain)
    is_biopython = isinstance(structure, BioStructure)
    final_ids = []
    for index, residue in enumerate(target):
        if is_biopython:
            if index in mapping:
                number, insertion_code = mapping[index]
                final_id = (residue.id[0], number, insertion_code or " ")
            else:
                final_id = residue.id
        else:
            if index in mapping:
                number, insertion_code = mapping[index]
                final_id = (residue.het_flag, number, insertion_code or " ")
            else:
                final_id = (
                    residue.het_flag,
                    residue.seqid.num,
                    residue.seqid.icode,
                )
        final_ids.append(final_id)
    if len(final_ids) != len(set(final_ids)):
        raise ValueError(
            "Renumbering would collide with an unchanged residue ID outside "
            "the selected range. Renumber the whole chain or choose a "
            "non-overlapping range."
        )


def apply_numbering(
    structure,
    chain: str,
    data: _ChainData,
    numbered: list,
    first_row: int,
):
    """Return a same-type clone with numbering applied to selected residues."""
    mapping = _new_residue_ids(data, numbered, first_row)
    _check_collisions(structure, chain, mapping)

    if isinstance(structure, BioStructure):
        result = copy.deepcopy(structure)
        target = _find_chain(result, chain)
        residues = list(target)
        for residue in residues:
            residue.detach_parent()
        for index, (number, insertion_code) in mapping.items():
            residue = residues[index]
            residue.id = (residue.id[0], number, insertion_code or " ")
        target.child_list = residues
        target.child_dict = {}
        for residue in residues:
            residue.set_parent(target)
            target.child_dict[residue.id] = residue
        return result

    result = structure.clone()
    target = _find_chain(result, chain)
    for index, (number, insertion_code) in mapping.items():
        if len(insertion_code.strip()) > 1:
            raise ValueError(
                "Gemmi structures cannot represent extended insertion codes; "
                "use a BioPython structure and mmCIF output."
            )
        target[index].seqid.num = number
        target[index].seqid.icode = insertion_code.strip() or " "
    return result
