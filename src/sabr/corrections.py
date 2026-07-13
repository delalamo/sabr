"""Deterministic CDR, light-chain DE-loop, and C-terminal corrections."""

import logging
import warnings

import numpy as np

from sabr import constants

LOGGER = logging.getLogger(__name__)


def _has_gap_in_region(
    gap_indices: frozenset[int], start_row: int, end_row: int
) -> bool:
    return any(index in gap_indices for index in range(start_row, end_row))


def _aligned_row_near(aln: np.ndarray, target_col: int) -> int | None:
    """Return the row aligned at or within two columns of a target."""
    for offset in (0, -1, 1, -2, 2):
        column = target_col + offset
        if 0 <= column < aln.shape[1]:
            rows = np.where(aln[:, column] == 1)[0]
            if len(rows):
                return int(rows[0])

    return None


def cdr_columns(residue_count: int, position_count: int) -> list:
    """Return IMGT's outside-in CDR columns, preserving central insertions."""
    columns = [None] * residue_count
    for index in range(min(residue_count, position_count)):
        if index % 2 == 0:
            row = column = index // 2
        else:
            offset = (index + 1) // 2
            row = residue_count - offset
            column = position_count - offset
        columns[row] = column
    return columns


def _skip_for_structural_gap(
    gap_indices: frozenset[int] | None,
    start_row: int,
    end_row: int,
    region_name: str,
) -> bool:
    """Warn and return true when a regional correction crosses a gap."""
    if gap_indices and _has_gap_in_region(gap_indices, start_row, end_row):
        message = (
            f"Skipping {region_name} deterministic correction: structural "
            f"gap detected between rows {start_row} and {end_row}; using "
            "the raw embedding alignment for this region."
        )
        warnings.warn(message, UserWarning, stacklevel=3)
        return True
    return False


def _move_or_clear_position(
    aln: np.ndarray,
    source_col: int,
    target_col: int,
    source_pos: str,
    target_pos: str,
) -> None:
    """Move a residue to an empty target column, then clear its source."""
    target_occupied = aln[:, target_col].sum() == 1
    if not target_occupied:
        LOGGER.info(
            f"Moving residue from position {source_pos} to position "
            f"{target_pos} (chain lacks position {source_pos})"
        )
        aln[:, target_col] = aln[:, source_col]
    else:
        LOGGER.info(
            f"Clearing position {source_pos} (chain lacks position "
            f"{source_pos}, but position {target_pos} already occupied)"
        )
    aln[:, source_col] = 0


def correct_light_chain_de_loop(
    aln: np.ndarray,
    gap_indices: frozenset[int] | None = None,
) -> np.ndarray:
    """Move light-chain residues aligned at 81/82 to IMGT 83/84."""
    pos81_col = constants.FR3_POS81_COL
    pos82_col = constants.FR3_POS82_COL
    pos83_col = constants.FR3_POS83_COL
    pos84_col = constants.FR3_POS84_COL

    if gap_indices:
        de_start_col = 78  # 0-indexed for position 79
        de_end_col = 83  # 0-indexed for position 84
        de_start_row = _aligned_row_near(aln, de_start_col)
        de_end_row = _aligned_row_near(aln, de_end_col)
        if (
            de_start_row is not None
            and de_end_row is not None
            and _skip_for_structural_gap(
                gap_indices, de_start_row, de_end_row, "FR3/DE loop"
            )
        ):
            return aln

    if aln[:, pos81_col].sum() == 1:
        _move_or_clear_position(aln, pos81_col, pos83_col, "81", "83")

    if aln[:, pos82_col].sum() == 1:
        _move_or_clear_position(aln, pos82_col, pos84_col, "82", "84")

    return aln


def correct_c_terminus(aln: np.ndarray) -> np.ndarray:
    """Assign trailing residues after IMGT 125, stopping at IMGT 128."""
    n_rows, n_cols = aln.shape
    assigned_rows = np.flatnonzero(aln.sum(axis=1))
    if not len(assigned_rows):
        return aln
    last_assigned_row = assigned_rows[-1]
    last_assigned_col = np.flatnonzero(aln.sum(axis=0))[-1]
    n_unassigned_trailing = n_rows - last_assigned_row - 1
    if n_unassigned_trailing <= 0:
        return aln
    if last_assigned_col < constants.C_TERMINUS_ANCHOR_POSITION:
        return aln

    LOGGER.info(
        f"Correcting C-terminus: {n_unassigned_trailing} unassigned "
        f"residues after row {last_assigned_row}, "
        f"last assigned col was {last_assigned_col}"
    )

    rows = range(last_assigned_row + 1, n_rows)
    columns = range(last_assigned_col + 1, n_cols)
    for row, column in zip(rows, columns):
        aln[row, :] = 0
        aln[row, column] = 1

    return aln


def correct_cdr_loop(
    aln: np.ndarray,
    loop_name: str,
    cdr_start: int,
    cdr_end: int,
    gap_indices: frozenset[int] | None = None,
) -> np.ndarray:
    """Place framework and CDR residues between one pair of anchors."""
    anchor_start, anchor_end = constants.CDR_ANCHORS[loop_name]
    anchor_start_col = anchor_start - 1
    anchor_end_col = anchor_end - 1

    anchor_start_row = _aligned_row_near(aln, anchor_start_col)
    anchor_end_row = _aligned_row_near(aln, anchor_end_col)
    if anchor_start_row is None or anchor_end_row is None:
        LOGGER.warning(
            "Skipping %s; missing anchor near IMGT %d or %d.",
            loop_name,
            anchor_start,
            anchor_end,
        )
        return aln

    if anchor_start_row >= anchor_end_row:
        LOGGER.warning(
            f"Skipping {loop_name}; anchor start row "
            f"({anchor_start_row}) >= end row ({anchor_end_row})"
        )
        return aln

    if _skip_for_structural_gap(
        gap_indices, anchor_start_row, anchor_end_row, loop_name
    ):
        return aln

    fw_before_cdr = list(range(anchor_start + 1, cdr_start))
    fw_after_cdr = list(range(cdr_end + 1, anchor_end))
    n_fw_before = len(fw_before_cdr)
    n_fw_after = len(fw_after_cdr)

    intermediate_rows = list(range(anchor_start_row + 1, anchor_end_row))
    n_residues = len(intermediate_rows)

    if n_residues < n_fw_before + n_fw_after:
        LOGGER.warning(
            f"Skipping {loop_name}; not enough residues "
            f"({n_residues}) between anchors for FW positions "
            f"({n_fw_before} + {n_fw_after})"
        )
        return aln

    n_cdr_residues = n_residues - n_fw_before - n_fw_after

    LOGGER.info(
        f"{loop_name}: anchors at {anchor_start} (row "
        f"{anchor_start_row}) and {anchor_end} (row "
        f"{anchor_end_row}). {n_residues} residues: "
        f"{n_fw_before} FW, {n_cdr_residues} CDR, {n_fw_after} FW"
    )

    aln[intermediate_rows, :] = 0

    for row, position in zip(intermediate_rows, fw_before_cdr):
        aln[row, position - 1] = 1

    fw_after_rows = intermediate_rows[-n_fw_after:] if n_fw_after > 0 else []
    for row, position in zip(fw_after_rows, fw_after_cdr):
        aln[row, position - 1] = 1

    cdr_rows = intermediate_rows[n_fw_before:]
    if n_fw_after > 0:
        cdr_rows = cdr_rows[:-n_fw_after]

    position_count = cdr_end - cdr_start + 1
    for row, column in zip(
        cdr_rows, cdr_columns(len(cdr_rows), position_count)
    ):
        if column is not None:
            aln[row, cdr_start - 1 + column] = 1

    return aln


def apply_corrections(
    aln: np.ndarray,
    chain_type: str,
    gap_indices: frozenset[int] | None = None,
) -> np.ndarray:
    """Apply CDR, light-chain DE-loop, then C-terminal corrections."""
    for loop_name, (cdr_start, cdr_end) in constants.IMGT_LOOPS.items():
        aln = correct_cdr_loop(
            aln, loop_name, cdr_start, cdr_end, gap_indices=gap_indices
        )

    if chain_type in ("K", "L"):
        aln = correct_light_chain_de_loop(aln, gap_indices=gap_indices)

    aln = correct_c_terminus(aln)

    return aln
