"""Convert an IMGT alignment into one of ANARCI's numbering schemes."""

import logging

import numpy as np

from sabr._anarci import schemes

LOGGER = logging.getLogger(__name__)


def alignment_to_states(matrix: np.ndarray) -> tuple:
    """Return ANARCI-compatible states and the first aligned query row."""
    if matrix.ndim != 2:
        raise ValueError("Alignment matrix must be two-dimensional.")
    path = sorted(np.argwhere(matrix.T == 1).tolist())
    if not path:
        raise ValueError("Alignment matrix contains no aligned residues.")

    column_rows = {}
    row_columns = {}
    for column, row in path:
        column_rows.setdefault(column, []).append(row)
        row_columns[row] = column

    first_column, first_row = path[0]
    last_column, last_row = path[-1]
    orphan_rows = {
        row for row in range(first_row, last_row + 1) if row not in row_columns
    }
    offset = first_column - first_row
    states = []

    for column in range(first_column, last_column + 1):
        imgt_position = column + 1
        if column not in column_rows:
            states.append(((imgt_position, "d"), None))
            continue

        rows = column_rows[column]
        states.append(((imgt_position, "m"), rows[0] + offset))
        for row in rows[1:]:
            states.append(((imgt_position, "i"), row + offset))

        next_row = None
        for later_column in range(column + 1, last_column + 1):
            if later_column in column_rows:
                next_row = column_rows[later_column][0]
                break
        if next_row is not None:
            for row in range(rows[-1] + 1, next_row):
                if row in orphan_rows:
                    states.append(((imgt_position, "i"), row + offset))

    return states, first_column, first_row


def _apply_scheme(states: list, sequence: str, scheme: str, chain_type: str):
    if scheme == "imgt":
        return schemes.number_imgt(states, sequence)
    if scheme == "aho":
        return schemes.number_aho(states, sequence, chain_type)

    heavy = chain_type == "H"
    functions = {
        "chothia": (
            schemes.number_chothia_heavy
            if heavy
            else schemes.number_chothia_light
        ),
        "kabat": (
            schemes.number_kabat_heavy if heavy else schemes.number_kabat_light
        ),
        "martin": (
            schemes.number_martin_heavy
            if heavy
            else schemes.number_martin_light
        ),
        "wolfguy": (
            schemes.number_wolfguy_heavy
            if heavy
            else schemes.number_wolfguy_light
        ),
    }
    return functions[scheme](states, sequence)


def number_alignment(
    alignment: np.ndarray,
    sequence: str,
    scheme: str,
    chain_type: str,
) -> tuple:
    """Return explicit query-row-to-number records for one alignment."""
    states, imgt_start, first_row = alignment_to_states(alignment)
    padded_sequence = "-" * imgt_start + sequence[first_row:]
    numbered, start, end = _apply_scheme(
        states, padded_sequence, scheme, chain_type
    )
    numbered = [entry for entry in numbered if entry[1] != "-"]
    if not 0 <= start <= end < len(padded_sequence):
        raise ValueError(
            f"ANARCI returned invalid sequence bounds {start}:{end}."
        )

    padded_rows = [None] * imgt_start + list(range(first_row, len(sequence)))
    query_rows = [
        row for row in padded_rows[start : end + 1] if row is not None
    ]
    if len(query_rows) != len(numbered):
        raise ValueError(
            "ANARCI numbering length does not match its query-row span."
        )

    records = []
    for query_row, ((number, insertion_code), amino_acid) in zip(
        query_rows, numbered
    ):
        if sequence[query_row] != amino_acid:
            raise ValueError(
                f"ANARCI expected {amino_acid} at query row {query_row}, "
                f"but the structure contains {sequence[query_row]}."
            )
        records.append((query_row, number, insertion_code, amino_acid))
    LOGGER.info(
        "Numbered %d residues as %s using %s.",
        len(records),
        chain_type,
        scheme,
    )
    return records
