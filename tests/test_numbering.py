import json
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser

from sabr import constants, renumber_structure
from sabr.numbering import alignment_to_states, number_alignment

DATA = Path(__file__).parent / "data"


def test_auto_selection_and_all_numbering_schemes_match_goldens(aligned_cases):
    expected = json.loads((DATA / "numbering_baseline.json").read_text())
    for chain_type, result in aligned_cases.items():
        assert result["selected"] == chain_type
        for scheme in constants.NUMBERING_SCHEMES:
            numbered, first_row = number_alignment(
                result["alignment"],
                result["data"].sequence,
                scheme,
                chain_type,
            )
            golden = expected[chain_type]["schemes"][scheme]
            actual = [
                [number, insertion_code, amino_acid]
                for (number, insertion_code), amino_acid in numbered
            ]
            assert first_row == golden["first_row"]
            assert actual == golden["numbered"]


def test_alignment_state_conversion_preserves_orphan_insertions():
    alignment = np.zeros((5, 8), dtype=int)
    alignment[0, 1] = 1
    alignment[1, 2] = 1
    alignment[4, 5] = 1
    states, start, first_row = alignment_to_states(alignment)
    assert start == 1
    assert first_row == 0
    assert ((3, "i"), 3) in states
    assert ((3, "i"), 4) in states
    assert ((4, "d"), None) in states


def test_long_cdr3_retains_extended_imgt_insertion_codes():
    structure = PDBParser(QUIET=True).get_structure("long", DATA / "8sve_L.pdb")
    with pytest.warns(UserWarning, match="structural gap"):
        result = renumber_structure(
            structure,
            "M",
            scheme="imgt",
            chain_type="K",
        )
    insertion_codes = [
        residue.id[2].strip()
        for residue in result[0]["M"]
        if not residue.id[0].strip()
    ]
    assert any(len(code) > 1 for code in insertion_codes)
