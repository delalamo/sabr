import json
from pathlib import Path

import gemmi
import numpy as np
import pytest
from Bio.PDB import PDBParser

import sabr.numbering as numbering
from sabr import constants, renumber_structure
from sabr.numbering import alignment_to_states, number_alignment

DATA = Path(__file__).parent / "data"


def test_auto_selection_and_all_numbering_schemes_match_goldens(aligned_cases):
    expected = json.loads((DATA / "numbering_baseline.json").read_text())
    for chain_type, result in aligned_cases.items():
        assert result["selected"] == chain_type
        for scheme in constants.NUMBERING_SCHEMES:
            numbered = number_alignment(
                result["alignment"],
                result["data"].sequence,
                scheme,
                chain_type,
            )
            golden = expected[chain_type]["schemes"][scheme]
            actual = [
                [number, insertion_code, amino_acid]
                for _, number, insertion_code, amino_acid in numbered
            ]
            assert numbered[0][0] == golden["first_row"]
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


def test_8sve_huge_cdr1_matches_the_accepted_full_mapping():
    golden = json.loads((DATA / "8sve_cdr1_imgt.json").read_text())
    structure = PDBParser(QUIET=True).get_structure("long", DATA / "8sve_L.pdb")
    source_residues = [
        residue for residue in structure[0]["M"] if not residue.id[0].strip()
    ]
    with pytest.warns(UserWarning) as warnings:
        result = renumber_structure(
            structure,
            "M",
            scheme="imgt",
            chain_type="K",
        )
    assert [str(item.message) for item in warnings] == golden["warning"]
    target_residues = [
        residue for residue in result[0]["M"] if not residue.id[0].strip()
    ]
    actual = [
        [
            source.id[1],
            source.id[2].strip(),
            target.id[1],
            target.id[2].strip(),
            source.resname,
        ]
        for source, target in zip(source_residues, target_residues, strict=True)
    ]
    assert actual == golden["mapping"]
    assert sum(item[2] == 32 for item in actual) == 69
    assert sum(item[2] == 33 for item in actual) == 70


def test_8sve_gemmi_reports_its_extended_insertion_code_limit():
    structure = gemmi.read_structure(str(DATA / "8sve_L.pdb"))
    with pytest.warns(UserWarning, match="CDR1"):
        with pytest.raises(ValueError, match="extended insertion codes"):
            renumber_structure(structure, "M", scheme="imgt", chain_type="K")


def test_number_alignment_retains_original_query_rows(monkeypatch):
    alignment = np.zeros((5, 8), dtype=int)
    alignment[2, 2] = alignment[3, 3] = alignment[4, 4] = 1

    def fake_scheme(states, sequence, scheme, chain_type):
        assert sequence == "--DEF"
        return [((3, ""), "D"), ((4, ""), "E"), ((5, ""), "F")], 2, 4

    monkeypatch.setattr(numbering, "_apply_scheme", fake_scheme)
    assert number_alignment(alignment, "ACDEF", "imgt", "H") == [
        (2, 3, "", "D"),
        (3, 4, "", "E"),
        (4, 5, "", "F"),
    ]


@pytest.mark.parametrize(
    ("records", "start", "end", "message"),
    [
        ([((1, ""), "A"), ((2, ""), "C")], -1, 0, "bounds"),
        ([((1, ""), "A")], 0, 1, "length"),
        ([((1, ""), "X"), ((2, ""), "C")], 0, 1, "expected"),
    ],
)
def test_number_alignment_rejects_invalid_anarci_output(
    monkeypatch, records, start, end, message
):
    monkeypatch.setattr(
        numbering,
        "_apply_scheme",
        lambda states, sequence, scheme, chain_type: (records, start, end),
    )
    with pytest.raises(ValueError, match=message):
        number_alignment(np.eye(2, dtype=int), "AC", "imgt", "H")
