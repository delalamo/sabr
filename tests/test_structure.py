import copy
from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from sabr import renumber_structure
from sabr.structure import extract_chain

DATA = Path(__file__).parent / "data"


def test_numeric_residue_range_includes_insertion_codes():
    structure = PDBParser(QUIET=True).get_structure(
        "insertions", DATA / "test_insertion_codes.pdb"
    )
    data = extract_chain(structure, "A", (52, 52))
    assert data.residue_ids == ((52, ""), (52, "A"), (52, "B"))


def test_missing_backbone_is_reported_with_the_residue():
    structure = PDBParser(QUIET=True).get_structure(
        "missing", DATA / "test_no_backbone.pdb"
    )
    with pytest.raises(ValueError, match="missing required backbone"):
        renumber_structure(structure, "A")


@pytest.mark.parametrize(
    ("chain", "message"),
    [("", "non-empty"), ("missing", "not found")],
)
def test_invalid_or_missing_chain_is_reported(chain, message):
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    with pytest.raises(ValueError, match=message):
        renumber_structure(structure, chain)


def test_missing_backbone_outside_selected_range_is_ignored():
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    list(structure[0]["F"])[-1].detach_child("N")
    data = extract_chain(structure, "F", (2, 127))
    assert data.residue_ids[0] == (2, "")
    assert data.residue_ids[-1] == (127, "")


def test_multiple_models_are_rejected():
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    second = copy.deepcopy(structure[0])
    second.id = 1
    structure.add(second)
    with pytest.raises(ValueError, match="exactly one model"):
        renumber_structure(structure, "F")


def test_unsupported_structure_type_raises_type_error():
    with pytest.raises(TypeError, match="Bio.PDB"):
        renumber_structure(object(), "A")
