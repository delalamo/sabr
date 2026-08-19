import copy
from pathlib import Path

import gemmi
import numpy as np
import pytest
from Bio.PDB import Atom, Chain, PDBParser, Residue
from Bio.PDB.Structure import Structure as BioStructure

import sabr.api as api
from sabr import renumber_structure
from sabr.structure import extract_chain

DATA = Path(__file__).parent / "data"


def _bio_ids(structure, chain):
    return [
        (residue.id[1], residue.id[2].strip())
        for residue in structure[0][chain]
        if not residue.id[0].strip()
    ]


def _gemmi_ids(structure, chain):
    target = next(
        candidate for candidate in structure[0] if candidate.name == chain
    )
    return [
        (residue.seqid.num, residue.seqid.icode.strip())
        for residue in target
        if residue.het_flag == "A"
    ]


def _stub_pipeline(monkeypatch, captured_modes=None):
    def fake_encode(coords, mode):
        if captured_modes is not None:
            captured_modes["encoder"] = mode
        return np.zeros((len(coords), 64))

    def fake_align(embeddings, gaps, chain_type, noise, mode, scfv):
        if captured_modes is not None:
            captured_modes["alignment"] = mode
        return (
            np.zeros((len(embeddings), 128), dtype=int),
            "H" if chain_type == "auto" else chain_type,
            0.0,
        )

    monkeypatch.setattr(
        api,
        "encode",
        fake_encode,
    )
    monkeypatch.setattr(
        api,
        "align",
        fake_align,
    )
    monkeypatch.setattr(
        api,
        "number_alignment",
        lambda alignment, sequence, scheme, chain_type: (
            [
                (index, index + 1, "", amino_acid)
                for index, amino_acid in enumerate(sequence)
            ]
        ),
    )


def test_biopython_and_gemmi_are_non_mutating_and_equivalent():
    path = DATA / "test_heavy_chain.pdb"
    bio = PDBParser(QUIET=True).get_structure("bio", path)
    gemmi_structure = gemmi.read_structure(str(path))
    original_bio = _bio_ids(bio, "F")
    original_gemmi = _gemmi_ids(gemmi_structure, "F")

    bio_result = renumber_structure(bio, "F", chain_type="H")
    gemmi_result = renumber_structure(gemmi_structure, "F", chain_type="H")

    assert isinstance(bio_result, BioStructure)
    assert isinstance(gemmi_result, gemmi.Structure)
    assert _bio_ids(bio, "F") == original_bio
    assert _gemmi_ids(gemmi_structure, "F") == original_gemmi
    assert _bio_ids(bio_result, "F") == _gemmi_ids(gemmi_result, "F")


def test_non_target_content_and_out_of_range_residues_are_preserved(
    monkeypatch,
):
    _stub_pipeline(monkeypatch)
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    other = Chain.Chain("Z")
    other_residue = Residue.Residue(("W", 1, " "), "HOH", " ")
    other_residue.add(
        Atom.Atom(
            "O",
            np.asarray((1.0, 2.0, 3.0)),
            1.0,
            1.0,
            " ",
            " O  ",
            1,
            element="O",
        )
    )
    other.add(other_residue)
    structure[0].add(other)
    original_other = copy.deepcopy(structure[0]["Z"])
    last_id = list(structure[0]["F"])[-1].id

    result = renumber_structure(
        structure,
        "F",
        chain_type="H",
        residue_range=(2, 127),
    )

    assert list(result[0]["F"])[-1].id == last_id
    assert list(structure[0]["F"])[-1].id == last_id
    assert result[0]["Z"][other_residue.id]["O"].coord.tolist() == (
        original_other[other_residue.id]["O"].coord.tolist()
    )


def test_numeric_residue_range_includes_insertion_codes():
    structure = PDBParser(QUIET=True).get_structure(
        "insertions", DATA / "test_insertion_codes.pdb"
    )
    data = extract_chain(structure, "A", (52, 52))
    assert data.residue_ids == ((52, ""), (52, "A"), (52, "B"))


def test_range_that_would_collide_with_unchanged_ids_is_rejected(
    monkeypatch,
):
    _stub_pipeline(monkeypatch)
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    unchanged = copy.deepcopy(list(structure[0]["F"])[0])
    unchanged.detach_parent()
    unchanged.id = (" ", 1, " ")
    structure[0]["F"].add(unchanged)
    with pytest.raises(ValueError, match="collide"):
        renumber_structure(
            structure,
            "F",
            chain_type="H",
            residue_range=(2, 128),
        )


def test_arbitrary_gemmi_chain_ids_are_supported(monkeypatch):
    _stub_pipeline(monkeypatch)
    structure = gemmi.read_structure(str(DATA / "test_heavy_chain.pdb"))
    structure[0][0].name = "heavy_chain"
    result = renumber_structure(structure, "heavy_chain", chain_type="H")
    assert result[0][0].name == "heavy_chain"


def test_softalign_mode_is_forwarded_to_encoder_and_alignment(monkeypatch):
    captured_modes = {}
    _stub_pipeline(monkeypatch, captured_modes)
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    renumber_structure(structure, "F", mode="softalign")
    assert captured_modes == {
        "encoder": "softalign",
        "alignment": "softalign",
    }


def test_chain_type_choices_come_from_reference_asset(monkeypatch):
    monkeypatch.setattr(
        api,
        "load_references",
        lambda noise_level, mode: {"X": (None, None)},
    )
    options = api._RenumberOptions(
        scheme="imgt",
        chain_type="x",
        noise_level=0.0,
        residue_range=None,
        mode="sabr",
        scfv=False,
    )
    assert options.chain_type == "X"

    with pytest.raises(ValueError, match="auto, X"):
        api._RenumberOptions(
            scheme="imgt",
            chain_type="H",
            noise_level=0.0,
            residue_range=None,
            mode="sabr",
            scfv=False,
        )


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


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"scheme": "unknown"}, "scheme"),
        ({"chain_type": "heavy"}, "chain_type"),
        ({"noise_level": 0.3}, "noise_level"),
        ({"residue_range": [1, 2]}, "residue_range"),
        ({"residue_range": (2, 1)}, "start"),
        ({"chain_type": None}, "chain_type"),
        ({"noise_level": "invalid"}, "noise_level"),
        ({"noise_level": False}, "noise_level"),
        ({"mode": "other"}, "mode"),
        ({"mode": None}, "mode"),
        ({"residue_range": (False, 10)}, "residue_range"),
        ({"scfv": "yes"}, "scfv"),
        ({"scfv": True, "chain_type": "H"}, "auto"),
    ],
)
def test_invalid_public_options_fail_before_model_execution(kwargs, message):
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    with pytest.raises(ValueError, match=message):
        renumber_structure(structure, "F", **kwargs)


def test_unsupported_structure_type_raises_type_error():
    with pytest.raises(TypeError, match="Bio.PDB"):
        renumber_structure(object(), "A")
