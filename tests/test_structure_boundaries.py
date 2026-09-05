import copy

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, PDBParser, Residue
from Bio.PDB.Atom import DisorderedAtom
from Bio.PDB.Residue import DisorderedResidue

import sabr.api as api
from sabr import constants, renumber_structure
from sabr.structure import _detect_gaps, extract_chain
from tests.helpers import DATA, assert_atomic_content_equal, residue_ids


@pytest.fixture
def heavy():
    return PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )


@pytest.fixture
def forbid_encoder(monkeypatch):
    monkeypatch.setattr(
        api,
        "encode",
        lambda *args, **kwargs: pytest.fail("invalid input reached encoder"),
    )


@pytest.mark.parametrize("value", (float("nan"), float("inf"), -float("inf")))
def test_nonfinite_coordinates_fail_before_encoding(
    heavy, forbid_encoder, value
):
    next(heavy.get_atoms()).coord[0] = value
    with (
        np.errstate(invalid="raise", divide="raise"),
        pytest.raises(ValueError, match="non-finite coordinates"),
    ):
        renumber_structure(heavy, "F")


@pytest.mark.parametrize("empty", ("range", "chain", "models"))
def test_empty_inputs_fail_before_encoding(heavy, forbid_encoder, empty):
    kwargs = {}
    if empty == "range":
        kwargs["residue_range"] = (500, 600)
    elif empty == "chain":
        for residue in list(heavy[0]["F"]):
            heavy[0]["F"].detach_child(residue.id)
    else:
        heavy.detach_child(0)
    with pytest.raises(ValueError, match="No polymer|exactly one model"):
        renumber_structure(heavy, "F", **kwargs)


@pytest.mark.parametrize("ambiguous", ("disordered", "duplicate"))
def test_ambiguous_polymer_identity_fails_before_encoding(
    heavy, forbid_encoder, ambiguous
):
    first = next(heavy[0]["F"].get_residues())
    other = copy.deepcopy(first)
    other.detach_parent()
    other.resname = "SER"
    if ambiguous == "duplicate":
        other.id = ("H_SER", first.id[1], first.id[2])
        heavy[0]["F"].add(other)
    else:
        heavy[0]["F"].detach_child(first.id)
        disordered = DisorderedResidue(first.id)
        disordered.disordered_add(first)
        disordered.disordered_add(other)
        heavy[0]["F"].add(disordered)
    with pytest.raises(ValueError, match="ambiguous.*microheterogeneous"):
        renumber_structure(heavy, "F")


def test_modified_hetero_residue_is_sequenced_without_chemistry_changes(heavy):
    first = next(heavy.get_residues())
    first.resname = "MSE"
    first.id = ("H_MSE", first.id[1], first.id[2])
    original = copy.deepcopy(heavy)
    data = extract_chain(heavy, "F", (2, 2))
    assert data.sequence == "M"
    assert_atomic_content_equal(original, heavy)
    assert residue_ids(original, "F") == residue_ids(heavy, "F")


def test_selection_limit_accepts_boundary_and_rejects_next_residue(
    heavy, monkeypatch, forbid_encoder
):
    monkeypatch.setattr(constants, "MAX_SELECTED_RESIDUES", 3)
    assert len(extract_chain(heavy, "F", (2, 4)).sequence) == 3
    with pytest.raises(ValueError, match="safety limit"):
        renumber_structure(heavy, "F", residue_range=(2, 5))


@pytest.mark.parametrize(
    "distance, expected",
    [(2.659999, frozenset()), (2.66, frozenset()), (2.660001, frozenset({0}))],
)
def test_gap_threshold_is_strict_and_reports_left_residue(distance, expected):
    coords = np.zeros((2, 4, 3), dtype=np.float64)
    coords[1, 0, 0] = distance
    assert _detect_gaps(coords) == expected
    assert _detect_gaps(coords[:1]) == frozenset()


def test_full_copy_preserves_atoms_conformers_metadata_and_partial_range(
    heavy, monkeypatch
):
    # Number a middle span into a separate block; preserve both sides.
    first = list(heavy[0]["F"])[10]
    ca = first["CA"]
    first.detach_child("CA")
    disordered = DisorderedAtom("CA")
    for name, occupancy in (("B", 0.6), ("A", 0.4)):
        atom = copy.deepcopy(ca)
        atom.detach_parent()
        atom.altloc = name
        atom.occupancy = occupancy
        disordered.disordered_add(atom)
    first.add(disordered)
    first.xtra = {"residue": {"tags": ["keep"]}}
    heavy.header["audit"] = {"nested": ["keep"]}
    heavy.xtra = {"structure": {"tags": ["keep"]}}
    next(heavy.get_atoms()).xtra = {"atom": ["keep"]}
    next(heavy.get_atoms()).set_anisou(np.arange(6, dtype=float))
    for chain_id, hetero, name in (("Z", "W", "HOH"), ("F", "H_LIG", "LIG")):
        if chain_id not in heavy[0]:
            heavy[0].add(Chain.Chain(chain_id))
        residue = Residue.Residue((hetero, 900, " "), name, " ")
        residue.add(
            Atom.Atom(
                "O",
                np.array([1.0, 2.0, 3.0]),
                12.34,
                0.75,
                " ",
                " O  ",
                9000,
                element="O",
            )
        )
        heavy[0][chain_id].add(residue)
    before = copy.deepcopy(heavy)
    monkeypatch.setattr(
        api, "encode", lambda coords, mode: np.zeros((len(coords), 64))
    )
    monkeypatch.setattr(api, "align", lambda *args, **kwargs: (None, "H", 0.0))
    monkeypatch.setattr(
        api,
        "number_alignment",
        lambda alignment, sequence, scheme, kind: [
            (i, 200 + i, "", aa) for i, aa in enumerate(sequence)
        ],
    )
    result = renumber_structure(heavy, "F", residue_range=(5, 120))
    assert_atomic_content_equal(before, result)
    assert_atomic_content_equal(before, heavy)
    assert heavy.header == result.header == before.header
    assert heavy.xtra == result.xtra == before.xtra
    assert residue_ids(heavy, "F") == residue_ids(before, "F")
    for old, new in zip(before[0]["F"], result[0]["F"], strict=True):
        if not 5 <= old.id[1] <= 120:
            assert old.id == new.id
        assert new.get_parent() is result[0]["F"]
        assert result[0]["F"][new.id] is new
    assert len({r.id for r in result[0]["F"]}) == len(result[0]["F"])
    next(result.get_atoms()).coord[0] += 100
    next(result.get_atoms()).xtra["atom"].append("changed")
    result.header["audit"]["nested"].append("changed")
    result.xtra["structure"]["tags"].append("changed")
    next(r for r in result[0]["F"] if r.xtra).xtra["residue"]["tags"].append(
        "changed"
    )
    assert_atomic_content_equal(before, heavy)
    assert heavy.header == before.header
    assert heavy.xtra == before.xtra
