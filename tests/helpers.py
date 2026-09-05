"""Assertions shared by integration tests; no expected values are generated."""

from pathlib import Path

import numpy as np

DATA = Path(__file__).parent / "data"
CASES = {
    "H": ("test_heavy_chain.pdb", "F"),
    "K": ("12e8_L.pdb", "L"),
    "L": ("1bjm_A.pdb", "A"),
}


def residue_ids(structure, chain):
    return [(r.id[1], r.id[2].strip()) for r in structure[0][chain]]


def assert_atomic_content_equal(before, after, *, serialized=False):
    """Compare every atom/conformer in hierarchy order, excluding residue IDs.

    PDB/mmCIF writers normalize atom serials/full names and round coordinates
    and scalar fields. In-memory copies must preserve these values exactly.
    """
    assert [m.id for m in before] == [m.id for m in after]
    assert [c.id for c in before.get_chains()] == [
        c.id for c in after.get_chains()
    ]
    left = [r for c in before.get_chains() for r in c.get_unpacked_list()]
    right = [r for c in after.get_chains() for r in c.get_unpacked_list()]
    assert len(left) == len(right)
    for a, b in zip(left, right, strict=True):
        assert (a.resname, a.id[0]) == (b.resname, b.id[0])
        if not serialized:
            assert a.segid == b.segid
            assert a.xtra == b.xtra
        aa, bb = a.get_unpacked_list(), b.get_unpacked_list()
        assert len(aa) == len(bb)
        for x, y in zip(aa, bb, strict=True):
            assert (x.name, x.altloc, x.element) == (
                y.name,
                y.altloc,
                y.element,
            )
            assert (x.occupancy is None) == (y.occupancy is None)
            np.testing.assert_allclose(
                x.coord, y.coord, rtol=0, atol=0.00051 if serialized else 0
            )
            if x.occupancy is not None:
                np.testing.assert_allclose(
                    x.occupancy,
                    y.occupancy,
                    rtol=0,
                    atol=0.0051 if serialized else 0,
                )
            np.testing.assert_allclose(
                x.bfactor, y.bfactor, rtol=0, atol=0.0051 if serialized else 0
            )
            if not serialized:
                assert (x.fullname, x.serial_number) == (
                    y.fullname,
                    y.serial_number,
                )
                assert x.xtra == y.xtra
                for getter in ("get_anisou", "get_sigatm", "get_siguij"):
                    np.testing.assert_equal(
                        getattr(x, getter)(), getattr(y, getter)()
                    )
