from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from sabr.alignment import align
from sabr.model import encode
from sabr.structure import extract_chain

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def aligned_cases():
    parser = PDBParser(QUIET=True)
    cases = {
        "H": (DATA / "test_heavy_chain.pdb", "F"),
        "K": (DATA / "12e8_L.pdb", "L"),
        "L": (DATA / "1bjm_A.pdb", "A"),
    }
    results = {}
    for expected, (path, chain) in cases.items():
        structure = parser.get_structure(expected, path)
        data = extract_chain(structure, chain, None)
        alignment, selected, score = align(
            encode(data.coords), data.gap_indices, "auto", 0.0
        )
        results[expected] = {
            "alignment": alignment,
            "chain": chain,
            "data": data,
            "path": path,
            "selected": selected,
        }
    return results
