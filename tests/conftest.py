import functools

import pytest
from Bio.PDB import PDBParser

from sabr.alignment import align
from sabr.model import encode
from sabr.structure import extract_chain
from tests.helpers import CASES, DATA


@pytest.fixture(scope="session")
def aligned_case():
    @functools.cache
    def prepare(chain_type, mode="sabr"):
        filename, chain = CASES[chain_type]
        path = DATA / filename
        structure = PDBParser(QUIET=True).get_structure(chain_type, path)
        data = extract_chain(structure, chain, None)
        embeddings = encode(data.coords, mode)
        alignment, selected, score = align(
            embeddings, data.gap_indices, "auto", 0.0, mode=mode
        )
        return {
            "alignment": alignment,
            "embeddings": embeddings,
            "chain": chain,
            "data": data,
            "path": path,
            "selected": selected,
            "score": score,
        }

    return prepare
