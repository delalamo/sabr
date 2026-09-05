"""Capture candidate regression data manually; never imported by pytest.

Run in the pinned CPU environment. Output must be a new directory. Review
candidate mappings before copying individual files into tests/data.
"""

import argparse
import hashlib
import importlib.metadata
import io
import json
import platform
import subprocess
from pathlib import Path
from urllib.request import urlopen

import jax
import numpy as np
from Bio.PDB import MMCIFIO, MMCIFParser, PDBParser, Select

from sabr import constants, renumber_structure
from sabr.alignment import _align_reference, align, load_references
from sabr.model import encode
from sabr.numbering import number_alignment
from sabr.structure import extract_chain


class ChainSelection(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        return chain.id == self.chain


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if jax.default_backend() != "cpu":
        parser.error("capture regressions with JAX_PLATFORMS=cpu")
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=False)
    data_dir = args.repo / "tests" / "data"
    manifest = {
        "producing_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=args.repo, text=True
        ).strip(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": {
            name: importlib.metadata.version(name)
            for name in ("biopython", "jax", "jaxlib", "numpy", "dm-haiku")
        },
        "backend": jax.default_backend(),
        "parameters": {"noise_level": 0.0, "scfv_schemes": ["imgt"]},
        "fixtures": {},
        "inputs": {},
        "review": "Candidate captures; require review before acceptance.",
    }
    soft_numbering = {}
    cases = {
        "H": ("test_heavy_chain.pdb", "F"),
        "K": ("12e8_L.pdb", "L"),
        "L": ("1bjm_A.pdb", "A"),
    }
    for kind, (filename, chain) in cases.items():
        path = data_dir / filename
        manifest["inputs"][filename] = digest(path)
        structure = PDBParser(QUIET=True).get_structure(kind, path)
        data = extract_chain(structure, chain, None)
        embeddings = encode(data.coords, "softalign")
        alignment, selected, _ = align(
            embeddings, data.gap_indices, "auto", 0.0, mode="softalign"
        )
        assert selected == kind
        soft_numbering[kind] = {"schemes": {}}
        for scheme in constants.NUMBERING_SCHEMES:
            numbered = number_alignment(alignment, data.sequence, scheme, kind)
            soft_numbering[kind]["schemes"][scheme] = {
                "first_row": numbered[0][0],
                "numbered": [
                    [number, code, aa] for _, number, code, aa in numbered
                ],
            }
        if kind == "H":
            reference, positions = load_references(0.0, "softalign")[kind]
            reduced, similarity, score = _align_reference(
                embeddings, reference, "softalign"
            )
            np.savez_compressed(
                output / "softalign_math_baseline.npz",
                embeddings=embeddings,
                reduced_alignment=reduced,
                similarity=similarity,
                score=score,
                positions=positions,
            )
    save_json(output / "softalign_numbering_baseline.json", soft_numbering)
    multi = {}
    for code, chain in (("8DY0", "A"), ("8DY1", "C")):
        url = f"https://files.rcsb.org/download/{code}.cif"
        with urlopen(url, timeout=60) as response:
            raw = response.read()
        structure = MMCIFParser(QUIET=True).get_structure(
            code, io.StringIO(raw.decode())
        )
        path = output / f"{code.lower()}_{chain}.cif"
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(str(path), ChainSelection(chain))
        extracted = MMCIFParser(QUIET=True).get_structure(code, path)
        data = extract_chain(extracted, chain, None)
        manifest["fixtures"][path.name] = {
            "source": url,
            "entry": code,
            "author_chain": chain,
            "source_sha256": hashlib.sha256(raw).hexdigest(),
            "extracted_sha256": digest(path),
            "extraction": (
                "Biopython MMCIFIO, first/only model, selected author chain "
                "including hetero residues and all atom conformers; "
                "atomic categories only."
            ),
            "polymer_residues": len(data.sequence),
            "gap_indices": sorted(data.gap_indices),
        }
        multi[path.stem] = {}
        for mode in constants.MODES:
            embeddings = encode(data.coords, mode)
            _, selected, score = align(
                embeddings, data.gap_indices, "auto", 0.0, mode=mode, scfv=True
            )
            result = renumber_structure(
                extracted,
                chain,
                mode=mode,
                scfv=True,
                dangerously_allow_structural_gaps=bool(data.gap_indices),
            )
            original = [
                list(extracted[0][chain])[i] for i in data.residue_indices
            ]
            numbered = [list(result[0][chain])[i] for i in data.residue_indices]
            mapping = [
                [a.id[1], a.id[2].strip(), b.id[1], b.id[2].strip(), a.resname]
                for a, b in zip(original, numbered, strict=True)
            ]
            multi[path.stem][mode] = {
                "selected": selected,
                "score": score,
                "mapping": mapping,
            }
            print(
                path.stem,
                mode,
                selected,
                len(mapping),
                mapping[0],
                mapping[-1],
                flush=True,
            )
    save_json(output / "multidomain_baseline.json", multi)
    manifest["outputs"] = {p.name: digest(p) for p in sorted(output.iterdir())}
    save_json(output / "regression_provenance.json", manifest)


if __name__ == "__main__":
    main()
