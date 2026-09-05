"""Validate an installed wheel outside the checkout in both parameter modes."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser

import sabr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    args = parser.parse_args()
    checkout = Path(__file__).resolve().parents[1]
    assert (
        not Path(sabr.__file__).resolve().is_relative_to(checkout)
    ), sabr.__file__
    data = args.data_dir.resolve()
    with tempfile.TemporaryDirectory(prefix="sabr-wheel-") as directory:
        working = Path(directory)
        source = PDBParser(QUIET=True).get_structure(
            "input", data / "test_heavy_chain.pdb"
        )
        for residue in source[0]["F"]:
            residue.id = (residue.id[0], residue.id[1] + 500, residue.id[2])
        input_path = working / "input.pdb"
        writer = PDBIO()
        writer.set_structure(source)
        writer.save(str(input_path))
        for filename in (
            "numbering_baseline.json",
            "softalign_numbering_baseline.json",
        ):
            shutil.copyfile(data / filename, working / filename)
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env["JAX_PLATFORMS"] = "cpu"
        for mode in ("sabr", "softalign"):
            output = working / f"{mode}.pdb"
            subprocess.run(
                [
                    str(Path(sys.executable).with_name("sabr")),
                    "-i",
                    str(input_path),
                    "-c",
                    "F",
                    "-o",
                    str(output),
                    "--mode",
                    mode,
                ],
                cwd=working,
                env=env,
                check=True,
            )
            filename = (
                "numbering_baseline.json"
                if mode == "sabr"
                else "softalign_numbering_baseline.json"
            )
            golden = json.loads((working / filename).read_text())["H"][
                "schemes"
            ]["imgt"]
            parsed = PDBParser(QUIET=True).get_structure("output", output)
            residues = list(parsed[0]["F"])
            assert [(r.id[1], r.id[2].strip()) for r in residues] == [
                (n, code.strip()) for n, code, _ in golden["numbered"]
            ]
            assert [r.resname for r in residues] == [
                r.resname for r in source[0]["F"]
            ]
            left, right = list(source.get_atoms()), list(parsed.get_atoms())
            assert len(left) == len(right)
            for a, b in zip(left, right, strict=True):
                assert (
                    a.name,
                    a.altloc,
                    a.element,
                    a.occupancy,
                    a.bfactor,
                ) == (b.name, b.altloc, b.element, b.occupancy, b.bfactor)
                np.testing.assert_allclose(
                    a.coord, b.coord, rtol=0, atol=0.00051
                )
            print(f"Installed wheel {mode}: numbering and atoms verified.")
        # Exercise the packaged CCD mapping through the public API as well.
        first = next(source.get_residues())
        first.resname = "MSE"
        result = sabr.renumber_structure(source, "F", chain_type="H")
        assert next(result.get_residues()).resname == "MSE"
        assert next(source.get_residues()).id[1] == 502


if __name__ == "__main__":
    main()
