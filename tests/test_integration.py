import json
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser
from click.testing import CliRunner

import sabr.cli as cli

DATA = Path(__file__).parent / "data"


def test_real_cli_round_trip(tmp_path):
    output = tmp_path / "numbered.pdb"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(output),
            "-t",
            "H,K",
        ],
    )
    assert result.exit_code == 0, result.output
    parsed = PDBParser(QUIET=True).get_structure("output", output)
    assert list(parsed[0]["F"])[-1].id[1] == 128


def test_real_softalign_cli_round_trip(tmp_path):
    output = tmp_path / "softalign-numbered.pdb"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(output),
            "-t",
            "H",
            "--mode",
            "softalign",
        ],
    )
    assert result.exit_code == 0, result.output
    parsed = PDBParser(QUIET=True).get_structure("output", output)
    assert list(parsed[0]["F"])[-1].id[1] == 128


def test_8sve_cli_mmcif_round_trip_matches_full_golden(tmp_path):
    golden = json.loads((DATA / "8sve_cdr1_imgt.json").read_text())
    output = tmp_path / "8sve-numbered.cif"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "8sve_L.pdb"),
            "-c",
            "M",
            "-o",
            str(output),
            "-t",
            "K",
            "--dangerously-allow-structural-gaps",
        ],
    )
    assert result.exit_code == 0, result.output
    assert result.output.startswith(
        "WARNING: Structural-gap safety check disabled by "
        "--dangerously-allow-structural-gaps"
    )
    parsed = MMCIFParser(QUIET=True).get_structure("output", output)
    actual = [
        [residue.id[1], residue.id[2].strip()]
        for residue in parsed[0]["M"]
        if not residue.id[0].strip()
    ]
    assert actual == [[item[2], item[3]] for item in golden["mapping"]]
