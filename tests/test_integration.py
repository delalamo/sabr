import json
from pathlib import Path

import pytest
from Bio.PDB import PDBIO, MMCIFParser, PDBParser
from click.testing import CliRunner

import sabr.cli as cli
from tests.helpers import assert_atomic_content_equal, residue_ids

DATA = Path(__file__).parent / "data"
pytestmark = pytest.mark.integration


@pytest.fixture
def offset_heavy_input(tmp_path):
    source = PDBParser(QUIET=True).get_structure(
        "input", DATA / "test_heavy_chain.pdb"
    )
    for residue in source[0]["F"]:
        residue.id = (residue.id[0], residue.id[1] + 500, residue.id[2])
    path = tmp_path / "offset-input.pdb"
    writer = PDBIO()
    writer.set_structure(source)
    writer.save(str(path))
    return path


def test_real_cli_round_trip(tmp_path, offset_heavy_input):
    output = tmp_path / "numbered.pdb"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(offset_heavy_input),
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
    golden = json.loads((DATA / "numbering_baseline.json").read_text())["H"][
        "schemes"
    ]["imgt"]
    assert golden["first_row"] == 0
    assert residue_ids(parsed, "F") == [
        (n, code.strip()) for n, code, _ in golden["numbered"]
    ]
    source = PDBParser(QUIET=True).get_structure("input", offset_heavy_input)
    assert residue_ids(parsed, "F") != residue_ids(source, "F")
    assert_atomic_content_equal(source, parsed, serialized=True)


def test_real_softalign_cli_round_trip(tmp_path, offset_heavy_input):
    output = tmp_path / "softalign-numbered.pdb"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(offset_heavy_input),
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
    golden = json.loads(
        (DATA / "softalign_numbering_baseline.json").read_text()
    )["H"]["schemes"]["imgt"]
    assert golden["first_row"] == 0
    assert residue_ids(parsed, "F") == [
        (n, code.strip()) for n, code, _ in golden["numbered"]
    ]
    source = PDBParser(QUIET=True).get_structure("input", offset_heavy_input)
    assert residue_ids(parsed, "F") != residue_ids(source, "F")
    assert_atomic_content_equal(source, parsed, serialized=True)


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
