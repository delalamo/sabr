import logging
from pathlib import Path

import pytest
from Bio.PDB import PDBParser
from click.testing import CliRunner

import sabr.cli as cli

DATA = Path(__file__).parent / "data"


def _passthrough(captured=None, fail=False):
    def fake(structure, chain, **kwargs):
        logging.getLogger("sabr.test").info("pipeline details")
        if captured is not None:
            captured.update({"chain": chain, **kwargs})
        if fail:
            raise ValueError("planned failure")
        return structure

    return fake


def test_cli_maps_the_complete_compact_interface(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(cli, "renumber_structure", _passthrough(captured))
    output = tmp_path / "numbered.cif"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(output),
            "-n",
            "kabat",
            "-t",
            "L",
            "--noise-level",
            "0.5",
            "--residue-range",
            "2",
            "127",
            "--random-seed",
            "7",
            "--overwrite",
            "-v",
        ],
    )
    assert result.exit_code == 0, result.output
    assert output.exists()
    assert captured == {
        "chain": "F",
        "scheme": "kabat",
        "chain_type": "L",
        "noise_level": 0.5,
        "residue_range": (2, 127),
        "random_seed": 7,
    }
    assert "pipeline details" in result.output


def test_cli_defaults_are_deterministic_and_quiet(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(cli, "renumber_structure", _passthrough(captured))
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
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["scheme"] == "imgt"
    assert captured["chain_type"] == "auto"
    assert captured["noise_level"] == 0.0
    assert captured["random_seed"] == 0
    assert captured["residue_range"] is None
    assert "pipeline details" not in result.output


@pytest.mark.parametrize("noise_level", ["0.0", "0.2", "0.5", "1.0", "2.0"])
def test_cli_accepts_every_reference_noise_level(
    monkeypatch, tmp_path, noise_level
):
    captured = {}
    monkeypatch.setattr(cli, "renumber_structure", _passthrough(captured))
    output = tmp_path / f"noise-{noise_level}.pdb"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(output),
            "--noise-level",
            noise_level,
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["noise_level"] == float(noise_level)


def test_cli_overwrite_protection(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough())
    output = tmp_path / "exists.pdb"
    output.write_text("existing")
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(output),
        ],
    )
    assert result.exit_code != 0
    assert "--overwrite" in result.output
    assert output.read_text() == "existing"


def test_cli_failure_leaves_no_output_or_temporary_file(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough(fail=True))
    output = tmp_path / "failed.pdb"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(output),
        ],
    )
    assert result.exit_code != 0
    assert "planned failure" in result.output
    assert not output.exists()
    assert list(tmp_path.iterdir()) == []


def test_real_cli_round_trip(monkeypatch, tmp_path):
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
            "H",
        ],
    )
    assert result.exit_code == 0, result.output
    parsed = PDBParser(QUIET=True).get_structure("output", output)
    assert list(parsed[0]["F"])[-1].id[1] == 128


def test_pdb_output_rejects_extended_insertion_codes():
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    residue = list(structure[0]["F"])[0]
    residue.detach_parent()
    residue.id = (" ", residue.id[1], "AA")
    with pytest.raises(ValueError, match="extended insertion codes"):
        cli._validate_pdb_output(structure)


@pytest.mark.parametrize(
    "removed_option",
    [
        "--reference-chain-type",
        "--disable-custom-gap-penalties",
        "--disable-deterministic-renumbering",
    ],
)
def test_removed_compatibility_options_are_rejected(removed_option):
    result = CliRunner().invoke(cli.main, [removed_option])
    assert result.exit_code == 2
    assert "No such option" in result.output


def test_help_and_version_are_available():
    runner = CliRunner()
    help_result = runner.invoke(cli.main, ["--help"])
    version_result = runner.invoke(cli.main, ["--version"])
    assert help_result.exit_code == 0
    assert "--noise-level" in help_result.output
    assert version_result.exit_code == 0
