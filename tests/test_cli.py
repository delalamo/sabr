import json
import logging
import shutil
from pathlib import Path

import pytest
from Bio.PDB import MMCIFParser, PDBParser
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
        "scfv": False,
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
    assert captured["residue_range"] is None
    assert captured["scfv"] is False
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


def test_cli_forwards_scfv_mode(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(cli, "renumber_structure", _passthrough(captured))
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(tmp_path / "scfv.pdb"),
            "--scfv",
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["chain_type"] == "auto"
    assert captured["scfv"] is True


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
        ],
    )
    assert result.exit_code == 0, result.output
    parsed = MMCIFParser(QUIET=True).get_structure("output", output)
    actual = [
        [residue.id[1], residue.id[2].strip()]
        for residue in parsed[0]["M"]
        if not residue.id[0].strip()
    ]
    assert actual == [[item[2], item[3]] for item in golden["mapping"]]


def test_cli_accepts_mmcif_input(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough())
    output = tmp_path / "numbered.cif"
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_minimal.cif"),
            "-c",
            "A",
            "-o",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    assert result.output == (
        "WARNING: Non-atomic mmCIF categories may not be preserved by CLI "
        "conversion; use the in-memory Gemmi API when metadata preservation "
        "is required.\n"
    )
    parsed = MMCIFParser(QUIET=True).get_structure("output", output)
    assert parsed[0]["A"]


@pytest.mark.parametrize("invalid_side", ["input", "output"])
def test_cli_rejects_unsupported_file_extensions(
    monkeypatch, tmp_path, invalid_side
):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough())
    input_path = DATA / "test_heavy_chain.pdb"
    output_path = tmp_path / "numbered.pdb"
    if invalid_side == "input":
        input_path = tmp_path / "structure.txt"
        shutil.copyfile(DATA / "test_heavy_chain.pdb", input_path)
    else:
        output_path = tmp_path / "numbered.txt"
    result = CliRunner().invoke(
        cli.main,
        ["-i", str(input_path), "-c", "F", "-o", str(output_path)],
    )
    assert result.exit_code != 0
    assert ".pdb, .cif, or .mmcif" in result.output
    assert not output_path.exists()


def test_writer_failure_removes_partial_temporary_file(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough())

    def fail_after_partial_write(structure, path):
        path.write_text("partial")
        raise OSError("planned writer failure")

    monkeypatch.setattr(cli, "_write_structure", fail_after_partial_write)
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
    assert "planned writer failure" in result.output
    assert list(tmp_path.iterdir()) == []


def test_cleanup_failure_does_not_replace_original_error(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough())

    def fail_after_partial_write(structure, path):
        path.write_text("partial")
        raise OSError("original writer failure")

    original_unlink = Path.unlink

    def fail_temporary_cleanup(path, *args, **kwargs):
        if path.name.startswith(".failed."):
            raise OSError("cleanup failure")
        return original_unlink(path, *args, **kwargs)

    monkeypatch.setattr(cli, "_write_structure", fail_after_partial_write)
    monkeypatch.setattr(Path, "unlink", fail_temporary_cleanup)
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(tmp_path / "failed.pdb"),
        ],
    )
    assert result.exit_code != 0
    assert "original writer failure" in result.output
    assert "cleanup failure" in result.output


def test_pdb_output_rejects_extended_insertion_codes():
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    residue = list(structure[0]["F"])[0]
    residue.detach_parent()
    residue.id = (" ", residue.id[1], "AA")
    with pytest.raises(ValueError, match="printable ASCII insertion"):
        cli._validate_pdb_output(structure)


def test_verbose_mode_reports_backend_and_traceback(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "renumber_structure", _passthrough(fail=True))
    result = CliRunner().invoke(
        cli.main,
        [
            "-i",
            str(DATA / "test_heavy_chain.pdb"),
            "-c",
            "F",
            "-o",
            str(tmp_path / "failed.pdb"),
            "--verbose",
        ],
    )
    assert result.exit_code != 0
    assert "JAX backend:" in result.output
    assert "Traceback (most recent call last)" in result.output


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda chain, residue: setattr(chain, "id", "AB"), "chain character"),
        (
            lambda chain, residue: setattr(residue, "id", (" ", 10000, " ")),
            "residue number",
        ),
        (
            lambda chain, residue: setattr(residue, "id", (" ", -1000, " ")),
            "residue number",
        ),
        (
            lambda chain, residue: setattr(residue, "id", (" ", 1, "\x01")),
            "insertion character",
        ),
    ],
)
def test_every_pdb_field_limit_is_validated(mutation, message):
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    chain = structure[0]["F"]
    residue = list(chain)[0]
    mutation(chain, residue)
    with pytest.raises(ValueError, match=message):
        cli._validate_pdb_output(structure)


def test_pdb_atom_count_limit_is_validated(monkeypatch):
    structure = PDBParser(QUIET=True).get_structure(
        "heavy", DATA / "test_heavy_chain.pdb"
    )
    atom = next(structure.get_atoms())
    monkeypatch.setattr(
        structure,
        "get_atoms",
        lambda: (atom for _ in range(100_000)),
    )
    with pytest.raises(ValueError, match="99,999 atoms"):
        cli._validate_pdb_output(structure)


def test_help_and_version_are_available():
    runner = CliRunner()
    help_result = runner.invoke(cli.main, ["--help"])
    version_result = runner.invoke(cli.main, ["--version"])
    assert help_result.exit_code == 0
    assert "--noise-level" in help_result.output
    assert "--scfv" in help_result.output
    assert version_result.exit_code == 0
