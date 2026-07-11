"""The SAbR command-line interface."""

import logging
import os
import uuid
from pathlib import Path

import click
from Bio.PDB import MMCIFIO, PDBIO, MMCIFParser, PDBParser

from sabr import constants, renumber_structure


def _read_structure(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        return PDBParser(QUIET=True).get_structure(path.stem, path)
    if suffix in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True).get_structure(path.stem, path)
    raise ValueError("Input must use the .pdb, .cif, or .mmcif extension.")


def _validate_pdb_output(structure) -> None:
    for chain in structure.get_chains():
        if len(chain.id) != 1:
            raise ValueError(
                f"PDB output requires one-character chain IDs; '{chain.id}' "
                "requires mmCIF output."
            )
        for residue in chain:
            if len(residue.id[2].strip()) > 1:
                raise ValueError(
                    "PDB output cannot represent extended insertion codes; "
                    "use .cif or .mmcif output."
                )


def _write_structure(structure, path: Path) -> None:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        _validate_pdb_output(structure)
        writer = PDBIO()
    elif suffix in (".cif", ".mmcif"):
        writer = MMCIFIO()
    else:
        raise ValueError("Output must use the .pdb, .cif, or .mmcif extension.")
    writer.set_structure(structure)
    writer.save(str(path))


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(package_name="sabr-kit")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Input PDB or mmCIF file.",
)
@click.option(
    "-c", "--chain", required=True, help="Chain identifier to renumber."
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output PDB or mmCIF file.",
)
@click.option(
    "-n",
    "--scheme",
    default="imgt",
    show_default=True,
    type=click.Choice(constants.NUMBERING_SCHEMES, case_sensitive=False),
)
@click.option(
    "-t",
    "--chain-type",
    default="auto",
    show_default=True,
    type=click.Choice(("auto", *constants.CHAIN_TYPES), case_sensitive=True),
    help="Reference and ANARCI chain type.",
)
@click.option(
    "--noise-level",
    default="0.0",
    show_default=True,
    type=click.Choice(
        tuple(f"{level:.1f}" for level in constants.NOISE_LEVELS)
    ),
)
@click.option("--residue-range", nargs=2, type=int, metavar="START END")
@click.option("--overwrite", is_flag=True)
@click.option("-v", "--verbose", is_flag=True)
def main(
    input_path: Path,
    chain: str,
    output_path: Path,
    scheme: str,
    chain_type: str,
    noise_level: str,
    residue_range: tuple | None,
    overwrite: bool,
    verbose: bool,
) -> None:
    """Renumber one antibody chain in a PDB or mmCIF structure."""
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
        force=True,
    )
    if output_path.exists() and not overwrite:
        raise click.ClickException(
            f"Output already exists: {output_path}. "
            "Use --overwrite to replace it."
        )

    temporary = output_path.with_name(
        f".{output_path.stem}.{uuid.uuid4().hex}.tmp{output_path.suffix}"
    )
    try:
        structure = _read_structure(input_path)
        result = renumber_structure(
            structure,
            chain,
            scheme=scheme,
            chain_type=chain_type,
            noise_level=float(noise_level),
            residue_range=residue_range,
        )
        _write_structure(result, temporary)
        os.replace(temporary, output_path)
    except click.ClickException:
        raise
    except Exception as error:
        if temporary.exists():
            temporary.unlink()
        raise click.ClickException(str(error)) from error


if __name__ == "__main__":
    main()
