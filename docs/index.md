# SAbR

SAbR assigns antibody residue numbers from backbone coordinates. It accepts
PDB and mmCIF structures, preserves the input object, and supports IMGT,
Chothia, Kabat, Martin, AHo, and Wolfguy numbering.

## Installation

SAbR requires Python 3.11 or newer.

```bash
pip install sabr-kit
```

## Command line

Renumber chain `H` using the default IMGT scheme and automatic heavy,
kappa, or lambda reference selection:

```bash
sabr -i antibody.pdb -c H -o numbered.pdb
```

Select a numbering scheme and chain type explicitly:

```bash
sabr -i antibody.cif -c light_chain -o numbered.cif \
  --scheme chothia --chain-type K
```

Use the complete SoftAlign parameter set (encoder, reference embeddings, and
gap penalties):

```bash
sabr -i antibody.pdb -c H -o numbered.pdb --mode softalign
```

Renumber only residues whose source PDB numbers are 1 through 130:

```bash
sabr -i antibody.pdb -c H -o numbered.pdb \
  --residue-range 1 130
```

Insertion-coded residues are included when their numeric component falls in
the range. Existing files are protected unless `--overwrite` is supplied.

The complete interface is:

```text
sabr -i INPUT -c CHAIN -o OUTPUT
     [-n imgt|chothia|kabat|martin|aho|wolfguy]
     [-t auto|H|K|L]
     [--noise-level 0.0|0.2|0.5|1.0|2.0]
     [-m sabr|softalign]
     [--residue-range START END]
     [--overwrite] [-v]
```

The defaults are IMGT, automatic chain selection, noise level `0.0`, and
`sabr` mode. SoftAlign mode uses its own fixed references, so `noise_level` is
ignored in that mode. Normal output contains only warnings and errors. Use
`--verbose` to show the JAX backend, chain-selection scores, and a traceback
on failure.

## Python API

`renumber_structure` accepts either a BioPython or Gemmi `Structure`. It
returns the same concrete type and does not mutate its input.

### BioPython

```python
from Bio.PDB import PDBParser
from sabr import renumber_structure

structure = PDBParser(QUIET=True).get_structure("antibody", "antibody.pdb")
numbered = renumber_structure(
    structure,
    chain="H",
    scheme="imgt",
    chain_type="auto",
)
```

### Gemmi

```python
import gemmi
from sabr import renumber_structure

structure = gemmi.read_structure("antibody.cif")
numbered = renumber_structure(
    structure,
    chain="heavy_chain",
    scheme="kabat",
    chain_type="H",
    noise_level=0.0,
    mode="softalign",
    residue_range=(1, 130),
)
```

The function signature is:

```python
renumber_structure(
    structure,
    chain: str,
    scheme: str = "imgt",
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
    mode: str = "sabr",
)
```

`mode="softalign"` selects the SoftAlign encoder weights, reference
embeddings, and exact gap penalties stored in `softalign_gap.npz`. The default
`mode="sabr"` preserves existing behavior.

Non-target chains, waters, ligands, metadata represented by the input object,
and residues outside the selected range are preserved in the returned clone.
Only single-model structures are supported.

## PDB and mmCIF output

Use mmCIF when a structure has:

- a chain ID longer than one character;
- multi-character insertion codes;
- residue numbers outside the PDB range `-999` through `9999`; or
- more than 99,999 atoms.

Gemmi structure objects can represent only one-character insertion codes. For
exceptionally long CDR insertions, use a BioPython structure and mmCIF output.

CLI conversion preserves atomic structure content but may not preserve every
non-atomic mmCIF category. Use the in-memory Gemmi API when those categories
are required.

## Structural gaps and modified residues

A C–N distance above 2.66 Å is treated as a structural gap. If a gap crosses a
CDR, SAbR warns and retains the learned alignment for that region while
continuing corrections elsewhere.

Supported modified peptide residues are translated to their canonical parent
only for sequence generation. Original residue names and atoms are preserved.
Unknown or ambiguous polymer chemistry is rejected rather than converted to
`X`.

## Errors

Common errors include:

- selecting a missing chain;
- selecting more than one structural model;
- residues missing N, CA, or C atoms;
- a range containing multiple or discontinuous domains;
- output that cannot be represented in PDB format; and
- selections above the 1,024-residue safety limit.

For long or multi-domain chains, pass `--residue-range` on the command line or
`residue_range=(start, end)` in Python.
