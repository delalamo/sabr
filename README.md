# SAbR

SAbR (Structure-based Antibody Renumbering) assigns antibody residue numbers
from backbone coordinates. It combines the original trained Haiku encoder with
the original affine Smith–Waterman alignment and ANARCI numbering rules.

SAbR is intentionally small and feature-complete. It provides one Python API
and one command-line program.

The complete usage guide is available in the
[SAbR documentation](https://sabr.readthedocs.io/).

## Installation

SAbR requires Python 3.11 or newer.

```bash
pip install sabr-kit
```

## Command line

```bash
sabr -i antibody.pdb -c H -o numbered.pdb
```

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

Defaults are IMGT numbering, automatic H/K/L selection, noise level `0.0`,
`sabr` mode, and the entire selected chain. Existing outputs are never
replaced unless `--overwrite` is given. Normal output contains only warnings
and errors; `-v` reports reference scores and pipeline decisions.

Use `--mode softalign` to select the original SoftAlign encoder weights,
reference embeddings, and affine gap penalties together. SoftAlign references
do not vary with `--noise-level`, so that option is ignored in this mode.

Input and output may be PDB (`.pdb`) or mmCIF (`.cif` or `.mmcif`). Use mmCIF
when chain names or ANARCI insertion codes exceed PDB's one-character fields.
Writes are atomic, so a failed run does not leave a partial output.

CLI conversion guarantees preservation of atomic structure content, not
arbitrary non-atomic mmCIF categories. It warns for every mmCIF input. When
those categories matter, load a Gemmi structure and use the in-memory API.

## Python API

```python
from Bio.PDB import PDBParser
from sabr import renumber_structure

structure = PDBParser(QUIET=True).get_structure("antibody", "antibody.pdb")
numbered = renumber_structure(structure, chain="H")
```

Gemmi structures use the same function:

```python
import gemmi
from sabr import renumber_structure

structure = gemmi.read_structure("antibody.cif")
numbered = renumber_structure(
    structure,
    chain="heavy_chain",
    scheme="chothia",
    chain_type="auto",
    noise_level=0.0,
    mode="softalign",
    residue_range=None,
)
```

`renumber_structure` never mutates its input and returns the same concrete
structure type. Non-target chains, hetero residues, waters, and residues
outside an inclusive `residue_range` are preserved. SAbR rejects multi-model
structures rather than silently modifying only one model. If a partial range
would create duplicate residue IDs with unchanged residues, the operation
fails with an explanation.

The same-type clone preserves metadata represented by the input BioPython or
Gemmi object. Alternate conformers are normalized deterministically: a
complete blank-altloc backbone is preferred, then the complete conformer with
the greatest summed occupancy, with altloc name as the final tie-breaker.
Selections above 1,024 polymer residues are rejected before quadratic model
work; use `residue_range` to select the antibody domain.

Modified peptide residues are translated only for sequence generation. Their
original names and atoms remain unchanged. The committed mapping was generated
from the wwPDB Chemical Component Dictionary snapshot dated 2026-07-11
(`components.cif.gz` SHA-256
`0b3323123ec10b997afe1c530b4cad30306e60b451b2b062c59bc9bb5cbe0679`) and
contains only peptide-linking components with exactly one canonical amino-acid
parent. Unsupported or ambiguous polymer chemistry fails explicitly; no
runtime network access occurs.

Gemmi itself only represents one-character insertion codes. For unusually
long loops that need extended codes, use a BioPython structure in memory or
the CLI with mmCIF output.

## Scientific behavior

- The default `sabr` mode preserves the trained SAbR encoder weights,
  references, and gap penalties unchanged.
- The optional `softalign` mode uses `softalign_encoder.npz`,
  `softalign_embeddings.npz`, and the exact penalties in
  `softalign_gap.npz` as one parameter set.
- Alignment uses the original differentiable affine Smith–Waterman method.
- In `sabr` mode, gap extension is `-0.175027` and gap opening is `-2.525591`.
  In `softalign` mode, they are `0.1942468136548996` and
  `-2.5441808700561523`, respectively, as stored in the repository asset.
- CDR gap distribution is always applied.
- No deterministic light-chain DE-loop or C-terminal correction is applied.
- Automatic chain selection aligns against H, K, and L references and uses
  the highest score, with deterministic H/K/L tie order.

A structural gap is detected when the C–N distance between consecutive
residues exceeds 2.66 Å. A gap skips only the affected CDR correction and
emits a warning; other regions continue normally.

## Development

```bash
pip install -c constraints.txt -e '.[test]'
JAX_PLATFORMS=cpu pytest
pre-commit run --all-files
```

`constraints.txt` records the exact canonical development and CI environment.
Package metadata remains ranged for normal installation. SAbR does not force a
JAX backend; CPU is simply the canonical CI regression baseline.

The committed tests are self-contained and never download data. They verify
the fixed asset hashes, encoder and alignment baselines, all numbering schemes,
H/K/L selection, regional corrections, structure-object behavior, and CLI
failure handling.

## Deferred full benchmark

The historical pre-2021 SAbDab manifest contains approximately 1,012 chains.
The current method scores about 90% on that set, not 100%. The full corpus is
not bundled or downloaded by CI.

Future benchmark work should create a checksum-pinned corpus, verify residue
IDs, insertion codes, and coordinate parity, and compare a lossless archive
with Foldcomp before adding a separate manual or nightly workflow. This is a
benchmarking TODO, not a unit-test or release requirement.

## License and attribution

SAbR is distributed under the repository license. The vendored ANARCI
numbering code retains its original license in `src/sabr/_anarci/LICENSE`.
