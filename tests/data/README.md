# Regression data

Existing `math_baseline.npz`, `numbering_baseline.json`, and
`8sve_cdr1_imgt.json` are historical accepted regression baselines. This change
preserves their bytes. Their original capture environment is not reconstructed
or invented here.

The additional fixtures and captures are described in
`regression_provenance.json`: source URLs, author chain IDs, SHA-256 hashes of
original downloads and extracted atomic mmCIF files, input and output hashes,
producing commit, Python/backend, dependencies, and capture parameters.
8DY0 chain A has 244 polymer residues and no structural gaps; 8DY1 chain C has
228 polymer residues and a gap after query row 110. Extraction retains all
atomic content in the selected chain, including alternate conformers and
hetero residues. No coordinates or missing residues are fabricated.

The new goldens were captured from the numerical implementation at commit
2227a34, before production changes. Review checked complete query-row and
sequence correspondence, residue counts, unique output IDs, first/last
positions, domain boundaries, and differences between parameter modes.
They preserve observed behavior; they are not independent biological accuracy
labels. Synthetic tests compare domain composition against the existing
single-domain goldens and explicitly do not model protein geometry.

## Known mode difference

On 8DY0, both modes select HK and start the K domain at query row 137, numbered
1001. SAbR numbers query rows 119–136 as sequential linker residues 129–146.
SoftAlign instead places query rows 118–134 at heavy-domain insertions
127A–127Q, row 135 at 128, and only row 136 as the sequential linker residue
129. Consequently 19 residue mappings differ. This is an existing alignment
behavior, not evidence that the synthetic linker or annotation is wrong.
The regression intentionally preserves it; interpreting the entire annotated
linker as unassigned would require a separately reviewed scientific change.
8DY1's SAbR and SoftAlign mappings agree; neither reconstructs its missing
linker coordinates.

## Capturing candidates manually

From a checkout in the pinned CPU environment:

```sh
JAX_PLATFORMS=cpu python scripts/capture_regressions.py \
  --repo "$PWD" --output-dir /tmp/sabr-new-candidates
```

The output directory must not exist. This explicitly invoked script downloads
the two public structures and writes candidates separately; pytest never
imports it, downloads structures, or regenerates expected outputs. Review
numerical differences and provenance before copying accepted files to this
directory. Do not refresh existing baselines merely to make a failing test
pass. Record the review in the manifest and update output checksums only for
intentionally accepted new captures.

Numerical regression tolerances remain those of the historical baseline
(rtol 1e-5 and array-specific absolute tolerances). Rigid-transform checks use
rtol 2e-4 / atol 2e-5 for accumulated floating-point geometry, while requiring
identical discrete numbering. In-memory preservation is exact; serialization
allows 0.00051 Å coordinate and 0.0051 occupancy/B-factor rounding, with no
relative tolerance. Writer-normalized atom serials and full names are not
asserted across file formats.
