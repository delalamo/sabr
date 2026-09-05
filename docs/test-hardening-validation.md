# Test hardening validation

Implementation branch: `codex/test-suite-hardening`, based on `2227a34`.
Worktree: `/Users/diegoda/Documents/SAbR-test-hardening`.
Validation date: 2026-09-05. Platform: macOS ARM64, JAX CPU.

## Results

| Environment | Result | Runtime | Statement coverage | Branch coverage | Combined coverage |
| --- | --- | --- | --- | --- | --- |
| Python 3.11.16, pinned constraints | 355 passed | 71.22 s | 98.70% | 95.03% | 97.81% |
| Python 3.12.14, pinned constraints | 355 passed | 67.83 s | 98.70% | 95.03% | 97.81% |

Both environments were independently installed from `constraints.txt`. These
are local CPU results; the updated GitHub workflows have not been dispatched.
The final full runs used the optimized fixture implementation at `f4a816a`.
No test failures remain. Both runs report two upstream Haiku/JAX deprecation
warnings and the existing expected 8sve structural-gap warning.

Validation command in each environment:

```sh
JAX_PLATFORMS=cpu python -m pytest --cov=sabr \
  --cov-report=term-missing --cov-report=json --durations=10
```

Coverage measures branches through repository configuration, retains the 95%
aggregate gate, and introduces no new coverage exclusions. The numerical
model, reference assets, and the three original golden files are unchanged.
A git comparison against `2227a34` verified those files byte-for-byte.

## Behavioral defects verified before fixing

New regressions initially reproduced temporary-file leakage when a writer
raised `click.ClickException`, late validation of unsupported output extensions
and TCR schemes, and unnecessary JAX backend initialization in quiet mode.
Non-finite coordinates also produced NumPy arithmetic errors before the
intended validation error when callers enabled strict floating-point handling.
Targeted tests passed after the small production fixes; no scientific
algorithm or parameter changes were needed.

The heavy-chain input was already correctly numbered, so both CLI and wheel
checks now offset input IDs by 500 before requiring the complete expected
mapping. This prevents a passthrough implementation from satisfying even a
full mapping assertion.

## Mutation, lint, and packaging checks

All four subprocess-only mutations were detected by test failures:

- Bypass renumbering entirely.
- Reverse the equal-occupancy altloc name tie-break.
- Remove multi-domain numbering offsets.
- Suppress CDR2 correction while retaining the other regional corrections.

The mutation runner never edits repository sources. Run it with
`JAX_PLATFORMS=cpu python scripts/check_test_mutations.py`.

All pre-commit hooks passed, including the relocated architecture check.
Workflow YAML parsed successfully. Source-distribution and wheel builds passed.
The non-editable wheel was installed into a clean pinned environment; both
CLI modes restored full numbering from offset input, preserved atomic content,
and the public API exercised the packaged modified-residue mapping. CLI
subprocesses ran outside the checkout with `PYTHONPATH` cleared. Test and
publish workflows share this wheel check.

## Timing and remaining scientific limits

The slowest final case was the complete SAbR 8DY1 pipeline/CLI round trip
(9.92 s on both Python versions), followed by the SAbR 8DY0 case (7.24–7.33 s)
and first-use K/L noise-reference cases (6.39–6.64 s, including cached fixture
preparation). Timings include compilation and concurrent validation workloads;
they are observations, not performance thresholds.

An initial timing pass exposed repeated deep copies of Biopython parent
structures during synthetic fixture composition. Detached residue copies
removed that cost without weakening the assertions. Both full suites were
rerun after this change.

New captures preserve observed historical output, not independent biological
accuracy labels. SoftAlign assigns part of the annotated 8DY0 linker to
heavy-domain insertions; this pre-existing difference is retained explicitly.
See [fixture provenance and review](../tests/data/README.md) for the exact
query rows, capture procedure, checksums, and numerical tolerances. The full
historical SAbDab benchmark remains deferred.

No changes were merged or published, and the original worktree and performance
branch were not modified by this work.
