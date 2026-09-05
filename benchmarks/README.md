# Inference performance experiments

These experiments use the four bundled structures, without downloading a
benchmark corpus. Run every stage sequentially on the same machine. The default
benchmark backend is CPU, matching the repository's numerical regression baseline.

## Reproduce

```sh
python benchmarks/inference.py --output /tmp/sabr-baseline
# Apply one optimization, then check against the saved outputs:
python benchmarks/inference.py --output /tmp/sabr-after --reference /tmp/sabr-baseline
# Alternate the previous revision and working-tree implementations:
python benchmarks/paired.py --before PREVIOUS_REVISION --output /tmp/sabr-paired
```

`inference.py` starts a fresh process for every structure and parameter mode.
It records first-use pipeline latency (including JAX initialization, model asset
loading and compilation), seven warm stage samples, and seven warm public-API
samples. Parsing and Python imports are outside the pipeline timer; import time
is recorded separately. JAX results are synchronized before timing stops.
It checks embeddings and scores against the original numerical tolerances,
and requires identical corrected alignments, chain selections and numbering.
The long-loop fixture explicitly opts into its known structural gaps.

`paired.py` alternates the previous and new implementations for eleven pairs
after warming both. This limits bias from changing desktop load. It loads only
`model.py` and `alignment.py` from the requested Git revision; it is intended for
these experiments, where extraction, correction and numbering are unchanged.
Its pipeline timing excludes extraction and includes numbering and structure copy.

## Reusing compilation

For repeated CLI processes, JAX supports an opt-in persistent cache:

```sh
export JAX_COMPILATION_CACHE_DIR="$HOME/.cache/sabr/jax"
export JAX_PERSISTENT_CACHE_MIN_COMPILE_TIME_SECS=0
```

The second setting includes compilations below JAX's normal one-second cache
threshold. Use a private, writable cache directory. The library leaves these
settings under the caller's control. A persistent Python process also reuses
in-memory compilations and cached model assets. Cache reuse requires compatible
input shapes, dtypes, JAX versions and devices; a new residue count can still
require compilation. The first population of a cache does not avoid compilation.

The benchmark disables persistent caching by default. To measure it, run
`inference.py` twice with the same `--cache-dir /tmp/sabr-jax-cache`, saving each
run to a separate output directory.

## Measurement environment

Apple M2 MacBook Air, 8 GB RAM, macOS, Python 3.12.8; dependencies installed from
`constraints.txt` (JAX/jaxlib 0.10.2, NumPy 2.5.1). The desktop had competing CPU
and memory activity. Absolute first-use and between-process timings varied
substantially. Small improvements should be treated as inconclusive; paired
measurements are more useful for judging warm performance.

Results and the final retention decisions are recorded below after each stage.

## Results by change

Changes were implemented, checked and measured individually. A speedup above
1 means faster; small wall-time differences on this desktop are inconclusive.

| Step | Comparison | Observation | Decision |
| --- | --- | --- | --- |
| 1. Whole-encoder JIT | `f629a0c` → `098f242` | Paired warm pipeline speedups: H 1.04×, K 0.96×, L 1.20×, long loop 1.11×. First-use latency generally improved. | Keep, primarily for compilation reuse; no large warm speedup claim. |
| 2. Persistent cache | Same JIT code, cache population vs reuse | Cache files were generated. First-use timing improved for H and long loop but regressed for K/L amid changing load. | Document opt-in configuration; no reliable wall-time speedup claim from this run. |
| 3. Neighbor-only distances | JIT encoder vs gather-before-distance | Inconsistent performance. One of 7,808 baseline embedding elements exceeded the existing tolerance (absolute difference 2.1644e-6). 91 tests passed, one failed. | Revert; patch retained in `experiments/neighbor-distances.patch`. |
| 4. Device embeddings | `613e056` → `184bd6d` | Essentially neutral CPU cost. Preserve the NumPy helper while avoiding the intermediate host conversion in the API. | Keep; accelerator benefit is unmeasured. |
| 5. Reference batching | `184bd6d` → `512c86f` | H/K/L alignment wall-time speedups 1.22× / 1.27× / 1.92×; CPU-time speedups 1.44× / 1.35× / 1.65×. Overall CPU savings approximately 4–8%. | Keep. |
| 6. Score-only selection | Batched full alignment vs scores then winner gradient | H/K/L pipeline wall-time speedups 0.99× / 1.03× / 1.00×; CPU changes approximately 0–2%. 100 relevant tests passed. | Revert; extra selection logic does not justify the marginal overall benefit. Patch retained in `experiments/score-first.patch`. |
| 7. Encoder length buckets | `2f2291d` → `25f057a` | Ten nearby lengths need one encoder executable instead of ten. Warm pipeline CPU cost increased by up to approximately 6% on the fixtures. | Keep as a startup/variable-length throughput tradeoff, isolated in its own commit. |

The long-loop fixture explicitly selects K, so reference batching and score-only
selection do not change its alignment path. Its wall-time fluctuations are a
useful control demonstrating desktop measurement noise.

### Length reuse

These are **encoder-only microbenchmarks** on progressively cropped fixture
coordinates, not ten independent complete antibody systems. The two paths run
in one process, alternating which path goes first at successive lengths. The
first baseline call includes backend initialization; executable counts provide
independent evidence of the reduction in compilation. All embeddings passed the
same numerical tolerances at each length.

| Fixture crops | Residue counts | Exact-shape total | Bucketed total | Speedup | Compiled variants |
| --- | --- | ---: | ---: | ---: | --- |
| Heavy | 113–122 | 30.03 s | 4.84 s | 6.21× | 10 → 1 |
| Kappa | 205–214 | 5.75 s | 1.27 s | 4.52× | 10 → 1 |

Reproduce with:

```sh
python benchmarks/paired.py --before 2f2291d --output /tmp/sabr-lengths --case heavy --length-sweep
python benchmarks/paired.py --before 2f2291d --output /tmp/sabr-lengths --case kappa --length-sweep
```

For inputs above 64 residues, the encoder pads to a multiple of 32 and trims its
output to the actual length. Padded residues are excluded from the neighbor
search with infinite distance. Smaller graphs retain their exact size to
preserve the original neighbor count. Alignment uses actual query lengths;
only reference lengths are padded for single-domain candidate batching.
Multi-domain candidates retain their existing scalar alignment and terminal-gap
scoring path.

A steady workload with one fixed, already-compiled length can pay the padding
cost without gaining compilation reuse. The separate bucketing commit makes
that tradeoff reviewable. No model weights, learned feature ordering, gap
parameters, precision defaults, or historical encoder argument ordering changed.

### Saved evidence

`results/00-baseline` contains original embeddings, alignments and scores plus
timing samples. Subsequent result directories contain timings and parity status;
numbering is compared through a SHA-256 digest of its canonical JSON representation.
`paired_*.json` holds all eleven alternating warm samples and medians. Starting
with step 4, it also includes process CPU time, summed across process threads,
which helps distinguish compute changes from scheduling delays. CPU time and
wall time are different quantities and are not interchangeable.

The failed neighbor experiment's test log is preserved. Its error refers to the
original, unchanged test tolerances; no baseline or tolerance was relaxed.

## Final comparison with the original code

The final pipeline is compared against `f629a0c` (unchanged inference code from
`2227a34`). Each row is the median of eleven alternating, synchronized warm
runs. Pipeline time includes encoding, alignment, numbering and structure copy;
it excludes extraction. Both parameter modes pass embedding/score tolerances
and exact alignment, selection and numbering parity against the original.

| Fixture | Mode | Original pipeline | Final pipeline | Wall speedup | CPU time reduction |
| --- | --- | ---: | ---: | ---: | ---: |
| heavy (122 residues) | sabr | 68.3 ms | 63.0 ms | 1.08× | 20.8% |
| heavy (122 residues) | softalign | 67.8 ms | 63.4 ms | 1.07× | 19.3% |
| kappa (214 residues) | sabr | 99.9 ms | 96.2 ms | 1.04× | 16.8% |
| kappa (214 residues) | softalign | 99.8 ms | 96.3 ms | 1.04× | 16.8% |
| lambda (215 residues) | sabr | 106.8 ms | 103.0 ms | 1.04× | 19.2% |
| lambda (215 residues) | softalign | 100.4 ms | 96.4 ms | 1.04× | 17.4% |
| long_loop (250 residues) | sabr | 111.6 ms | 106.9 ms | 1.04× | 18.4% |
| long_loop (250 residues) | softalign | 111.5 ms | 107.5 ms | 1.04× | 18.9% |

The results support a modest warm wall-time improvement (about 4–8%), with a
larger reduction in CPU work (about 17–21%). The 4.5–6.2× length-sweep results
above apply to encoder compilation reuse and must not be read as whole-pipeline
steady-state speedups. No GPU measurements or full historical corpus were run.

### First use in a fresh process

The original implementation was rechecked immediately after the final one
because the initial baseline run had severe desktop-load variation. Source
hashes in `09-original-recheck` confirm that it ran the original model and
alignment code from `2227a34`. Every run disables persistent caching and starts
a fresh Python process. Parsing/import time is excluded; initialization, asset
loading, compilation and the first pipeline execution are included. These are
single first-use samples, so the range is more useful than small differences.

| Fixture | Mode | Original recheck | Final | Speedup |
| --- | --- | ---: | ---: | ---: |
| heavy | sabr | 3.863 s | 1.438 s | 2.69× |
| heavy | softalign | 3.575 s | 1.123 s | 3.18× |
| kappa | sabr | 3.724 s | 1.206 s | 3.09× |
| kappa | softalign | 3.610 s | 1.178 s | 3.07× |
| lambda | sabr | 3.544 s | 1.126 s | 3.15× |
| lambda | softalign | 3.848 s | 1.140 s | 3.37× |
| long_loop | sabr | 2.845 s | 1.184 s | 2.40× |
| long_loop | softalign | 2.861 s | 1.210 s | 2.36× |

## Validation

- Original code: **199 tests passed**, 96.44% coverage.
- Final code: **227 tests passed**, **96.51% coverage**, using
  `JAX_PLATFORMS=cpu JAX_ENABLE_COMPILATION_CACHE=false pytest --cov=sabr`.
- All eight full-structure/parameter-mode comparisons matched the original
  embeddings and scores within the original tolerances, with exact corrected
  alignments, chain selections and numbering.
- Added direct scalar-vs-batched alignment checks for all four structures in
  both modes, encoder parity around bucket and neighbor-count boundaries, and
  a tied-coordinate test excluding every padded neighbor.
- Repository-wide Black, isort, flake8 (including bugbear/comprehensions), and
  `git diff --check` passed. Tool versions match the pre-commit configuration.

All timing samples, rejected experiment patches and baseline arrays are kept
with this report so that each retention decision can be reviewed independently.
