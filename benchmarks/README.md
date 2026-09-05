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
