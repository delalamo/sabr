"""Run isolated, synchronized inference comparisons on bundled structures.

Example: python benchmarks/inference.py --output benchmarks/results/baseline
Later: add --reference benchmarks/results/baseline to check output parity.
Each case/mode runs in a fresh process. Run stages sequentially, never in
parallel. Persistent JAX caching is disabled unless --cache-dir is supplied.
"""

import argparse
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "heavy": ("test_heavy_chain.pdb", "F", "auto"),
    "kappa": ("12e8_L.pdb", "L", "auto"),
    "lambda": ("1bjm_A.pdb", "A", "auto"),
    "long_loop": ("8sve_L.pdb", "M", "K"),
}


def worker(args):
    start = time.perf_counter()
    import jax
    import numpy as np
    from Bio.PDB import PDBParser

    from sabr import renumber_structure
    from sabr.alignment import align
    from sabr.model import encode
    from sabr.numbering import number_alignment
    from sabr.structure import apply_numbering, extract_chain

    import_seconds = time.perf_counter() - start
    filename, chain, candidate = CASES[args.case]
    structure = PDBParser(QUIET=True).get_structure(
        args.case, ROOT / "tests" / "data" / filename
    )
    # This fixture intentionally exercises the documented structural-gap opt-in.
    warnings.filterwarnings("ignore", category=UserWarning)

    def timed(function, *values, **kwargs):
        start = time.perf_counter()
        result = function(*values, **kwargs)
        jax.block_until_ready(result)
        return result, time.perf_counter() - start

    def pipeline():
        durations = {}
        data, durations["extract"] = timed(extract_chain, structure, chain, None)
        embedding, durations["encode"] = timed(encode, data.coords, args.mode)
        alignment_result, durations["align"] = timed(
            align, embedding, data.gap_indices, candidate, 0.0, mode=args.mode
        )
        matrix, selected, score = alignment_result
        numbers, durations["number"] = timed(
            number_alignment, matrix, data.sequence, "imgt", selected
        )
        _, durations["copy"] = timed(
            apply_numbering, structure, chain, data, numbers
        )
        durations["total"] = sum(durations.values())
        return durations, embedding, matrix, selected, score, numbers, data

    first, embedding, matrix, selected, score, numbers, data = pipeline()
    samples = [pipeline()[0] for _ in range(args.repeat)]
    api_samples = []
    for _ in range(args.repeat):
        _, duration = timed(
            renumber_structure,
            structure,
            chain,
            chain_type=candidate,
            mode=args.mode,
            dangerously_allow_structural_gaps=True,
        )
        api_samples.append(duration)

    key = f"{args.case}_{args.mode}"
    args.output.mkdir(parents=True, exist_ok=True)
    if args.reference:
        with np.load(args.reference / f"{key}.npz") as baseline:
            np.testing.assert_allclose(
                embedding, baseline["embedding"], rtol=1e-5, atol=2e-6
            )
            np.testing.assert_array_equal(matrix, baseline["matrix"])
            np.testing.assert_allclose(
                score, baseline["score"], rtol=1e-5, atol=1e-6
            )
        prior = json.loads((args.reference / f"{key}.json").read_text())
        assert selected == prior["selected"]
        assert json.loads(json.dumps(numbers)) == prior["numbering"]
    np.savez_compressed(
        args.output / f"{key}.npz",
        embedding=np.asarray(embedding),
        matrix=matrix,
        score=score,
    )
    result = {
        "case": args.case,
        "mode": args.mode,
        "residues": len(data.coords),
        "candidate": candidate,
        "selected": selected,
        "numbering": numbers,
        "python": platform.python_version(),
        "jax": jax.__version__,
        "devices": [str(device) for device in jax.devices()],
        "platform": platform.platform(),
        "import_seconds": import_seconds,
        "first": first,
        "warm_median": {
            name: statistics.median(sample[name] for sample in samples)
            for name in first
        },
        "warm_samples": samples,
        "api_median": statistics.median(api_samples),
        "api_samples": api_samples,
        "parity": "passed" if args.reference else "baseline",
        "cache_dir": str(args.cache_dir) if args.cache_dir else None,
    }
    (args.output / f"{key}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"{key}: N={len(data.coords)} first={first['total']:.3f}s "
        f"encode={result['warm_median']['encode']:.4f}s "
        f"align={result['warm_median']['align']:.4f}s "
        f"API={result['api_median']:.4f}s parity={result['parity']}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--case", choices=CASES)
    parser.add_argument("--mode", choices=("sabr", "softalign"))
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    os.environ["JAX_PLATFORMS"] = "cpu"
    os.environ["JAX_ENABLE_COMPILATION_CACHE"] = (
        "true" if args.cache_dir else "false"
    )
    if args.cache_dir:
        os.environ["JAX_COMPILATION_CACHE_DIR"] = str(args.cache_dir.resolve())
        os.environ["JAX_PERSISTENT_CACHE_MIN_COMPILE_TIME_SECS"] = "0"
    if args.case and args.mode:
        worker(args)
        return
    for case in ([args.case] if args.case else CASES):
        for mode in ([args.mode] if args.mode else ("sabr", "softalign")):
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--output", str(args.output),
                "--repeat", str(args.repeat),
                "--case", case,
                "--mode", mode,
            ]
            if args.reference:
                command.extend(("--reference", str(args.reference)))
            if args.cache_dir:
                command.extend(("--cache-dir", str(args.cache_dir)))
            subprocess.run(command, check=True, cwd=ROOT)


if __name__ == "__main__":
    main()
