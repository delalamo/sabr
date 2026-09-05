"""Alternate old/new implementations to reduce machine-load timing bias.

Run after each isolated stage: python benchmarks/paired.py --before REV
--output PATH [--case heavy] [--mode sabr]. Uses the same fixtures as inference.
Only model.py and alignment.py are loaded from REV; structure handling and
numbering are unchanged in these experiments. No persistent cache is used.
"""

import argparse
import importlib.util
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

from inference import CASES, ROOT


def worker(args):
    import jax
    import numpy as np
    from Bio.PDB import PDBParser

    from sabr import alignment, model
    from sabr.numbering import number_alignment
    from sabr.structure import apply_numbering, extract_chain

    def load_module(name, directory):
        path = Path(directory) / f"{name}.py"
        path.write_bytes(
            subprocess.check_output(
                ["git", "show", f"{args.before}:src/sabr/{name}.py"], cwd=ROOT
            )
        )
        spec = importlib.util.spec_from_file_location(f"before_{name}", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module

    with tempfile.TemporaryDirectory(prefix="sabr-paired-") as directory:
        old_model = load_module("model", directory)
        old_alignment = load_module("alignment", directory)
    filename, chain, candidates = CASES[args.case]
    structure = PDBParser(QUIET=True).get_structure(
        args.case, ROOT / "tests" / "data" / filename
    )
    warnings.filterwarnings("ignore", category=UserWarning)
    data = extract_chain(structure, chain, None)
    implementations = {
        "before": (
            getattr(old_model, "_encode_device", old_model.encode),
            old_alignment.align,
        ),
        "after": (
            getattr(model, "_encode_device", model.encode),
            alignment.align,
        ),
    }

    if args.length_sweep:
        samples = {"before": [], "after": []}
        for index, length in enumerate(
            range(len(data.coords) - 9, len(data.coords) + 1)
        ):
            outputs = {}
            order = (
                ("before", "after") if index % 2 == 0 else ("after", "before")
            )
            for label in order:
                start = time.perf_counter()
                outputs[label] = implementations[label][0](
                    data.coords[:length], args.mode
                )
                jax.block_until_ready(outputs[label])
                samples[label].append(
                    {"length": length, "seconds": time.perf_counter() - start}
                )
            np.testing.assert_allclose(
                outputs["after"], outputs["before"], rtol=1e-5, atol=2e-6
            )
        result = {
            "before_revision": args.before,
            "case": args.case,
            "mode": args.mode,
            "samples": samples,
            "total": {
                label: sum(sample["seconds"] for sample in values)
                for label, values in samples.items()
            },
            "compiled_encoder_variants": {
                "before": old_model._APPLY_ENCODER._cache_size(),
                "after": model._APPLY_ENCODER._cache_size(),
            },
            "parity": "passed",
        }
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / f"length_sweep_{args.case}_{args.mode}.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(result, flush=True)
        return

    def run(label):
        encode, align = implementations[label]
        start = time.perf_counter()
        start_cpu = time.process_time()
        embedding = encode(data.coords, args.mode)
        jax.block_until_ready(embedding)
        encoded = time.perf_counter()
        encoded_cpu = time.process_time()
        matrix, selected, score = align(
            embedding, data.gap_indices, candidates, 0.0, mode=args.mode
        )
        jax.block_until_ready(matrix)
        aligned = time.perf_counter()
        aligned_cpu = time.process_time()
        numbers = number_alignment(matrix, data.sequence, "imgt", selected)
        result = apply_numbering(structure, chain, data, numbers)
        del result
        end = time.perf_counter()
        end_cpu = time.process_time()
        return (
            {
                "encode": encoded - start,
                "align": aligned - encoded,
                "pipeline": end - start,
                "cpu_encode": encoded_cpu - start_cpu,
                "cpu_align": aligned_cpu - encoded_cpu,
                "cpu_pipeline": end_cpu - start_cpu,
            },
            embedding,
            matrix,
            selected,
            score,
            numbers,
        )

    before = run("before")
    after = run("after")
    np.testing.assert_allclose(after[1], before[1], rtol=1e-5, atol=2e-6)
    np.testing.assert_array_equal(after[2], before[2])
    assert after[3] == before[3]
    np.testing.assert_allclose(after[4], before[4], rtol=1e-5, atol=1e-6)
    assert after[5] == before[5]
    samples = {"before": [], "after": []}
    # Warm both paths before the measured, alternating pairs.
    run("after")
    run("before")
    for index in range(args.repeat):
        order = ("before", "after") if index % 2 == 0 else ("after", "before")
        for label in order:
            samples[label].append(run(label)[0])
    medians = {
        label: {
            stage: statistics.median(sample[stage] for sample in timings)
            for stage in timings[0]
        }
        for label, timings in samples.items()
    }
    result = {
        "case": args.case,
        "mode": args.mode,
        "before_revision": args.before,
        "residues": len(data.coords),
        "medians": medians,
        "samples": samples,
        "speedup": {
            stage: medians["before"][stage] / medians["after"][stage]
            for stage in medians["before"]
        },
        "parity": "passed",
    }
    args.output.mkdir(parents=True, exist_ok=True)
    key = f"paired_{args.case}_{args.mode}"
    (args.output / f"{key}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(f"{key}: {result['speedup']}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", choices=CASES)
    parser.add_argument("--mode", choices=("sabr", "softalign"), default="sabr")
    parser.add_argument("--repeat", type=int, default=11)
    parser.add_argument("--length-sweep", action="store_true")
    args = parser.parse_args()
    if args.length_sweep and not args.case:
        parser.error("--length-sweep requires --case")
    os.environ["JAX_PLATFORMS"] = "cpu"
    os.environ["JAX_ENABLE_COMPILATION_CACHE"] = "false"
    if args.case:
        worker(args)
    else:
        for case in CASES:
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--before",
                    args.before,
                    "--output",
                    str(args.output),
                    "--case",
                    case,
                    "--mode",
                    args.mode,
                    "--repeat",
                    str(args.repeat),
                ],
                check=True,
                cwd=ROOT,
            )


if __name__ == "__main__":
    main()
