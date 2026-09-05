"""Run four focused mutations in subprocess memory, never editing sources."""

import os
import subprocess
import sys
from pathlib import Path

MUTATIONS = {
    "bypass renumbering": (
        "import sabr.cli as m; "
        "m.renumber_structure = lambda structure, *a, **k: structure",
        [
            "tests/test_integration.py::test_real_cli_round_trip",
            "tests/test_integration.py::test_real_softalign_cli_round_trip",
        ],
    ),
    "reverse altloc tie": (
        "import inspect, sabr.structure as m; "
        "source = inspect.getsource(m); old = '(-item[0], item[1])'; "
        "assert source.count(old) == 1; "
        "replacement = '(-item[0], tuple(-ord(c) for c in item[1]))'; "
        "exec(compile(source.replace(old, replacement), m.__file__, 'exec'), "
        "m.__dict__)",
        [
            "tests/test_structure_hardening.py::"
            "test_altloc_selection_prefers_blank_then_occupancy_then_name"
        ],
    ),
    "remove domain offsets": (
        "from sabr import constants; constants.DOMAIN_NUMBERING_STRIDE = 0",
        ["tests/test_multidomain.py", "-k", "synthetic and HHK and 1"],
    ),
    "suppress CDR2 correction": (
        "import sabr.corrections as m; original = m.correct_cdr_loop; "
        "m.correct_cdr_loop = lambda aln, name, *a, **k: "
        "aln if name == 'CDR2' else original(aln, name, *a, **k)",
        [
            "tests/test_corrections.py::"
            "test_gap_in_cdr1_does_not_suppress_cdr2_repair"
        ],
    ),
}


def main():
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["JAX_PLATFORMS"] = "cpu"
    failures = []
    for name, (mutation, tests) in MUTATIONS.items():
        code = (
            mutation
            + "; import pytest, sys; sys.exit(pytest.main(sys.argv[1:]))"
        )
        result = subprocess.run(
            [sys.executable, "-c", code, "-q", "--tb=short", *tests],
            cwd=root,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode != 1:
            failures.append(name)
            print(result.stdout)
            print(
                f"FAIL: {name}: expected test failure, "
                f"got exit {result.returncode}"
            )
        else:
            print(f"Detected: {name}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
