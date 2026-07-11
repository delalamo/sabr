import ast
from pathlib import Path

SOURCE = Path(__file__).parents[1] / "src" / "sabr"


def _walk_with_scope(node, inside_function=False):
    for child in ast.iter_child_nodes(node):
        child_inside = inside_function or isinstance(
            child, (ast.FunctionDef, ast.AsyncFunctionDef)
        )
        yield child, inside_function
        yield from _walk_with_scope(child, child_inside)


def test_no_lazy_imports_or_variable_annotations():
    for path in SOURCE.glob("*.py"):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node, inside_function in _walk_with_scope(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                assert (
                    not inside_function
                ), f"Lazy import in {path}: {node.lineno}"
            if isinstance(node, ast.AnnAssign):
                assert (
                    inside_function
                ), f"Variable/class annotation in {path}: {node.lineno}"


def test_removed_legacy_surface_is_absent():
    source = "\n".join(
        path.read_text()
        for path in SOURCE.glob("*.py")
        if path.name != "constants.py"
    )
    for removed in (
        "correct_fr1",
        "detect_chain_type",
        "reference_chain_type",
        "disable_custom_gap_penalties",
        "disable_deterministic_renumbering",
    ):
        assert removed not in source
