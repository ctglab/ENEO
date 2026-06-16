"""
Static analysis of setup/download_res.py: syntax, function presence, and CLI contract.

Uses ast.parse rather than importing the module to avoid heavy dependencies
(pandas, rich) in the CI environment.
"""
import ast
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).parent.parent
SCRIPT_PATH = ROOT / "setup" / "download_res.py"

EXPECTED_FUNCTIONS = {
    "parse_arguments",
    "download_resource",
    "decompress_file",
    "convert_notations",
    "generate_allele_frequency",
    "create_sequence_dictionary",
    "download_deepvariant_model_files",
    "download_sortmerna_db",
    "convert_REDI",
    "main",
}

EXPECTED_CLI_ARGS = {"--outfolder", "--config", "--resources", "--dry-run", "--workers"}


@pytest.fixture(scope="module")
def script_src():
    return SCRIPT_PATH.read_text()


@pytest.fixture(scope="module")
def script_ast(script_src):
    return ast.parse(script_src)


def test_script_exists():
    assert SCRIPT_PATH.exists(), f"Script not found: {SCRIPT_PATH}"


def test_script_has_valid_python_syntax(script_src):
    try:
        ast.parse(script_src)
    except SyntaxError as exc:
        pytest.fail(f"Syntax error in download_res.py: {exc}")


def test_expected_functions_are_defined(script_ast):
    defined = {
        node.name
        for node in ast.walk(script_ast)
        if isinstance(node, ast.FunctionDef)
    }
    missing = EXPECTED_FUNCTIONS - defined
    assert not missing, f"Functions missing from download_res.py: {missing}"


def test_cli_args_are_declared(script_src):
    for arg in EXPECTED_CLI_ARGS:
        assert arg in script_src, (
            f"CLI argument '{arg}' not found in parse_arguments"
        )


def test_dry_run_uses_store_true(script_src):
    assert "store_true" in script_src, (
        "--dry-run should use action='store_true'"
    )


def test_outfolder_is_required(script_src):
    assert "required=True" in script_src, (
        "--outfolder should be a required argument"
    )


def test_filetype_handlers_cover_resources_json(script_src):
    """Every filetype in resources.json must have a corresponding handler branch."""
    resources = json.loads((ROOT / "setup" / "resources.json").read_text())
    filetypes = {entry["filetype"].lower() for entry in resources.values()}
    for ft in filetypes:
        assert ft in script_src, (
            f"Filetype '{ft}' from resources.json has no handler in download_res.py"
        )


def test_main_guarded_by_name_check(script_src):
    assert '__name__' in script_src and '__main__' in script_src, (
        "Script should use 'if __name__ == \"__main__\"' guard"
    )
