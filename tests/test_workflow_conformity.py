"""
Structural conformity tests for the Snakemake workflow.

Parses .smk files as plain text — no Snakemake installation required.
Checks that every conda env, script, and include path referenced in rules
actually exists on disk, and that no rule names are duplicated.
"""
import pathlib
import re
from collections import Counter

import pytest

ROOT = pathlib.Path(__file__).parent.parent
WORKFLOW_DIR = ROOT / "workflow"
RULES_DIR = WORKFLOW_DIR / "rules"
ENVS_DIR = WORKFLOW_DIR / "envs"
SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
SNAKEFILE = WORKFLOW_DIR / "Snakefile"

# Matches both inline and indented quoted string after a directive keyword:
#   conda: "../envs/star.yml"
#   conda:
#       "../envs/star.yml"
_DIRECTIVE_RE = re.compile(r'{keyword}:\s*\n?\s*"([^"]+)"')


def _extract(text, keyword):
    pattern = re.compile(rf'{keyword}:\s*\n?\s*"([^"]+)"')
    return pattern.findall(text)


def test_snakefile_exists():
    assert SNAKEFILE.exists(), f"Snakefile not found: {SNAKEFILE}"


def test_rules_dir_exists():
    assert RULES_DIR.is_dir(), f"Rules directory not found: {RULES_DIR}"


def test_all_includes_exist():
    text = SNAKEFILE.read_text()
    paths = _extract(text, "include")
    missing = [p for p in paths if not (WORKFLOW_DIR / p).exists()]
    assert not missing, f"Included files not found: {missing}"


def test_all_rule_files_are_included():
    text = SNAKEFILE.read_text()
    included_names = {p.split("/")[-1] for p in _extract(text, "include")}
    rule_file_names = {f.name for f in RULES_DIR.glob("*.smk")}
    missing = rule_file_names - included_names
    assert not missing, f"Rule files not included in Snakefile: {missing}"


def test_all_conda_envs_exist(rule_files):
    missing = []
    for rule_file in rule_files:
        text = rule_file.read_text()
        for env_path in _extract(text, "conda"):
            # Paths in rules are relative to workflow/rules/
            full = (RULES_DIR / env_path).resolve()
            if not full.exists():
                missing.append(f"{rule_file.name}: {env_path}")
    assert not missing, "Missing conda env files:\n" + "\n".join(missing)


def test_all_scripts_exist(rule_files):
    missing = []
    for rule_file in rule_files:
        text = rule_file.read_text()
        for script_path in _extract(text, "script"):
            full = (RULES_DIR / script_path).resolve()
            if not full.exists():
                missing.append(f"{rule_file.name}: {script_path}")
    assert not missing, "Missing script files:\n" + "\n".join(missing)


def test_no_duplicate_rule_names(rule_files):
    all_names = []
    for rule_file in rule_files:
        names = re.findall(r'^rule\s+(\w+)\s*:', rule_file.read_text(), re.MULTILINE)
        all_names.extend(names)
    # Also check the top-level Snakefile
    top_names = re.findall(r'^rule\s+(\w+)\s*:', SNAKEFILE.read_text(), re.MULTILINE)
    all_names.extend(top_names)

    counts = Counter(all_names)
    duplicates = {name: n for name, n in counts.items() if n > 1}
    assert not duplicates, f"Duplicate rule names: {duplicates}"


def test_env_files_are_valid_yaml(rule_files):
    """Every conda env .yml referenced in rules must parse as valid YAML."""
    import yaml

    broken = []
    seen = set()
    for rule_file in rule_files:
        for env_path in _extract(rule_file.read_text(), "conda"):
            full = (RULES_DIR / env_path).resolve()
            if full in seen or not full.exists():
                continue
            seen.add(full)
            try:
                yaml.safe_load(full.read_text())
            except yaml.YAMLError as exc:
                broken.append(f"{full.name}: {exc}")
    assert not broken, "Invalid YAML in conda env files:\n" + "\n".join(broken)
