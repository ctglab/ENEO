"""
Validate config/config_main.yaml structure and cross-references with resources.json.
"""
import pathlib

import pytest
import yaml

ROOT = pathlib.Path(__file__).parent.parent

REQUIRED_TOP_LEVEL_KEYS = {
    "OUTPUT_FOLDER", "TEMP_DIR", "datadirs", "params", "resources", "execution_mode"
}
VALID_EXECUTION_MODES = {"full", "CI"}
REQUIRED_PARAM_SECTIONS = {
    "BQSR", "deepvariant", "fastp", "gatk", "MarkDuplicates", "pMHC",
    "STAR", "SplitNCigarReads", "salmon", "samtools", "strelka2", "t1k", "vcfanno", "vep",
}


def test_config_parses():
    data = yaml.safe_load((ROOT / "config" / "config_main.yaml").read_text())
    assert isinstance(data, dict) and len(data) > 0


def test_required_top_level_keys(config):
    missing = REQUIRED_TOP_LEVEL_KEYS - set(config.keys())
    assert not missing, f"Missing top-level keys: {missing}"


def test_execution_mode_is_valid(config):
    mode = config.get("execution_mode")
    assert mode in VALID_EXECUTION_MODES, (
        f"execution_mode '{mode}' is not one of {VALID_EXECUTION_MODES}"
    )


def test_params_sections_present(config):
    params = config.get("params", {})
    missing = REQUIRED_PARAM_SECTIONS - set(params.keys())
    assert not missing, f"Missing params sections: {missing}"


def test_all_downloadable_resources_in_config(config, resources_json):
    config_resources = set(config.get("resources", {}).keys())
    json_resources = set(resources_json.keys())
    missing = json_resources - config_resources
    assert not missing, (
        f"Resources defined in resources.json but absent from config.resources: {missing}"
    )


def test_datadirs_has_logs_section(config):
    assert "logs" in config.get("datadirs", {}), "datadirs is missing 'logs' section"


def test_in_repo_resources_exist(config):
    """Resources whose paths start with 'workflow/' must exist on disk."""
    resources = config.get("resources", {})
    missing = []
    for name, path in resources.items():
        if isinstance(path, str) and path.startswith("workflow/"):
            full_path = ROOT / path
            if not full_path.exists():
                missing.append(f"{name}: {path}")
    assert not missing, "In-repo resources missing from disk:\n" + "\n".join(missing)
