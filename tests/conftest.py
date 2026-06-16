import json
import pathlib

import pytest
import yaml

ROOT = pathlib.Path(__file__).parent.parent


@pytest.fixture(scope="session")
def root():
    return ROOT


@pytest.fixture(scope="session")
def resources_json():
    return json.loads((ROOT / "setup" / "resources.json").read_text())


@pytest.fixture(scope="session")
def config():
    return yaml.safe_load((ROOT / "config" / "config_main.yaml").read_text())


@pytest.fixture(scope="session")
def rule_files():
    return list((ROOT / "workflow" / "rules").glob("*.smk"))
