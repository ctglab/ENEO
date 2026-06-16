"""
Validate setup/resources.json structure and completeness.
"""
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).parent.parent
RESOURCES_PATH = ROOT / "setup" / "resources.json"

VALID_FILETYPES = {"vcf", "fasta", "gtf", "table", "archive", "model", "sortmerna"}


def test_resources_json_parses():
    data = json.loads(RESOURCES_PATH.read_text())
    assert isinstance(data, dict) and len(data) > 0


def test_each_entry_has_required_fields(resources_json):
    for name, entry in resources_json.items():
        assert "filetype" in entry, f"{name}: missing 'filetype'"
        assert "url" in entry, f"{name}: missing 'url'"


def test_filetypes_are_known(resources_json):
    for name, entry in resources_json.items():
        ft = entry["filetype"].lower()
        assert ft in VALID_FILETYPES, f"{name}: unknown filetype '{ft}'"


def test_urls_are_nonempty_strings_or_lists(resources_json):
    for name, entry in resources_json.items():
        url = entry["url"]
        if isinstance(url, list):
            assert len(url) > 0, f"{name}: url list is empty"
            for u in url:
                assert isinstance(u, str) and u.strip(), (
                    f"{name}: url list contains empty or non-string entry"
                )
        else:
            assert isinstance(url, str) and url.strip(), (
                f"{name}: url is empty or not a string"
            )


def test_sortmerna_entries_have_keep_file(resources_json):
    for name, entry in resources_json.items():
        if entry["filetype"].lower() == "sortmerna":
            assert "keep_file" in entry, f"{name}: sortmerna entry missing 'keep_file'"
            assert entry["keep_file"].strip(), f"{name}: keep_file is empty"


def test_model_entries_have_list_of_urls(resources_json):
    for name, entry in resources_json.items():
        if entry["filetype"].lower() == "model":
            assert isinstance(entry["url"], list), (
                f"{name}: model filetype should have a list of urls"
            )
