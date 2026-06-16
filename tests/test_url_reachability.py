"""
Network-gated URL reachability tests.

Excluded from the default CI run. To run: pytest tests/test_url_reachability.py -m network
Or via a scheduled GitHub Actions workflow with '-m network'.
"""
import json
import pathlib

import pytest
import requests

ROOT = pathlib.Path(__file__).parent.parent
RESOURCES_PATH = ROOT / "setup" / "resources.json"

pytestmark = pytest.mark.network

ACCEPTABLE_CODES = {200, 206, 301, 302, 303, 307, 308}
TIMEOUT = 15


def _all_urls():
    data = json.loads(RESOURCES_PATH.read_text())
    pairs = []
    for name, entry in data.items():
        url = entry["url"]
        if isinstance(url, list):
            for u in url:
                pairs.append((name, u))
        else:
            pairs.append((name, url))
    return pairs


@pytest.mark.parametrize("name,url", _all_urls(), ids=[f"{n}" for n, _ in _all_urls()])
def test_url_is_reachable(name, url):
    try:
        resp = requests.head(url, timeout=TIMEOUT, allow_redirects=True)
        assert resp.status_code in ACCEPTABLE_CODES, (
            f"{name} ({url}): unexpected HTTP status {resp.status_code}"
        )
    except requests.RequestException as exc:
        pytest.fail(f"{name} ({url}): connection error — {exc}")
