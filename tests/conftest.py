import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = REPO_ROOT / "example"
BASELINE_TAG = "v1.4.1"


@pytest.fixture(scope="session")
def baseline_dir(tmp_path_factory):
    """Materializes the pre-packaging release (v1.4.1) scripts into an isolated dir."""
    dest = tmp_path_factory.mktemp("baseline")
    for name in ("disco.py", "ca_disco.py"):
        result = subprocess.run(
            ["git", "show", f"{BASELINE_TAG}:{name}"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            pytest.skip(f"could not extract {BASELINE_TAG}:{name} from git history")
        (dest / name).write_text(result.stdout)
    return dest


def run_baseline(baseline_dir, script, *args):
    # Run with cwd=REPO_ROOT (so relative paths in fixture files, e.g. the
    # alignment list CSV, resolve the same way they do for the new code) but
    # invoke the script by its absolute path: Python still puts the script's
    # own directory first on sys.path, so ca_disco.py's "from disco import *"
    # resolves to the sibling disco.py script rather than the installed
    # disco package.
    subprocess.run(
        [sys.executable, str(baseline_dir / script), *args],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
