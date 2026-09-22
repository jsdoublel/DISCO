import subprocess
import sys

from tests.conftest import EXAMPLE_DIR

GTREES = EXAMPLE_DIR / "g_100.trees"


def run_cli(*args):
    return subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from disco.cli import main; sys.exit(main())",
            *args,
        ],
        capture_output=True,
        text=True,
    )


def test_nth_delimiter_without_delimiter_errors():
    result = run_cli("decomp", "-i", str(GTREES), "-n", "2")
    assert result.returncode != 0
    assert "Cannot set -n without a delimiter" in result.stderr


def test_keep_labels_without_delimiter_errors():
    result = run_cli("decomp", "-i", str(GTREES), "--keep-labels")
    assert result.returncode != 0
    assert "Cannot use --keep-labels without a delimiter" in result.stderr


def test_single_tree_and_no_decomp_together_errors():
    result = run_cli(
        "decomp", "-i", str(GTREES), "--single_tree", "--no-decomp"
    )
    assert result.returncode != 0
    assert "Cannot combine --single_tree and --no-decomp" in result.stderr
