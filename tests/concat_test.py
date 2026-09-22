import re

import treeswift as ts
from Bio import AlignIO

from disco.concat import concat_main
from tests.conftest import EXAMPLE_DIR, run_baseline

ALIGNMENT_LIST = EXAMPLE_DIR / "seq_list.csv"
GTREES = EXAMPLE_DIR / "g_100.trees"

PARTITION_LINE_RE = re.compile(r"^.+, \S+=\d+-\d+$")


# --- Equivalence with v1.4.1 (last released version) ---


def test_concat_matches_v1_4_1(baseline_dir, tmp_path):
    old_prefix = tmp_path / "old"
    run_baseline(
        baseline_dir,
        "ca_disco.py",
        "-i",
        str(GTREES),
        "-o",
        str(old_prefix),
        "-a",
        str(ALIGNMENT_LIST),
        "-f",
        "phylip",
        "-p",
    )

    new_prefix = tmp_path / "new"
    # explicit filter=4 to match the old script's -m default (concat_main's own
    # default is 2 -- only the CLI, which always passes -m explicitly, avoids
    # this discrepancy; see disco/cli.py's -m/--filter default of 4).
    concat_main(
        str(GTREES),
        str(ALIGNMENT_LIST),
        str(new_prefix),
        "_",
        "phylip",
        filter=4,
        partition=True,
    )

    old_aln = tmp_path / "old-aln.phy"
    new_aln = tmp_path / "new-aln.phy"
    assert new_aln.read_text() == old_aln.read_text()

    old_partitions = tmp_path / "old-partitions.txt"
    new_partitions = tmp_path / "new-partitions.txt"
    assert new_partitions.read_text() == old_partitions.read_text()


# --- Unit tests ---


def test_concat_main_alignment_has_one_row_per_taxon(tmp_path):
    label_to_species = lambda x: x.split("_")[0]
    trees = ts.read_tree_newick(str(GTREES))
    expected_taxa = set(
        label_to_species(label)
        for tree in trees
        for label in tree.labels(True, False)
    )

    prefix = tmp_path / "out"
    concat_main(str(GTREES), str(ALIGNMENT_LIST), str(prefix), "_", "phylip", filter=4)

    aln = AlignIO.read(str(prefix) + "-aln.phy", "phylip")
    assert len(aln) == len(expected_taxa)


def test_concat_main_partition_file_format(tmp_path):
    prefix = tmp_path / "out"
    concat_main(
        str(GTREES),
        str(ALIGNMENT_LIST),
        str(prefix),
        "_",
        "phylip",
        filter=4,
        partition=True,
    )

    with open(str(prefix) + "-partitions.txt") as f:
        partition_lines = f.read().splitlines()

    assert partition_lines
    for line in partition_lines:
        assert PARTITION_LINE_RE.match(line), line


def test_concat_main_filter_excludes_small_decomposed_trees(tmp_path):
    strict_prefix = tmp_path / "strict"
    concat_main(
        str(GTREES), str(ALIGNMENT_LIST), str(strict_prefix), "_", "phylip", filter=4
    )
    strict_aln = AlignIO.read(str(strict_prefix) + "-aln.phy", "phylip")

    lenient_prefix = tmp_path / "lenient"
    concat_main(
        str(GTREES), str(ALIGNMENT_LIST), str(lenient_prefix), "_", "phylip", filter=2
    )
    lenient_aln = AlignIO.read(str(lenient_prefix) + "-aln.phy", "phylip")

    # a lower minimum-taxa filter admits more (or equally many) decomposed
    # trees per gene, so the concatenated alignment can only grow.
    assert lenient_aln.get_alignment_length() > strict_aln.get_alignment_length()
