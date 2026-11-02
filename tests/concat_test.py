import re

import treeswift as ts
from Bio import AlignIO

from disco.concat import concat_main
from disco.decomp import decomp_main, decompose, get_min_root, reroot_on_edge, tag
from tests.conftest import EXAMPLE_DIR, REPO_ROOT, run_baseline

ALIGNMENT_LIST = EXAMPLE_DIR / "seq_list.csv"
GTREES = EXAMPLE_DIR / "g_100.trees"

PARTITION_LINE_RE = re.compile(r"^.+, \S+=\d+-\d+$")


# --- Equivalence with v1.4.1 (last released version) ---


def test_concat_partitions_match_v1_4_1(baseline_dir, tmp_path):
    # Only the partition file is compared: v1.4.1 rerooted *at* the min-score
    # node rather than on its incident edge, which occasionally gave a
    # non-optimal tagging, so its alignments are not a valid reference (see
    # test_concat_uses_same_decomposition_as_decomp below).
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
    concat_main(
        str(GTREES),
        str(ALIGNMENT_LIST),
        str(new_prefix),
        "_",
        "phylip",
        filter=4,
        partition=True,
    )

    old_partitions = tmp_path / "old-partitions.txt"
    new_partitions = tmp_path / "new-partitions.txt"
    assert new_partitions.read_text() == old_partitions.read_text()


def test_concat_uses_same_decomposition_as_decomp(tmp_path):
    decomp_out = tmp_path / "decomp.trees"
    decomp_main(str(GTREES), str(decomp_out), delimiter="_", minimum=4)
    expected = [
        set(t.labels(True, False)) for t in ts.read_tree_newick(str(decomp_out))
    ]

    prefix = tmp_path / "out"
    concat_main(str(GTREES), str(ALIGNMENT_LIST), str(prefix), "_", "phylip", filter=4)
    aln = AlignIO.read(str(prefix) + "-aln.phy", "phylip")

    # Each decomposed tree contributes one block of columns as wide as its gene
    # alignment; the taxa with non-gap data in that block are the tree's leaves.
    with open(ALIGNMENT_LIST) as f:
        aln_files = [line.split(",")[0] for line in f if line.strip()]
    actual, col = [], 0
    for aln_file, n_trees in zip(aln_files, _trees_per_gene()):
        gene_len = AlignIO.read(str(REPO_ROOT / aln_file), "phylip").get_alignment_length()
        for _ in range(n_trees):
            block = aln[:, col : col + gene_len]
            actual.append({r.id for r in block if set(str(r.seq)) != {"-"}})
            col += gene_len
    assert col == aln.get_alignment_length()
    assert actual == expected


def _trees_per_gene():
    counts = []
    label_to_species = lambda x: x.split("_")[0]
    for tree in ts.read_tree_newick(str(GTREES)):
        reroot_on_edge(tree, get_min_root(tree, label_to_species)[0])
        tag(tree, label_to_species)
        counts.append(
            sum(1 for t in decompose(tree) if t.num_nodes(internal=False) >= 4)
        )
    return counts


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
