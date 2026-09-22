import treeswift

from disco.decomp import decompose, decomp_main, get_min_root, reroot_on_edge, tag
from tests.conftest import EXAMPLE_DIR, run_baseline


# --- Equivalence with v1.4.1 (last released version) ---


def test_decomp_matches_v1_4_1_no_delimiter(baseline_dir, tmp_path):
    old_out = tmp_path / "old.trees"
    run_baseline(
        baseline_dir,
        "disco.py",
        "-i",
        str(EXAMPLE_DIR / "g_100-mult.trees"),
        "-o",
        str(old_out),
    )

    new_out = tmp_path / "new.trees"
    decomp_main(str(EXAMPLE_DIR / "g_100-mult.trees"), output_file_name=str(new_out))

    assert new_out.read_text() == old_out.read_text()


def test_decomp_matches_v1_4_1_with_delimiter(baseline_dir, tmp_path):
    old_out = tmp_path / "old.trees"
    run_baseline(
        baseline_dir,
        "disco.py",
        "-i",
        str(EXAMPLE_DIR / "g_100.trees"),
        "-o",
        str(old_out),
        "-d",
        "_",
    )

    new_out = tmp_path / "new.trees"
    decomp_main(
        str(EXAMPLE_DIR / "g_100.trees"), output_file_name=str(new_out), delimiter="_"
    )

    assert new_out.read_text() == old_out.read_text()


# --- Core algorithm unit tests ---


def test_tag_and_decompose_detects_duplication():
    gene_to_species = lambda x: x.split("_")[0]
    tree = treeswift.read_tree_newick("((A_1,A_2),B_1);")
    root, _, _ = get_min_root(tree, gene_to_species)
    reroot_on_edge(tree, root)
    tag(tree, gene_to_species)

    assert tree.n_dup == 1

    trees = decompose(tree)
    sizes = sorted(t.num_nodes(internal=False) for t in trees)
    assert sizes == [1, 2]


def test_decomp_main_minimum_filters_small_trees(tmp_path):
    input_file = tmp_path / "in.trees"
    input_file.write_text("((A_1,A_2),B_1);\n")

    out_min4 = tmp_path / "out_min4.trees"
    decomp_main(
        str(input_file), output_file_name=str(out_min4), delimiter="_", minimum=4
    )
    assert out_min4.read_text() == ""

    out_min1 = tmp_path / "out_min1.trees"
    decomp_main(
        str(input_file), output_file_name=str(out_min1), delimiter="_", minimum=1
    )
    assert len(out_min1.read_text().splitlines()) == 2

    out_min2 = tmp_path / "out_min2.trees"
    decomp_main(
        str(input_file), output_file_name=str(out_min2), delimiter="_", minimum=2
    )
    assert len(out_min2.read_text().splitlines()) == 1


def test_decomp_main_no_decomp_skips_decomposition(tmp_path):
    input_file = tmp_path / "in.trees"
    input_file.write_text("((A_1,A_2),B_1);\n")

    out = tmp_path / "out.trees"
    decomp_main(
        str(input_file), output_file_name=str(out), delimiter="_", no_decomp=True
    )

    lines = out.read_text().splitlines()
    assert len(lines) == 1
    tree = treeswift.read_tree_newick(lines[0])
    assert len(list(tree.labels(True, False))) == 3


def test_decomp_main_single_tree_outputs_only_largest(tmp_path):
    input_file = tmp_path / "in.trees"
    input_file.write_text("((A_1,A_2),B_1);\n")

    out = tmp_path / "out.trees"
    decomp_main(
        str(input_file),
        output_file_name=str(out),
        delimiter="_",
        minimum=1,
        single_tree=True,
    )

    lines = out.read_text().splitlines()
    assert len(lines) == 1
    tree = treeswift.read_tree_newick(lines[0])
    assert len(list(tree.labels(True, False))) == 2


def test_decomp_main_keep_labels_preserves_original_leaf_names(tmp_path):
    input_file = tmp_path / "in.trees"
    input_file.write_text("((A_1,A_2),B_1);\n")

    relabeled = tmp_path / "relabeled.trees"
    decomp_main(
        str(input_file), output_file_name=str(relabeled), delimiter="_", minimum=1
    )
    assert "A_1" not in relabeled.read_text()
    assert "A" in relabeled.read_text()

    kept = tmp_path / "kept.trees"
    decomp_main(
        str(input_file),
        output_file_name=str(kept),
        delimiter="_",
        minimum=1,
        keep_labels=True,
    )
    assert "A_1" in kept.read_text()


def test_decomp_main_outgroups_writes_outgroup_file(tmp_path):
    input_file = tmp_path / "in.trees"
    input_file.write_text("((A_1,A_2),B_1);\n")

    decomp_main(
        str(input_file),
        output_file_name=str(tmp_path / "out.trees"),
        delimiter="_",
        minimum=1,
        outgroups=True,
    )

    outgroup_file = tmp_path / "in_outgroups.txt"
    assert outgroup_file.exists()
    assert outgroup_file.read_text() == "Tree 1:\n{A}\n"
