import argparse
import sys
import disco
from disco.decomp import decomp_main
from disco.concat import concat_main


def build_parser():
    parser = argparse.ArgumentParser(
        prog="disco", description=f"DISCO {disco.__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    # Decomp parser
    decomp = subparsers.add_parser("decomp", help="decomposes gene trees")
    decomp.add_argument(
        "-i", "--input", type=str, help="Input tree list file", required=True
    )
    decomp.add_argument("-o", "--output", type=str, help="Output tree list file")
    decomp.add_argument(
        "-d",
        "--delimiter",
        type=str,
        help="Delimiter separating species name from rest of leaf label",
    )
    decomp.add_argument(
        "-n",
        "--nth-delimiter",
        type=int,  # Default is 1 -- set below
        help="Split on nth delimiter (only works with -d)",
    )
    decomp.add_argument(
        "-m", "--minimum", type=int, help="Minimum tree size outputted", default=4
    )
    decomp.add_argument(
        "-v", "--verbose", action="store_true", help="Enables verbose output"
    )
    decomp.add_argument(
        "--keep-labels",
        action="store_true",
        help="Keep original leaf labels instead of relabeling them with their species labels (only relevant with delimiter)",
    )
    decomp.add_argument(
        "--single_tree", action="store_true", help="Only output single large tree"
    )
    decomp.add_argument(
        "--no-decomp",
        action="store_true",
        help="Outputs rooted trees without decomposition",
    )
    decomp.add_argument(
        "--outgroups",
        action="store_true",
        help="Output outgroups to file (including ties)",
    )
    decomp.add_argument(
        "--remove_in_paralogs",
        action="store_true",
        help="Remove in-paralogs before rooting/scoring tree.",
    )

    # Concat parser
    concat = subparsers.add_parser(
        "concat", help="generate concatenation files from gene-family trees"
    )
    concat.add_argument(
        "-i", "--input", type=str, help="input tree list file", required=True
    )
    concat.add_argument(
        "-o", "--output-prefix", type=str, required=True, help="output tree list file"
    )
    concat.add_argument(
        "-a", "--alignment", required=True, type=str, help="alignment files list"
    )
    concat.add_argument(
        "-f",
        "--format",
        choices=["phylip", "fasta"],
        required=True,
        help="alignment file format",
    )
    concat.add_argument(
        "-d",
        "--delimiter",
        type=str,
        default="_",
        help="delimiter separating taxon label from the rest of the leaf label.",
    )
    concat.add_argument(
        "-m",
        "--filter",
        type=int,
        default=4,
        help="exclude decomposed trees with less then m taxa",
    )
    concat.add_argument(
        "-p", "--partition", action="store_true", help="generate partition file"
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    match args.command:
        case "concat":
            concat_main(args)
        case "decomp", _:
            if args.delimiter is None:
                if args.nth_delimiter is not None:
                    parser.error("Cannot set -n without a delimiter")
                if args.keep_labels:
                    parser.error("Cannot use --keep-labels without a delimiter")
            elif args.nth_delimiter is None:
                args.nth_delimiter = 1
            if args.single_tree and args.no_decomp:
                parser.error("Cannot combine --single_tree and --no-decomp")
            if not args.verbose and args.remove_in_paralogs:
                print(
                    "--remove_in_paralogs is meaningless without --verbose, as it does not change the optimal rooting. "
                    + "It may also slow the program."
                )
            decomp_main(args)
