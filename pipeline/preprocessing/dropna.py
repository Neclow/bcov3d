"""Drop sequences with a high proportion of gaps from a FASTA file."""

from argparse import ArgumentParser

from src.fasta import dropna


def parse_args():
    parser = ArgumentParser(
        description="Filter sequences in a FASTA file based on the proportion of gaps."
    )
    parser.add_argument("infile", type=str, help="Input FASTA file")
    parser.add_argument(
        "outfile", type=str, help="Output FASTA file with filtered sequences"
    )
    parser.add_argument(
        "--thresh",
        type=float,
        default=0.9,
        help="Threshold for the proportion of gaps (0.0 to 1.0)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    dropna(args.infile, args.outfile, args.thresh)
