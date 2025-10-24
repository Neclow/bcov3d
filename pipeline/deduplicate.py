"""Remove duplicate sequences from a FASTA file."""

from argparse import ArgumentParser

from src.fasta import deduplicate


def parse_args():
    parser = ArgumentParser(description="Remove duplicate sequences from a FASTA file.")
    parser.add_argument("infile", help="Input FASTA file")
    parser.add_argument("outfile", help="Output FASTA file")
    parser.add_argument(
        "--verbose", action="store_true", help="If True, print additional information."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    deduplicate(args.infile, args.outfile, verbose=args.verbose)
