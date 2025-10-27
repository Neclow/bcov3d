"""Intersect two FASTA files to find common sequences."""

from argparse import ArgumentParser

from Bio import SeqIO


def parse_args():
    parser = ArgumentParser(
        description="Find common sequences between two FASTA files."
    )
    parser.add_argument("input1", help="First input FASTA file")
    parser.add_argument("input2", help="Second input FASTA file")
    parser.add_argument("output1", help="Output FASTA file for sequences from input1")
    parser.add_argument("output2", help="Output FASTA file for sequences from input2")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    keys = intersect(
        input1=args.input1,
        input2=args.input2,
        output1=args.output1,
        output2=args.output2,
    )
