"""Sample a subset of sequences from two FASTA files and write them to separate output files."""

import os

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import pandas as pd

from src.fasta import pgap, filter_taxa, filter_sites, unalign

MAX_CHARS = {
    "dna": 4,
    "rna": 4,
    "aa": 20,
    "3di": 20,
}


def parse_args():
    parser = ArgumentParser(
        description=(
            "Sample a subset of sequences from two FASTA files "
            "and write them to separate output files."
        ),
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "n", type=int, help="Number of sequences to sample from each input file"
    )
    parser.add_argument("input1", help="First input FASTA file")
    parser.add_argument("input2", help="Second input FASTA file")
    parser.add_argument("metadata", help="Metadata file (.csv)")
    # parser.add_argument("output1", help="Output FASTA file for sequences from input1")
    # parser.add_argument("output2", help="Output FASTA file for sequences from input2")
    parser.add_argument(
        "--by",
        default="Virus,MajorVariant",
        help="Comma-separated column name(s) to group by in the metadata file",
    )
    parser.add_argument(
        "--col",
        default="PDB",
        help="Column name for sequence IDs in the metadata file",
    )
    parser.add_argument(
        "--filtersites",
        type=int,
        help="Filter out highly variable sites in the sequences before sampling",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (used with --random)",
    )
    parser.add_argument(
        "--dtype",
        type=str,
        default="aa",
        choices=("dna", "rna", "aa", "3di"),
        help="Type of sequences in the input files",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    start_df = pd.read_csv(args.metadata).set_index(args.col)

    # Calculate the proportion of gaps in each sequence
    start_df["pgap1"] = pgap(args.input1)
    start_df["pgap2"] = pgap(args.input2)
    start_df["mpgap"] = start_df[["pgap1", "pgap2"]].mean(axis=1)
    start_df["has_receptor_or_antibody"] = (
        start_df["Receptor/Other"].notna() | start_df["Antibody"].notna()
    )

    start_df.dropna(subset=["mpgap"], inplace=True)

    # Sort by absence of receptor + structure resolution
    start_df.sort_values(
        by=["has_receptor_or_antibody", "Resolution"],
        ascending=True,
        inplace=True,
    )
    # elif args.random:
    #     # Shuffle the DataFrame
    #     start_df.sample(frac=1, random_state=args.seed, inplace=True)

    subset_df = start_df.groupby(args.by.split(","), sort=False, dropna=False).head(
        args.n
    )

    # Save the subset metadata to a CSV file
    subset_df.to_csv(Path(args.metadata).parent / f"subset{args.n}.csv")

    # Filter sequences from both input files
    ids_to_sample = subset_df.index.tolist()
    dir1 = os.path.dirname(args.input1)
    dir2 = os.path.dirname(args.input2)
    output1 = f"{dir1}/subset{args.n}.fa"
    output2 = f"{dir2}/subset{args.n}.fa"
    filter_taxa(args.input1, output1, ids_to_sample)
    filter_taxa(args.input2, output2, ids_to_sample)
    if args.filtersites is not None:
        assert 1 < args.filtersites <= MAX_CHARS[args.dtype]
        filter_sites(
            output1,
            output1,
            colfile=f"{dir1}/subset{args.n}.colnumbering",
            max_unique=args.filtersites,
        )
        filter_sites(
            output2,
            output2,
            colfile=f"{dir2}/subset{args.n}.colnumbering",
            max_unique=args.filtersites,
        )

    print(f"(subset) Remaining sequences: {len(ids_to_sample)}")

    stem1, ext1 = os.path.splitext(output1)
    unalign(output1, f"{stem1}_unaligned{ext1}", f"{stem1}_gaps.json")
    stem2, ext2 = os.path.splitext(output2)
    unalign(output2, f"{stem2}_unaligned{ext2}", f"{stem2}_gaps.json")
