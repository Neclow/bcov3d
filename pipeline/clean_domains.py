import json

from argparse import ArgumentParser
from functools import partial

import pandas as pd


def get_aligned_coords(x, gaps):
    gap_coords = gaps[x.id]
    aligned_start = x.start + sum(1 for g in gap_coords if g < x.start)
    aligned_end = x.end + sum(1 for g in gap_coords if g < x.end)
    return [aligned_start, aligned_end]


COLUMNS = [
    "id",
    "md5",
    "length",
    "signature_library",
    "signature_accession",
    "signature_name",
    "start",
    "end",
    "evalue",
    "t",
    "date",
    "entry_accession",
    "entry_description",
    "goXRefs",
    "pathwayXRefs",
]


def parse_args():
    parser = ArgumentParser(description="Clean domain data")
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input file containing domain data.",
    )
    parser.add_argument(
        "gaps",
        type=str,
        help="Path to the file containing gaps data.",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to save the cleaned domain data.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Load the domain data
    doms = (
        pd.read_csv(args.input, sep="\t", header=None)
        .rename(columns=dict(enumerate(COLUMNS)))
        .drop(["t", "date", "goXRefs", "pathwayXRefs"], axis=1)
    )

    # Load the gaps data
    with open(args.gaps, "r", encoding="utf-8") as f:
        gaps = json.load(f)

    print(len(set(gaps.keys()).intersection(set(doms.id))))

    coord_fn = partial(get_aligned_coords, gaps=gaps)

    doms[["start_aligned", "end_aligned"]] = doms.apply(coord_fn, axis=1).tolist()

    # Remove duplicate entry_accessions and save
    doms.sort_values(by=["id", "start", "signature_library"]).groupby(
        ["id", "entry_accession"]
    ).head(1).to_csv(args.output, sep="\t", index=False)
