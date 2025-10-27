import bisect
import json
import operator
from argparse import ArgumentParser
from functools import partial

import pandas as pd
from Bio import SeqIO

fadir_aa = "data/aa/fa"
fadir_3di = "data/3di/fa"

columns = [
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


def get_aligned_coords(x, start_col, end_col, gaps, op="+"):
    gap_coords = gaps[x.id]
    if callable(op):
        op_fn = op
    elif op == "+":
        op_fn = operator.add
    elif op == "-":
        op_fn = operator.sub
    else:
        raise ValueError("op must be '+' or '-' or a callable")
    aligned_start = op_fn(x[start_col], sum(1 for g in gap_coords if g <= x[start_col]))
    aligned_end = op_fn(x[end_col], sum(1 for g in gap_coords if g <= x[end_col]))
    return [aligned_start, aligned_end]


def get_aligned_coords_from_old(x, col, gaps):
    guess = x[col]
    while True:
        count_gaps = bisect.bisect_right(gaps[x.id], guess)
        new_guess = x[col] + count_gaps
        if new_guess == guess:
            return guess
        guess = new_guess


def old_interval_to_new_coords(x, start_col, end_col, old_positions_in_new):
    # Find all NEW indices where OLD position falls inside [start_old, end_old)
    relevant_new_indices = [
        i
        for i, old_pos in enumerate(old_positions_in_new)
        if x[start_col] <= old_pos < x[end_col]
    ]
    if not relevant_new_indices:
        return None, None  # no match

    new_start = min(relevant_new_indices)
    new_end = max(relevant_new_indices) + 1  # if end-exclusive

    return new_start, new_end


def parse_args():
    parser = ArgumentParser(description="Clean domain data")
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input file containing domain data.",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to save the cleaned domain data.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    doms = (
        pd.read_csv(args.input, sep="\t", header=None)
        .rename(columns=dict(enumerate(columns)))
        .drop(["t", "date", "goXRefs", "pathwayXRefs"], axis=1)
    )

    # coord_fn = partial(get_aligned_coords, gaps=gaps)
    with open(f"{fadir_aa}/subset1_gaps.json", "r", encoding="utf-8") as f:
        gaps = json.load(f)
    doms[["start_subset", "end_subset"]] = doms.apply(
        partial(
            get_aligned_coords, start_col="start", end_col="end", gaps=gaps, op="+"
        ),
        axis=1,
    ).tolist()

    with open(f"{fadir_aa}/subset1.colnumbering", "r", encoding="utf-8") as f:
        subset_colnumbering = [
            int(_) for _ in f.read().split("#ColumnsMap ")[1].strip().split(", ")
        ]
    doms["start_filtersites"] = doms["start_subset"].apply(
        lambda x: subset_colnumbering[x - 1]
    )
    doms["end_filtersites"] = doms["end_subset"].apply(
        lambda x: subset_colnumbering[x - 1]
    )

    with open(f"{fadir_aa}/trimmed.colnumbering", "r", encoding="utf-8") as f:
        trimmed_colnumbering = [
            int(_) for _ in f.read().split("#ColumnsMap")[1].strip().split(", ")
        ]
    doms["start_trimmed"] = doms["start_filtersites"].apply(
        lambda x: trimmed_colnumbering[x - 1]
    )
    doms["end_trimmed"] = doms["end_filtersites"].apply(
        lambda x: trimmed_colnumbering[x - 1]
    )

    with open(f"{fadir_aa}/subset1.fa", "r", encoding="utf-8") as f:
        subset1 = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    with open(f"{fadir_aa}/aligned.fa", "r", encoding="utf-8") as f:
        aligned1 = {
            r.id: str(r.seq) for r in SeqIO.parse(f, "fasta") if r.id in subset1
        }
    aligned_gaps = {}
    for k in subset1.keys():
        aligned_gaps[k] = [i for i, c in enumerate(aligned1[k]) if c == "-"]
    doms[["start_clean", "end_clean"]] = doms.apply(
        partial(
            get_aligned_coords,
            start_col="start_trimmed",
            end_col="end_trimmed",
            gaps=aligned_gaps,
            op="-",
        ),
        axis=1,
    ).tolist()

    with open(f"{fadir_aa}/cleandiff.json", "r", encoding="utf-8") as f:
        cleandiff = json.load(f)
    doms[["start_3di", "end_3di"]] = doms.apply(
        partial(
            get_aligned_coords,
            start_col="start_clean",
            end_col="end_clean",
            gaps=cleandiff,
            op="-",
        ),
        axis=1,
    ).tolist()

    with open(f"{fadir_3di}/aligned.fa", "r", encoding="utf-8") as f:
        aligned1_3di = {
            r.id: str(r.seq) for r in SeqIO.parse(f, "fasta") if r.id in subset1
        }
    aligned_gaps_3di = {}
    for k in subset1.keys():
        aligned_gaps_3di[k] = [i for i, c in enumerate(aligned1_3di[k]) if c == "-"]
    doms["start_aligned_3di"] = doms.apply(
        partial(get_aligned_coords_from_old, col="start_3di", gaps=aligned_gaps_3di),
        axis=1,
    )
    doms["end_aligned_3di"] = doms.apply(
        partial(get_aligned_coords_from_old, col="end_3di", gaps=aligned_gaps_3di),
        axis=1,
    )

    with open(f"{fadir_3di}/trimmed.colnumbering", "r", encoding="utf-8") as f:
        trimmed_colnumbering_3di = [
            int(_) for _ in f.read().split("#ColumnsMap")[1].strip().split(", ")
        ]
    doms[["start_trimmed_3di", "end_trimmed_3di"]] = doms.apply(
        partial(
            old_interval_to_new_coords,
            start_col="start_aligned_3di",
            end_col="end_aligned_3di",
            old_positions_in_new=trimmed_colnumbering_3di,
        ),
        axis=1,
    ).tolist()

    with open(f"{fadir_3di}/subset1.colnumbering", "r", encoding="utf-8") as f:
        subset_colnumbering_3di = [
            int(_) for _ in f.read().split("#ColumnsMap ")[1].strip().split(", ")
        ]
    doms[["start_subset_3di", "end_subset_3di"]] = doms.apply(
        partial(
            old_interval_to_new_coords,
            start_col="start_trimmed_3di",
            end_col="end_trimmed_3di",
            old_positions_in_new=subset_colnumbering_3di,
        ),
        axis=1,
    ).tolist()

    doms["start_subset_3di"] += 1
    doms["end_subset_3di"] += 1

    example = "8i3w"
    example_start = 600
    sub_doms_aa = (
        doms.sort_values(by="start")
        .query("id == @example & start_subset > @example_start")
        .iloc[:, 11:]
        .head()
    )
    print(sub_doms_aa)

    example_doms = sub_doms_aa.iloc[0]
    with open(f"{fadir_aa}/clean.fa", "r", encoding="utf-8") as f:
        clean1 = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta") if r.id in subset1}
    with open(f"{fadir_aa}/unique.fa", "r", encoding="utf-8") as f:
        unique1 = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta") if r.id in subset1}

    for dataset, col in [
        (aligned1, "trimmed"),
        (clean1, "clean"),
        (subset1, "subset"),
        (unique1, "filtersites"),
    ]:
        print(
            dataset[example][
                int(example_doms[f"start_{col}"]) : int(example_doms[f"end_{col}"])
            ]
        )

    sub_doms_3di = (
        doms.sort_values(by="start")
        .query("id == @example & start_subset_3di > @example_start")
        .iloc[:, 19:]
        .head()
    )
    print(sub_doms_3di)

    example_doms_3di = sub_doms_3di.iloc[0]
    with open(f"{fadir_3di}/subset1.fa", "r", encoding="utf-8") as f:
        subset1_3di = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    with open(f"{fadir_3di}/clean.fa", "r", encoding="utf-8") as f:
        clean1_3di = {
            r.id: str(r.seq) for r in SeqIO.parse(f, "fasta") if r.id in subset1
        }

    with open(f"{fadir_3di}/unique.fa", "r", encoding="utf-8") as f:
        unique1_3di = {
            r.id: str(r.seq) for r in SeqIO.parse(f, "fasta") if r.id in subset1
        }
    for dataset, col in [
        (clean1_3di, "3di"),
        (aligned1_3di, "aligned_3di"),
        (unique1_3di, "trimmed_3di"),
        (subset1_3di, "subset_3di"),
    ]:
        print(
            dataset[example][
                int(example_doms_3di[f"start_{col}"]) : int(
                    example_doms_3di[f"end_{col}"]
                )
            ]
        )

    doms.sort_values(by=["id", "start", "signature_library"]).groupby(
        ["id", "entry_accession"]
    ).head(1).to_csv(args.output, sep="\t", index=False)

    print("Done")
