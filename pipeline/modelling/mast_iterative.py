# pylint: disable=redefined-outer-name
import os
import subprocess
import tempfile

from argparse import ArgumentParser
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd


def grep_metrics_from_iqtree(file):
    with open(file, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()

    i = -1
    for i, line in enumerate(lines):
        if "Akaike" in line:
            break

    if i > -1:
        metrics_dict = {
            "AIC": lines[i].split("score:")[1].strip(),
            "AICc": lines[i + 1].split("score:")[1].strip(),
            "BIC": lines[i + 2].split("score:")[1].strip(),
        }
    else:
        metrics_dict = {"AIC": float("nan"), "AICc": float("nan"), "BIC": float("nan")}

    return metrics_dict


def parse_args():
    parser = ArgumentParser(description="Iteratively run MAST on a set of FASTA files.")
    parser.add_argument("-s", type=str, help="Input sequence", required=True)
    parser.add_argument("-m", type=str, help="Substitution model", required=True)
    parser.add_argument(
        "-te", "--trees", type=str, help="Input .trees file", required=True
    )
    parser.add_argument(
        "--prefix_dir",
        type=str,
        help="Directory prefix for IQ-TREE output files",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for IQ-TREE (default: 42)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for IQ-TREE (default: 1)",
    )
    parser.add_argument(
        "--metric",
        type=str,
        choices=["AIC", "AICc", "BIC"],
        default="BIC",
        help="Metric to use for comparison.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    indir = args.prefix_dir

    with open(args.trees, "r", encoding="utf-8") as f:
        sorted_trees = f.readlines()

    for i in range(1, len(sorted_trees)):
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".trees"
        ) as temp_trees:
            for j in range(i):
                temp_trees.write(sorted_trees[j])

            print(
                f"(mast_iterative/{args.s}) Running on best {i}/{len(sorted_trees)} trees..."
            )

        command = [
            "iqtree",
            "-s",
            args.s,
            "-m",
            f"{args.m}+T" if i > 1 else args.m,
            "-te",
            temp_trees.name,
            "-wspm" if i > 1 else "",
            "-wslm" if i > 1 else "",
            "-seed",
            f"{args.seed}",
            "-nt",
            f"{args.threads}",
            "--prefix",
            f"{indir}/{i}",
            "--quiet",
            "--redo",
        ]

        subprocess.run(" ".join(command), check=True, shell=True)

    print("MAST runs done. Collecting metrics...")

    mastiq_files = sorted(glob(f"{indir}/*.iqtree"))

    all_metrics = {}

    for f in mastiq_files:
        all_metrics[Path(f).stem] = grep_metrics_from_iqtree(f)

    metrics_df = pd.DataFrame.from_dict(all_metrics, orient="index")
    metrics_df.sort_values(by=args.metric, inplace=True)
    metrics_df.to_csv(f"{indir}.csv")

    print("Done")
