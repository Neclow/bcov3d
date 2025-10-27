from argparse import ArgumentParser

import numpy as np
import pandas as pd


def parse_args():
    parser = ArgumentParser(description="Sort MAST results by tree weight.")
    parser.add_argument("-wslm", type=str, help="MAST .sitelh file", required=True)
    parser.add_argument(
        "-te", "--trees", type=str, help="Input .trees file", required=True
    )
    parser.add_argument(
        "--ascending",
        action="store_true",
        help="Sort in ascending order (default: descending)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    sitelhs = pd.read_csv(args.wslm, comment="#", sep="\t").set_index("Site")

    # Exclude first column (weighted log-likelihood)
    treelhs = sitelhs.sum(0)[1:]

    print(
        treelhs.sort_values(ascending=args.ascending)
        .to_frame(name="Tree log-likelihood")
        .to_markdown()
    )

    idxs = np.argsort(treelhs.values)
    if not args.ascending:
        idxs = idxs[::-1]

    print(f"Sorting trees from {args.trees}")

    with open(args.trees, "r", encoding="utf-8") as f_in:
        lines = f_in.readlines()

    with open(args.trees + ".sorted", "w", encoding="utf-8") as f_out:
        for i in idxs:
            f_out.write(lines[i])

    # idxs = None

    # with open(args.sitelh, "r", encoding="utf-8") as f:
    #     for line in f:
    #         if "Tree weights: " in line:
    #             tree_weights = np.array(
    #                 [
    #                     float(x)
    #                     for x in line.split("Tree weights: ")[1].strip().split(", ")
    #                 ]
    #             )
    #             idxs = np.argsort(tree_weights)
    #             if not args.ascending:
    #                 idxs = idxs[::-1]
    #             break

    # assert idxs is not None, "No tree weights found in MAST .iqtree file."

    # sort_file(args.trees, args.trees + ".sorted", idxs)
