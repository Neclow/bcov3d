import json
import os
import warnings

from argparse import ArgumentParser
from itertools import chain
from pathlib import Path

from ete4.parser.newick import NewickError
from ete4 import Tree
from tqdm import tqdm


def rename(intree, mapping):
    for leaf in intree.leaves():
        # Ignore node annotations in mapping
        if "[" in leaf.name:
            leaf_name, rest = leaf.name.split("[", 2)
            leaf.name = mapping[leaf_name] + "[" + rest
        else:
            leaf.name = mapping[leaf.name]


def parse_args():
    parser = ArgumentParser(description="Annotate trees with metadata.")
    parser.add_argument(
        "indir",
        type=str,
        help="Directory containing tree files to annotate.",
    )
    parser.add_argument(
        "mapper",
        type=str,
        help="Path to the JSON file containing the mapping of names.",
    )
    parser.add_argument(
        "--outgroup",
        type=str,
        help="Taxon name to use as outgroup for rooting trees. If not provided, trees will not be rooted.",
    )
    parser.add_argument(
        "--errors",
        choices=("ignore", "warn", "raise"),
        default="raise",
        help="How to handle errors during annotation.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing annotated files.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove all existing annotated files in the directory before processing.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    files = chain(
        Path(args.indir).rglob("*.treefile"),
        Path(args.indir).rglob("*.contree"),
    )

    with open(args.mapper, "r", encoding="utf-8") as f:
        mapping = json.load(f)

    if args.clean:
        # Remove all *.annotated files in the directory
        old_annotated_files = chain(
            Path(args.indir).rglob("*.treefile.annotated"),
            Path(args.indir).rglob("*.contree.annotated"),
        )
        for oaf in old_annotated_files:
            os.remove(oaf)

    for file in tqdm(files):
        if file.suffix == ".annotated":
            continue
        outfile = f"{file}.annotated"
        if os.path.exists(outfile) and not args.overwrite:
            if args.errors != "ignore":
                warnings.warn(
                    f"Skipping {file}, already annotated. Use --overwrite to force overwrite.",
                    UserWarning,
                )
            continue
        with open(file, "r", encoding="utf-8") as f_in, open(
            outfile, "w", encoding="utf-8"
        ) as f_out:
            for line in f_in.readlines():
                try:
                    tr = Tree(line.strip(), parser=1)
                    if args.outgroup is not None and args.outgroup in tr.leaf_names():
                        tr.set_outgroup(args.outgroup)
                    rename(tr, mapping)
                    f_out.write(tr.write(parser=1) + "\n")
                except (KeyError, NewickError) as e:
                    if args.errors == "raise":
                        raise ValueError(f"Error found with {file}") from e
                    if args.errors == "warn":
                        warnings.warn(
                            f"Error renaming tree in {file}: {e}. Skipping file",
                            UserWarning,
                        )
                    elif args.errors == "ignore":
                        # Do nothing
                        pass
                    continue
            print(f"Annotated {file} to {outfile}")
