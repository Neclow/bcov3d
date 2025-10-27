"""Download and prepare metadata and protein structure files from the RCSB PDB."""

import json
import multiprocessing
import os

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from functools import partial
from pathlib import Path

import pandas as pd

from biotite.database import rcsb
from tqdm import tqdm


def parse_args():
    parser = ArgumentParser(
        description="Download PDB files in CIF format.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file", help="Metadata (.csv) file containing files to download."
    )
    parser.add_argument("variant_file", help="Variants (.json)")
    parser.add_argument("output_dir", help="Directory to save the downloaded files.")
    parser.add_argument(
        "--col",
        default="PDB",
        help="Column name in the input file containing PDB IDs.",
    )
    parser.add_argument(
        "--format", default="cif", help="Format of the downloaded files."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads to use for downloading files.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Read metadata
    start_df = pd.read_csv(args.input_file).set_index(args.col)

    # Read list of major variants from the JSON file
    with open(args.variant_file, "r", encoding="utf-8") as f:
        major_variants = json.load(f)

    print(f"# Input entries: {len(start_df)}")

    # Filter out rows with undefined variants
    start_df["MajorVariant"] = start_df["Variant"].apply(
        lambda x: major_variants.get(str(x).split(" ", maxsplit=1)[0].lower(), "Other")
    )

    # Push clean metadata
    clean_df = start_df.query("MajorVariant != 'NS' & Domain == 'full'").sort_index()
    print(f"# Filtered entries: {len(clean_df)}")
    clean_df.to_csv(Path(args.input_file).parent / "clean.csv")

    # Set a name for each entry and dump mapping to a JSON file
    clean_df["name"] = (
        clean_df["Virus"]
        + "/"
        + clean_df["MajorVariant"]
        + "/"
        + clean_df.index.astype(str)
    )

    with open(
        Path(args.input_file).parent / "names.json",
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(clean_df.name.to_dict(), f, indent=4)

    # Download PDB files
    os.makedirs(args.output_dir, exist_ok=True)
    with multiprocessing.Pool(processes=args.threads) as pool:
        ids = list(clean_df.index)
        download_fn = partial(
            rcsb.fetch, format=args.format, target_path=args.output_dir, verbose=False
        )

        work = pool.imap_unordered(download_fn, ids)

        for output in list(tqdm(work, total=len(ids), desc="Downloading CIF files")):
            pass
