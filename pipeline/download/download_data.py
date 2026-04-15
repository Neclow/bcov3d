"""Download and prepare metadata and protein structure files from the RCSB PDB."""

import gzip
import json
import multiprocessing
import os
import urllib.request
import xml.etree.ElementTree as ET

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from functools import partial

import pandas as pd

from biotite.database import rcsb
from tqdm import tqdm

VALIDATION_URL = (
    "https://files.wwpdb.org/pub/pdb/validation_reports/{mid}/{pdb}/{pdb}_validation.xml.gz"
)

VALIDATION_FIELDS = {
    "absolute-percentile-clashscore": "percentile_clashscore",
    "relative-percentile-clashscore": "percentile_clashscore_relative",
    "absolute-percentile-percent-rama-outliers": "percentile_rama",
    "relative-percentile-percent-rama-outliers": "percentile_rama_relative",
    "absolute-percentile-percent-rota-outliers": "percentile_rota",
    "relative-percentile-percent-rota-outliers": "percentile_rota_relative",
    "clashscore": "clashscore",
    "percent-rama-outliers": "rama_outliers",
    "percent-rota-outliers": "rota_outliers",
}


STABILISATION_TAGS = [
    "foldon", "fibritin", "t4 lysozyme", "t4 phage lysozyme", "gfp",
    "green fluorescent", "mcherry", "strep", "his-tag", "flag",
]


def check_chimeric(cif_path):
    """Check if a structure has a chimeric spike entity.

    Flags structures where the spike polymer entity description contains a
    comma, indicating a fusion construct. Excludes common stabilisation tags
    (foldon, fibritin, etc.) which are standard in spike ectodomain constructs.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    pdb_id = cif_path.stem
    try:
        d = MMCIF2Dict(str(cif_path))
        descriptions = d.get("_entity.pdbx_description", [])
        types = d.get("_entity.type", [])
        for desc, etype in zip(descriptions, types):
            if etype != "polymer" or "," not in desc:
                continue
            desc_lower = desc.lower()
            # Only flag if the spike entity itself is chimeric
            if "spike" not in desc_lower and "glycoprotein" not in desc_lower:
                continue
            # Check if the fusion partner is just a stabilisation tag
            parts = [p.strip().lower() for p in desc.split(",")]
            non_spike = [p for p in parts if "spike" not in p and "glycoprotein" not in p]
            if any(
                tag in partner for partner in non_spike for tag in STABILISATION_TAGS
            ):
                continue
            return pdb_id, True
    except Exception:
        pass
    return pdb_id, False


def fetch_validation(pdb_id):
    """Fetch wwPDB validation percentiles for a single PDB entry."""
    pdb = pdb_id.lower()
    mid = pdb[1:3]
    url = VALIDATION_URL.format(mid=mid, pdb=pdb)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            xml_bytes = gzip.decompress(resp.read())
        root = ET.fromstring(xml_bytes)
        entry = root.find("Entry")
        if entry is None:
            return pdb_id, {}
        return pdb_id, {
            col: float(entry.attrib[attr])
            for attr, col in VALIDATION_FIELDS.items()
            if attr in entry.attrib
        }
    except Exception:
        return pdb_id, {}


def parse_args():
    parser = ArgumentParser(
        description="Download PDB files in CIF format.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file", help="Metadata (.csv) file containing files to download."
    )
    parser.add_argument("variant_file", help="Variants (.json)")
    parser.add_argument("output_dir", help="Directory to save the downloaded CIF files.")
    parser.add_argument("output_metadata", help="Output path for cleaned metadata CSV.")
    parser.add_argument("output_names", help="Output path for names JSON.")
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

    # Filter to clean entries
    clean_df = start_df.query("MajorVariant != 'NS' & Domain == 'full'").sort_index()
    print(f"# Filtered entries: {len(clean_df)}")

    # Set a name for each entry and dump mapping to a JSON file
    clean_df["name"] = (
        clean_df["Virus"]
        + "/"
        + clean_df["MajorVariant"]
        + "/"
        + clean_df.index.astype(str)
    )

    with open(args.output_names, "w", encoding="utf-8") as f:
        json.dump(clean_df.name.to_dict(), f, indent=4)

    # Download PDB files
    os.makedirs(args.output_dir, exist_ok=True)
    ids = list(clean_df.index)
    with multiprocessing.Pool(processes=args.threads) as pool:
        download_fn = partial(
            rcsb.fetch, format=args.format, target_path=args.output_dir, verbose=False
        )
        for _ in tqdm(pool.imap_unordered(download_fn, ids), total=len(ids), desc="Downloading CIF files"):
            pass

    # Check for chimeric constructs
    from pathlib import Path

    print("Checking for chimeric constructs...")
    cif_paths = [Path(args.output_dir) / f"{pdb}.cif" for pdb in ids]
    cif_paths = [p for p in cif_paths if p.exists()]
    with multiprocessing.Pool(processes=args.threads) as pool:
        chimeric_results = list(
            tqdm(pool.imap_unordered(check_chimeric, cif_paths), total=len(cif_paths), desc="Chimeric check")
        )
    chimeric_set = {pdb_id for pdb_id, is_chimeric in chimeric_results if is_chimeric}
    clean_df["Chimeric"] = clean_df.index.isin(chimeric_set)
    n_chimeric = clean_df["Chimeric"].sum()
    print(f"Chimeric constructs: {n_chimeric}/{len(clean_df)}")

    # Fetch wwPDB validation reports (percentiles + quality metrics)
    print("Fetching wwPDB validation reports...")
    with multiprocessing.Pool(processes=args.threads) as pool:
        val_results = list(
            tqdm(pool.imap_unordered(fetch_validation, ids), total=len(ids), desc="Validation reports")
        )

    val_df = pd.DataFrame(
        {pdb_id: metrics for pdb_id, metrics in val_results if metrics}
    ).T
    if not val_df.empty:
        val_df.index.name = clean_df.index.name
        clean_df = clean_df.join(val_df)
        n_with_val = val_df.notna().any(axis=1).sum()
        print(f"Validation data retrieved for {n_with_val}/{len(ids)} structures")

    # Write final metadata with validation metrics
    clean_df.to_csv(args.output_metadata)
