"""Clean and filter extracted sequences.

Performs three operations in order:
1. Match raw 3DI sequences to spike chains identified by extract_aa (chain map).
2. Filter metadata to spike-only, non-chimeric structures.
3. Filter both AA and 3DI FASTAs to only include structures in filtered metadata.
"""

import json

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import pandas as pd

from Bio import SeqIO


def parse_args():
    parser = ArgumentParser(
        description="Clean and filter extracted AA/3DI sequences and metadata.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    # Inputs
    parser.add_argument("raw_3di", help="Raw 3DI FASTA from extract_3di.")
    parser.add_argument("raw_aa", help="Raw AA FASTA from extract_aa.")
    parser.add_argument("chains", help="Chain map JSON from extract_aa.")
    parser.add_argument("metadata", help="Cleaned metadata CSV from download_data.")

    # Outputs
    parser.add_argument("out_3di", help="Output filtered 3DI FASTA.")
    parser.add_argument("out_aa", help="Output filtered AA FASTA.")
    parser.add_argument("out_metadata", help="Output filtered metadata CSV.")

    # Options
    parser.add_argument(
        "--spike_only",
        action="store_true",
        help="Keep only spike-only structures (no antibody/receptor).",
    )
    parser.add_argument(
        "--exclude_chimeric",
        action="store_true",
        help="Exclude chimeric constructs (requires 'Chimeric' column in metadata).",
    )
    parser.add_argument(
        "--min_3di_len",
        type=int,
        default=800,
        help="Minimum 3DI sequence length. Shorter = incomplete foldseek output.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Load inputs
    meta = pd.read_csv(args.metadata)
    n_start = len(meta)

    with open(args.chains, encoding="utf-8") as f:
        chain_map = json.load(f)

    # --- Step 1: Match 3DI to chain map ---
    raw_3di = SeqIO.to_dict(SeqIO.parse(args.raw_3di, "fasta"))
    print(f"Raw 3DI records: {len(raw_3di)}")

    matched_3di = {}
    skipped_3di = []
    for pdb_id, chain in chain_map.items():
        key = f"{pdb_id}_{chain}"
        if key in raw_3di:
            matched_3di[pdb_id] = str(raw_3di[key].seq)
        else:
            skipped_3di.append(pdb_id)

    print(f"3DI matched to chain map: {len(matched_3di)} (skipped {len(skipped_3di)})")
    if skipped_3di:
        print(f"  Missing: {skipped_3di}")

    # Filter out structures with too-short 3DI (incomplete foldseek output)
    short_3di = {k for k, v in matched_3di.items() if len(v) < args.min_3di_len}
    if short_3di:
        print(f"3DI too short (< {args.min_3di_len}): {len(short_3di)} removed")
        for pdb_id in sorted(short_3di):
            print(f"  {pdb_id}: {len(matched_3di[pdb_id])} residues")
        matched_3di = {k: v for k, v in matched_3di.items() if k not in short_3di}

    # --- Step 2: Filter metadata ---
    pdb_col = "PDB" if "PDB" in meta.columns else meta.columns[0]

    if args.spike_only:
        meta = meta[meta["Antibody"].isna() & meta["Receptor/Other"].isna()]
        print(f"Spike-only: {n_start} -> {len(meta)}")

    if args.exclude_chimeric and "Chimeric" in meta.columns:
        n_before = len(meta)
        meta = meta[~meta["Chimeric"]]
        print(f"Exclude chimeric: {n_before} -> {len(meta)}")

    # --- Step 3: Intersect metadata with available sequences ---
    raw_aa = SeqIO.to_dict(SeqIO.parse(args.raw_aa, "fasta"))
    print(f"Raw AA records: {len(raw_aa)}")

    # Only keep structures that have both AA and 3DI sequences
    aa_pdbs = set(raw_aa.keys())
    tdi_pdbs = set(matched_3di.keys())
    meta_pdbs = set(meta[pdb_col].str.lower())

    common_pdbs = aa_pdbs & tdi_pdbs & meta_pdbs
    print(f"Structures with AA + 3DI + metadata: {len(common_pdbs)}")

    # Filter metadata to common set
    meta = meta[meta[pdb_col].str.lower().isin(common_pdbs)]

    # Write filtered FASTAs
    n_aa = 0
    with open(args.out_aa, "w", encoding="utf-8") as f:
        for pdb_id in sorted(common_pdbs):
            n_aa += 1
            f.write(f">{pdb_id}\n{raw_aa[pdb_id].seq}\n")

    n_3di = 0
    with open(args.out_3di, "w", encoding="utf-8") as f:
        for pdb_id in sorted(common_pdbs):
            n_3di += 1
            f.write(f">{pdb_id}\n{matched_3di[pdb_id]}\n")

    meta.to_csv(args.out_metadata, index=False)

    n_viruses = meta["Virus"].nunique()
    print(f"\nFinal: {n_aa} AA, {n_3di} 3DI, {len(meta)} metadata entries, {n_viruses} viruses")


if __name__ == "__main__":
    main()
