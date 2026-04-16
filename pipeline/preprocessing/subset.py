"""Select one representative structure per virus/variant group.

Uses a geometric mean rank across resolution, 3DI completeness (gap fraction),
and wwPDB model quality percentiles to pick the best representative.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import pandas as pd
from Bio import SeqIO
from scipy.stats import gmean


def parse_args():
    parser = ArgumentParser(
        description="Select best representative per virus/variant group.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("filtered_3di", help="Filtered 3DI FASTA (trimmed).")
    parser.add_argument("filtered_aa", help="Filtered AA FASTA (trimmed).")
    parser.add_argument("metadata", help="Filtered metadata CSV.")

    parser.add_argument("out_3di", help="Output subset 3DI FASTA.")
    parser.add_argument("out_aa", help="Output subset AA FASTA.")
    parser.add_argument("out_metadata", help="Output subset metadata CSV.")
    parser.add_argument("out_ranked", help="Output full ranked metadata CSV (all candidates).")

    parser.add_argument(
        "--by",
        default="Virus,MajorVariant",
        help="Comma-separated columns to group by for selection.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    group_cols = [c.strip() for c in args.by.split(",")]

    meta = pd.read_csv(args.metadata)

    # Gap fraction in trimmed 3DI (fewer gaps = more complete)
    tdi_gaps = {}
    for r in SeqIO.parse(args.filtered_3di, "fasta"):
        seq = str(r.seq)
        tdi_gaps[r.id] = seq.count("-") / len(seq) if len(seq) > 0 else 1.0
    meta["3di_gap_frac"] = meta["PDB"].map(tdi_gaps)

    # Geometric mean rank across quality metrics (lower = better for all)
    # NaN clashscores (wwPDB software issue) are skipped, not penalised
    quality_cols = ["Resolution", "clashscore", "rama_outliers", "rota_outliers", "3di_gap_frac"]
    ranks = meta[quality_cols].rank()
    meta["geo_rank"] = ranks.apply(lambda r: gmean(r.dropna()), axis=1)

    # Save full ranked metadata
    meta.sort_values("geo_rank").to_csv(args.out_ranked, index=False)

    # Select best (lowest geo_rank) per group
    best = meta.sort_values("geo_rank").groupby(group_cols).first().reset_index()
    selected_pdbs = set(best["PDB"])

    print(f"Selected: {len(best)} structures from {best['Virus'].nunique()} viruses")
    print(f"Grouped by: {group_cols}")

    # Write outputs
    aa_records = SeqIO.to_dict(SeqIO.parse(args.filtered_aa, "fasta"))
    tdi_records = SeqIO.to_dict(SeqIO.parse(args.filtered_3di, "fasta"))

    with open(args.out_aa, "w", encoding="utf-8") as f:
        for pdb in sorted(selected_pdbs):
            if pdb in aa_records:
                f.write(f">{pdb}\n{aa_records[pdb].seq}\n")

    with open(args.out_3di, "w", encoding="utf-8") as f:
        for pdb in sorted(selected_pdbs):
            if pdb in tdi_records:
                f.write(f">{pdb}\n{tdi_records[pdb].seq}\n")

    best.to_csv(args.out_metadata, index=False)

    # Summary
    for _, row in best.sort_values("geo_rank").iterrows():
        group = "/".join(str(row[c]) for c in group_cols)
        print(
            f"  {row['PDB']}  {group:40s}  "
            f"res={row['Resolution']:.2f}  "
            f"clash={row['clashscore']:.1f}  "
            f"rama={row['rama_outliers']:.2f}  "
            f"rota={row['rota_outliers']:.2f}  "
            f"gaps={row['3di_gap_frac']:.1%}  "
            f"rank={row['geo_rank']:.1f}"
        )


if __name__ == "__main__":
    main()
