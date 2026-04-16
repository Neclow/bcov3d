"""Convert InterProScan domain coordinates to AA and 3DI subset alignment space.

Maps raw SEQRES coordinates from InterProScan through:
  1. SEQRES → foldseek AA (pairwise alignment) → raw 3DI position
  2. Raw 3DI position → 3DI subset alignment column (non-gap counting)
  3. Raw SEQRES position → AA subset alignment column (non-gap counting)
"""

import json

from argparse import ArgumentParser
from functools import partial

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner


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


def build_seqres_to_foldseek_map(seqres_seq, foldseek_seq):
    """Pairwise align SEQRES and foldseek AA to map SEQRES positions to foldseek positions.

    Returns a dict: seqres_pos (0-based) -> foldseek_pos (0-based), or None if gapped.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(seqres_seq, foldseek_seq)[0]

    seqres_pos = -1
    foldseek_pos = -1
    mapping = {}

    for s_char, f_char in zip(*alignment):
        if s_char != "-":
            seqres_pos += 1
        if f_char != "-":
            foldseek_pos += 1
        if s_char != "-" and f_char != "-":
            mapping[seqres_pos] = foldseek_pos

    return mapping


def build_raw_to_aligned_map(aligned_seq):
    """Map raw (ungapped) position to alignment column.

    Returns a dict: raw_pos (0-based) -> aligned_col (0-based).
    """
    mapping = {}
    raw_pos = 0
    for col, c in enumerate(aligned_seq):
        if c != "-":
            mapping[raw_pos] = col
            raw_pos += 1
    return mapping


def map_coord(pos, mapping):
    """Map a 1-based InterProScan coordinate through a 0-based mapping, return 1-based."""
    mapped = mapping.get(pos - 1)
    return mapped + 1 if mapped is not None else None


def parse_args():
    parser = ArgumentParser(
        description="Convert InterProScan domain coords to AA and 3DI subset alignment space."
    )
    parser.add_argument("interproscan_tsv", help="InterProScan output TSV.")
    parser.add_argument("seqres_fa", help="SEQRES AA FASTA (filtered, unaligned).")
    parser.add_argument("foldseek_aa_fa", help="Foldseek AA FASTA (resolved residues).")
    parser.add_argument("chains_json", help="Chain map JSON from extract_aa.")
    parser.add_argument("aa_subset_fa", help="AA subset alignment FASTA.")
    parser.add_argument("tdi_subset_fa", help="3DI subset alignment FASTA.")
    parser.add_argument("output", help="Output TSV with mapped coordinates.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Load InterProScan results
    doms = (
        pd.read_csv(args.interproscan_tsv, sep="\t", header=None)
        .rename(columns=dict(enumerate(COLUMNS)))
        .drop(["t", "date", "goXRefs", "pathwayXRefs"], axis=1)
    )

    # Load sequences
    seqres = {r.id: str(r.seq) for r in SeqIO.parse(args.seqres_fa, "fasta")}

    with open(args.chains_json, encoding="utf-8") as f:
        chains = json.load(f)

    foldseek_aa = {}
    for r in SeqIO.parse(args.foldseek_aa_fa, "fasta"):
        parts = r.id.split("_", 1)
        if len(parts) == 2:
            pdb = parts[0]
            chain = parts[1].split()[0]
            foldseek_aa.setdefault(pdb, {})[chain] = str(r.seq)

    aa_subset = {r.id: str(r.seq) for r in SeqIO.parse(args.aa_subset_fa, "fasta")}
    tdi_subset = {r.id: str(r.seq) for r in SeqIO.parse(args.tdi_subset_fa, "fasta")}

    subset_pdbs = set(aa_subset.keys())

    # Filter domains to subset PDBs only
    doms = doms[doms["id"].isin(subset_pdbs)].copy()
    print(f"Domains for {doms['id'].nunique()} subset structures: {len(doms)} entries")

    # Build mappings for each subset structure
    seqres_to_fs = {}  # PDB -> {seqres_pos -> foldseek_pos}
    fs_to_tdi_aligned = {}  # PDB -> {foldseek_pos -> 3DI alignment col}
    seqres_to_aa_aligned = {}  # PDB -> {seqres_pos -> AA alignment col}

    for pdb in sorted(subset_pdbs):
        chain = chains.get(pdb)

        # SEQRES -> foldseek AA mapping
        if pdb in foldseek_aa and chain in foldseek_aa[pdb]:
            seqres_to_fs[pdb] = build_seqres_to_foldseek_map(
                seqres[pdb], foldseek_aa[pdb][chain]
            )
        else:
            print(f"  WARNING: {pdb} chain {chain} not in foldseek AA, skipping")
            continue

        # Foldseek AA pos -> 3DI subset alignment column
        if pdb in tdi_subset:
            fs_to_tdi_aligned[pdb] = build_raw_to_aligned_map(tdi_subset[pdb])

        # SEQRES pos -> AA subset alignment column
        if pdb in aa_subset:
            seqres_to_aa_aligned[pdb] = build_raw_to_aligned_map(aa_subset[pdb])

        print(
            f"  {pdb}: SEQRES={len(seqres[pdb])}, "
            f"foldseek={len(foldseek_aa[pdb][chain])}, "
            f"mapped={len(seqres_to_fs[pdb])}"
        )

    # Map coordinates
    def map_to_aa_subset(row):
        m = seqres_to_aa_aligned.get(row.id)
        if m is None:
            return None, None
        return map_coord(row.start, m), map_coord(row.end, m)

    def map_to_tdi_subset(row):
        s2f = seqres_to_fs.get(row.id)
        f2t = fs_to_tdi_aligned.get(row.id)
        if s2f is None or f2t is None:
            return None, None
        fs_start = s2f.get(row.start - 1)
        fs_end = s2f.get(row.end - 1)
        if fs_start is None or fs_end is None:
            return None, None
        tdi_start = f2t.get(fs_start)
        tdi_end = f2t.get(fs_end)
        if tdi_start is None or tdi_end is None:
            return None, None
        return tdi_start + 1, tdi_end + 1

    doms[["start_aa_subset", "end_aa_subset"]] = (
        doms.apply(map_to_aa_subset, axis=1).tolist()
    )
    doms[["start_3di_subset", "end_3di_subset"]] = (
        doms.apply(map_to_tdi_subset, axis=1).tolist()
    )

    # Deduplicate
    doms = (
        doms.sort_values(by=["id", "start", "signature_library"])
        .groupby(["id", "entry_accession"])
        .head(1)
    )

    doms.to_csv(args.output, sep="\t", index=False)

    # Summary
    n_mapped_aa = doms["start_aa_subset"].notna().sum()
    n_mapped_3di = doms["start_3di_subset"].notna().sum()
    print(f"\nOutput: {len(doms)} domains, {n_mapped_aa} mapped to AA, {n_mapped_3di} mapped to 3DI")
