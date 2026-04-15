"""Extract amino acid sequences from CIF files, identifying spike chains by length."""

import json
import multiprocessing
import warnings

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from functools import partial
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm


def parse_args():
    parser = ArgumentParser(
        description="Extract spike amino acid sequences from CIF files.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("cif_dir", help="Directory containing CIF files.")
    parser.add_argument("output_raw", help="Output FASTA file with amino acid sequences.")
    parser.add_argument(
        "output_chains",
        help="Output JSON mapping PDB ID to selected chain ID.",
    )
    parser.add_argument(
        "--min_len",
        type=int,
        default=1000,
        help="Minimum sequence length to consider as spike.",
    )
    parser.add_argument(
        "--metadata",
        help="Metadata CSV. If provided, only process PDBs listed in it.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads for processing.",
    )
    return parser.parse_args()


def extract_spike_chain(cif_path, min_len):
    """Extract the longest spike-length chain from a CIF file.

    Uses auth_asym_id (author chain ID) for the chain mapping, since foldseek
    uses auth_asym_id when generating 3DI sequences.

    Parameters
    ----------
    cif_path : Path
        Path to a CIF file.
    min_len : int
        Minimum length to consider as a spike chain.

    Returns
    -------
    tuple
        (pdb_id, auth_chain_id, sequence) or (pdb_id, None, None) if no spike chain found.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    pdb_id = cif_path.stem
    best_chain = None
    best_seq = None
    best_len = 0
    best_label_chain = None

    try:
        # First pass: find the longest spike chain using BioPython's seqres parser
        for record in SeqIO.parse(str(cif_path), "cif-seqres"):
            seq = str(record.seq)
            if len(seq) >= min_len and len(seq) > best_len:
                best_len = len(seq)
                best_label_chain = record.id.split(":")[-1] if ":" in record.id else record.id
                best_seq = seq

        if best_label_chain is None:
            return (pdb_id, None, None)

        # Second pass: map label_asym_id -> auth_asym_id
        d = MMCIF2Dict(str(cif_path))
        label_ids = d.get("_pdbx_poly_seq_scheme.asym_id", [])
        auth_ids = d.get("_pdbx_poly_seq_scheme.pdb_strand_id", [])

        label_to_auth = {}
        for label, auth in zip(label_ids, auth_ids):
            if label not in label_to_auth:
                label_to_auth[label] = auth

        best_chain = label_to_auth.get(best_label_chain, best_label_chain)

    except (KeyError, UnicodeDecodeError, ValueError) as err:
        warnings.warn(f"Failed to process {pdb_id}.cif: {err}")

    return (pdb_id, best_chain, best_seq)


if __name__ == "__main__":
    args = parse_args()

    cif_dir = Path(args.cif_dir)
    cifs = sorted(cif_dir.glob("*.cif"))

    if args.metadata:
        import pandas as pd

        meta = pd.read_csv(args.metadata)
        allowed_pdbs = set(meta["PDB"].str.lower())
        cifs = [c for c in cifs if c.stem.lower() in allowed_pdbs]
        print(f"Metadata filter: {len(cifs)} structures to process")

    extract_fn = partial(extract_spike_chain, min_len=args.min_len)

    with multiprocessing.Pool(processes=args.threads) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(extract_fn, cifs),
                total=len(cifs),
                desc="Extracting AA sequences",
            )
        )

    # Collect successful extractions
    aa_seqs = {}
    chain_map = {}
    for pdb_id, chain_id, seq in results:
        if seq is not None:
            aa_seqs[pdb_id] = seq
            chain_map[pdb_id] = chain_id

    print(f"Structures with spike chain found: {len(aa_seqs)} / {len(cifs)}")

    # Write FASTA
    with open(args.output_raw, "w", encoding="utf-8") as f_out:
        for pdb_id, seq in sorted(aa_seqs.items()):
            f_out.write(f">{pdb_id}\n{seq}\n")

    # Write chain mapping
    with open(args.output_chains, "w", encoding="utf-8") as f_out:
        json.dump(chain_map, f_out, indent=2, sort_keys=True)
