# pylint: disable=redefined-outer-name
"""Retrieve AA sequences from CIF files based on 3di FASTA input."""

import json
import multiprocessing
import warnings

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from functools import partial
from pathlib import Path

import pandas as pd

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import SeqIO
from tqdm import tqdm


def parse_args():
    parser = ArgumentParser(
        description="Retrieve amino acid sequences from 3di FASTA file.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_fa", help="Input 3di FASTA file.")
    parser.add_argument(
        "output_fa", help="Output FASTA file with amino acid sequences."
    )
    parser.add_argument(
        "output_diff",
        help="Output JSON file to note missing residues in the structure.",
    )
    parser.add_argument(
        "--min_len", type=int, default=0, help="Minimum length of sequences to keep."
    )
    parser.add_argument(
        "--chains",
        type=str,
        default="A",
        help="Comma-separated list of chains to keep.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads to use for processing.",
    )
    parser.add_argument(
        "--q",
        type=float,
        default=0.99,
        help="Quantile for the maximum length.",
    )
    return parser.parse_args()


def extract_from_homomultimer(id_, cif_dir, chains):
    """Process a single CIF file to extract amino acid sequences.

    Parameters
    ----------
    id_ : str
        Identifier for the CIF file (e.g., '6VXX').
    cif_dir : Path
        Directory where CIF files are stored.
    chains : list
        List of chains to consider for the amino acid sequence extraction.

    Returns
    -------
    tuple
        A tuple containing the identifier and the longest amino acid sequence found.
        If no valid sequence is found, returns an empty string and None.
    """
    aa_id_ = {}
    # We utilise the fact that spike is a homotrimer, so we can just take the longest chain
    max_len = 0
    max_id = ""
    try:
        for i, record in enumerate(SeqIO.parse(f"{cif_dir}/{id_}.cif", "cif-seqres")):
            # cif-atom yields tons of X, so not sure what to do with it
            # so using cif-seqres instead
            seq = str(record.seq)
            seq_id = f"{id_}_{chains[i]}"
            if len(seq) > max_len:
                max_id = seq_id
                max_len = len(seq)
            aa_id_[seq_id] = seq
            if i == 2:
                break
        return (max_id, aa_id_[max_id])  # Return the longest sequence
    except (KeyError, UnicodeDecodeError) as err:
        warnings.warn(
            f"Failed to process {id_}.cif: non-CIF format compatible entry. {err}"
        )
    return ("", None)


def diff_aa_to_struct(
    id_,
    cif_dir,
    aa_id_key="_pdbx_poly_seq_scheme.seq_id",
    structure_id_key="_pdbx_poly_seq_scheme.auth_seq_num",
):
    mmcif_dict = MMCIF2Dict(f"{cif_dir}/{id_}.cif")

    diff = []

    for v_aa, v_3di in zip(mmcif_dict[aa_id_key], mmcif_dict[structure_id_key]):
        if v_3di == "?":
            diff.append(int(v_aa))

    return (id_, list(set(diff)))


if __name__ == "__main__":
    args = parse_args()

    # Load sequence IDs
    unique_ids = {r.id for r in SeqIO.parse(args.input_fa, "fasta")}

    # CIF directory
    cif_dir = Path(args.input_fa).parents[2] / "cif"

    # Parse chains
    chains = args.chains.split(",")

    raw_cov3d_aa = {}
    with multiprocessing.Pool(processes=args.threads) as pool:
        extract_fn = partial(
            extract_from_homomultimer,
            cif_dir=cif_dir,
            chains=chains,
        )

        work = pool.imap_unordered(extract_fn, unique_ids)

        seqs = list(tqdm(work, total=len(unique_ids), desc="Processing CIF files"))

        for id_, seq in seqs:
            if seq:
                raw_cov3d_aa[id_.split("_")[0]] = seq

        diff_fn = partial(
            diff_aa_to_struct,
            cif_dir=cif_dir,
        )

        work2 = pool.imap_unordered(diff_fn, unique_ids)

        diffs = dict(tqdm(work2, total=len(unique_ids), desc="Getting sequence diffs"))

    len_aa_cov3d = {k: len(seq) for k, seq in raw_cov3d_aa.items()}
    max_len_aa = int(pd.Series(len_aa_cov3d).quantile(args.q))
    print(f"Max length of amino acid sequences: {max_len_aa}")

    clean_cov3d_aa = {
        k: v for k, v in raw_cov3d_aa.items() if args.min_len <= len(v) < max_len_aa
    }
    print(f"Number of raw sequences: {len(raw_cov3d_aa)}")
    print(f"Number of sequences with OK length: {len(clean_cov3d_aa)}")

    # Push clean amino acid sequences to output FASTA file
    with open(args.output_fa, "w", encoding="utf-8") as f_out:
        for id_, seq in clean_cov3d_aa.items():
            f_out.write(f">{id_}\n{seq}\n")

    # Save the differences to a JSON file
    with open(args.output_diff, "w", encoding="utf-8") as f_diff:
        json.dump(diffs, f_diff, indent=4)
