"""Data cleaning on 3di files: filtering out sequences based on length and chains."""

from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm


def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(description="Clean 3di files.")
    parser.add_argument(
        "input_fa",
        help="Input 3di FASTA file. Should be in the format ><cif_id>_<chain>.",
    )
    parser.add_argument("output_fa", help="Cleaned 3di FASTA file.")
    parser.add_argument(
        "--min_len", type=int, default=0, help="Minimum length of sequences to keep."
    )
    parser.add_argument(
        "--max_len",
        type=int,
        default=10000,
        help="Maximum length of sequences to keep.",
    )
    parser.add_argument(
        "--chains",
        type=str,
        default="A",
        help="Comma-separated list of chains to keep.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Load 3di fasta file
    records = SeqIO.to_dict(SeqIO.parse(args.input_fa, "fasta"))
    print("# Records loaded:", len(records))

    # Find related CIF files
    cif_dir = Path(args.input_fa).parents[2] / "cif"
    cifs = cif_dir.glob("*.cif")

    # Parse chains
    chains = args.chains.split(",")

    count = 0
    with open(args.output_fa, "w", encoding="utf-8") as f_out:
        for cif in tqdm(cifs):
            cif_id = Path(cif).stem

            cif_3di_seqs = []

            # Find all chain sequences for the current CIF ID
            for chain in chains:
                if f"{cif_id}_{chain}" in records:
                    cif_3di_seqs.append(str(records[f"{cif_id}_{chain}"].seq))

            # If we have sequences for all specified chains, continue processing
            if len(cif_3di_seqs) == len(chains):
                lens = {len(seq) for seq in cif_3di_seqs}

                if min(lens) > args.max_len or max(lens) < args.min_len:
                    # Ignore sequences that are all too short or too long
                    continue

                count += 1
                # Otherwise, write the sequences to the output file
                if len(lens) == 1:
                    # All sequences are the same length
                    f_out.write(f">{cif_id}\n{cif_3di_seqs[0]}\n")
                else:
                    # Sequences are of different lengths, keep the longest one
                    longest_seq = max(cif_3di_seqs, key=len)
                    f_out.write(f">{cif_id}\n{longest_seq}\n")

    print("# Records written:", count)
