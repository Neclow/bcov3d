import math
import os

from argparse import ArgumentParser
from pathlib import Path

from src.fasta import count, cut, maxlen


def parse_args():
    parser = ArgumentParser(description="Cut a FASTA file into smaller files.")
    parser.add_argument("infile", type=str, help="Input FASTA file")
    parser.add_argument(
        "min_size", type=int, help="Minimum size of each alignment chunk"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    inpath = Path(args.infile)

    # Make subdirectory for cut files
    cutdir = inpath.parent / "cut" / inpath.stem
    os.makedirs(cutdir, exist_ok=True)

    # Find sequence length and maximum number of chunks
    seqlen = maxlen(inpath, is_aligned=True)
    ntax = count(inpath)
    # Number of branch lengths ~ 2 * # taxa
    # For MAST to work, # sites must be > # parameters
    # But chunk sizes must also be reasonable, so we tune
    # with args.min_size
    nparam = 2 * ntax
    # min_n = 2
    max_n = math.floor(seqlen / max(nparam, args.min_size))

    print(f"Max n: {max_n}, sequence length: {seqlen}")

    chunksize = seqlen // max_n
    print(chunksize)
    cut(inpath, cutdir, chunksize, min_size=args.min_size)
