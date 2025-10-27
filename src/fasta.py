"""FASTA processing functions."""

import math
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd

from Bio import SeqIO
from scipy.stats import entropy


def count(input_file):
    """Count the number of sequences in a FASTA file.

    Parameters
    ----------
    input_file : str
        Input FASTA file

    Returns
    -------
    int
        Number of sequences in the FASTA file
    """
    return int(
        subprocess.run(
            f"grep -c '>' {input_file}",
            capture_output=True,
            text=True,
            check=True,
            shell=True,
        ).stdout
    )


def maxlen(infile, is_aligned=False):
    """Find the maximum sequence length in a FASTA file.

    If it is an alignment, it will simply parse the length
    of the first sequence.

    Parameters
    ----------
    infile : str
        Input FASTA file

    Returns
    -------
    int
        Maximum sequence length found in the file
    """
    records = SeqIO.parse(infile, "fasta")

    if is_aligned:
        return len(next(records).seq)

    return max(len(record.seq) for record in records)


def cut(infile, outidr, size, min_size=None):
    """Cut a FASTA file into smaller files of specified size.

    Parameters
    ----------
    infile : str
        Input FASTA file
    outstem : str
        Output directory
    size : int
        Number of sequences per output file
    min_size : int
        Minimum size of each alignment chunk
    """
    records = SeqIO.parse(infile, "fasta")

    chunked_records = defaultdict(dict)

    min_size = min_size or sys.maxsize

    for record in records:
        seq = str(record.seq)

        for j, i in enumerate(range(0, len(seq), size)):

            chunk = seq[i : i + size]
            if len(chunk) < min_size:
                break
            chunked_records[j][record.id] = chunk

    d = math.floor(len(chunked_records) / 10)
    for j in chunked_records:
        out_file = f"{outidr}/{j:0{d}}.fa"
        with open(out_file, "w", encoding="utf-8") as f_out:
            for id_, seq in chunked_records[j].items():
                f_out.write(f">{id_}\n{seq}\n")


def deduplicate(infile, outfile, verbose=False):
    """Remove duplicate sequences from a FASTA file.

    Saves a JSON file with key = unique sequences and values = duplicate IDs.

    Parameters
    ----------
    infile : str
        Input FASTA file
    outfile : str
        Output FASTA file
    verbose : bool, optional
        If True, print additional information, by default False
    """
    unique = {}

    for record in SeqIO.parse(infile, "fasta"):
        if str(record.seq) in unique:
            unique[str(record.seq)].append(str(record.id))
        else:
            unique[str(record.seq)] = [str(record.id)]

    if verbose:
        print(f"(deduplicate) Found {len(unique)} unique sequences")

    with open(
        f"{os.path.splitext(outfile)[0]}.json", "w", encoding="utf-8"
    ) as f_unique:
        json.dump({ids[0]: ids for ids in unique.values()}, f_unique, indent=2)

    with open(outfile, "w", encoding="utf-8") as f_out:
        for seq, ids in unique.items():
            f_out.write(f">{ids[0]}\n{seq}\n")


def dropna(infile, outfile, thresh):
    """Filter sequences in a FASTA file based on the proportion of gaps.

    Parameters
    ----------
    infile : str
        Input FASTA file
    outfile : str
        Output FASTA file with filtered sequences
    thresh : float
        Threshold for the proportion of gaps (0.0 to 1.0)
    """
    assert 0 <= thresh <= 1, "Threshold must be between 0 and 1"
    with open(outfile, "w", encoding="utf-8") as f_out:
        for record in SeqIO.parse(infile, "fasta"):
            pgap = record.seq.count("-") / len(record.seq)
            if pgap < thresh:
                f_out.write(f">{record.id}\n{str(record.seq)}\n")


def intersect(input1, input2, output1, output2):
    """Find common sequences between two FASTA files and write them to separate output files.

    Parameters
    ----------
    input1 : str
        First input FASTA file
    input2 : str
        Second input FASTA file
    output1 : str
        Output FASTA file for sequences from input1
    output2 : str
        Output FASTA file for sequences from input2

    Returns
    -------
    cross_keys : set
        Set of sequence IDs that are common to both input files
    """
    # Get unique sequences and intersection between 3di and aa
    with open(input1, "r", encoding="utf-8") as f:
        data1 = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    with open(input2, "r", encoding="utf-8") as f:
        data2 = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    cross_keys = sorted(set(data1.keys()).intersection(set(data2.keys())))

    print(f"(intersect) found {len(cross_keys)} common sequence IDs.")

    with (
        open(output1, "w", encoding="utf-8") as f1,
        open(output2, "w", encoding="utf-8") as f2,
    ):
        for id_ in cross_keys:
            f1.write(f">{id_}\n{str(data1[id_].seq)}\n")
            f2.write(f">{id_}\n{str(data2[id_].seq)}\n")

    return cross_keys


def pgap(infile):
    """Calculate the proportion of gaps in each sequence of a FASTA file.

    Parameters
    ----------
    infile : str
        Input FASTA file

    Returns
    -------
    dict
        Sequence IDs as index and proportion of gaps as values
    """
    pgaps = {}
    for record in SeqIO.parse(infile, "fasta"):
        seq = str(record.seq)
        pgaps[record.id] = seq.count("-") / len(seq)
    return pgaps


def filter_taxa(infile, outfile, ids):
    """Filter sequences in a FASTA file based on a list of IDs.

    Parameters
    ----------
    infile : str
        Input FASTA file
    outfile : str
        Output FASTA file with filtered sequences
    ids : list of str
        List of sequence IDs to keep
    """

    with open(outfile, "w", encoding="utf-8") as f_out:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in ids:
                f_out.write(f">{record.id}\n{str(record.seq)}\n")


def to_csv(infile, outfile):
    """Convert a FASTA file to a CSV file with sequence IDs and sequences.

    Parameters
    ----------
    infile : str
        Input FASTA file
    outfile : str
        Output CSV file
    """
    with open(outfile, "w", encoding="utf-8") as f_out:
        for record in SeqIO.parse(infile, "fasta"):
            f_out.write(f"{record.id},{','.join(str(record.seq))}\n")


def filter_sites(infile, outfile, colfile=None, max_unique=10):
    """
    Filter sites in a FASTA file based on
    the number of unique chars in a site.

    Parameters
    ----------
    infile : str
        Input FASTA file
    outfile : str
        Output FASTA file with filteredsites
    max_unique : int, optional
        Max. number of unique chars in a site, by default 10
    """
    # Read FASTA as csv
    with tempfile.NamedTemporaryFile() as fp:
        to_csv(infile, fp.name)
        df = pd.read_csv(fp.name, header=None, index_col=0)

    print(f"(subset/{infile}) # of sites before: {df.shape[1]}")

    # Find sites with OK number of unique chars
    unique_counts = df.apply(pd.Series.nunique)
    ok_counts = unique_counts[unique_counts < max_unique]
    ok_counts.value_counts().sort_index().to_json(
        f"{os.path.splitext(outfile)[0]}_okvaluecounts.json", indent=4
    )

    # Filter the OK sites
    sites_to_keep = ok_counts.index.tolist()
    new_df = df.loc[:, sites_to_keep]

    print(f"(subset/{infile}) # of sites after: {new_df.shape[1]}")
    with open(outfile, "w", encoding="utf-8") as f_out:
        for k, v in new_df.agg("".join, axis=1).items():
            f_out.write(f">{k}\n{v}\n")

    if colfile is not None:
        stem, _ = os.path.splitext(outfile)
        colfile = f"{stem}.colnumbering"
    with open(colfile, "w", encoding="utf-8") as f_col:
        f_col.write(f"#ColumnsMap {', '.join(map(str,sites_to_keep))} ")


def unalign(infile, outfile, gapfile):
    """Remove alignment from a FASTA file.

    Parameters
    ----------
    infile : str
        Input aligned FASTA file
    outfile : str
        Output unaligned FASTA file
    """
    gap_idxs = {}
    with (
        open(outfile, "w", encoding="utf-8") as f_out,
        open(gapfile, "w", encoding="utf-8") as f_gap,
    ):
        for record in SeqIO.parse(infile, "fasta"):
            gap_idxs[record.id] = [
                i for i, char in enumerate(str(record.seq)) if char == "-"
            ]
            f_out.write(f">{record.id}\n{str(record.seq).replace('-', '')}\n")
        json.dump(gap_idxs, f_gap, indent=4)


def entropy_str(col):
    counts = col.value_counts(normalize=True)
    return entropy(counts, base=2)


def site_entropy(seq_file, window):
    with tempfile.NamedTemporaryFile() as f:
        to_csv(seq_file, f.name)
        seq = pd.read_csv(f.name, header=None, index_col=0)

    ent = np.convolve(
        seq.apply(entropy_str, axis=0), np.ones(window) / window, mode="same"
    )

    return ent
