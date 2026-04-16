#!/usr/bin/env bash
# Extract 3DI sequences from CIF files using Foldseek.
set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") <cif_dir> <foldseek_dir> <output_fa> [threads]"
    echo
    echo "Arguments:"
    echo "  cif_dir       Directory containing CIF files"
    echo "  foldseek_dir  Directory for intermediate Foldseek database files"
    echo "  output_fa     Output FASTA file with 3DI sequences"
    echo "  threads       Number of threads (default: 16)"
    exit 1
}

[[ $# -lt 3 ]] && usage

CIF_DIR="$1"
FOLDSEEK_DIR="$2"
OUTPUT_FA="$3"
THREADS="${4:-16}"

mkdir -p "${FOLDSEEK_DIR}"

foldseek createdb "${CIF_DIR}" "${FOLDSEEK_DIR}/3di_raw.db" -v 3 --threads "${THREADS}"
foldseek lndb "${FOLDSEEK_DIR}/3di_raw.db_h" "${FOLDSEEK_DIR}/3di_raw.db_ss_h" -v 3
foldseek convert2fasta "${FOLDSEEK_DIR}/3di_raw.db_ss" "${OUTPUT_FA}" -v 3
foldseek convert2fasta "${FOLDSEEK_DIR}/3di_raw.db" "${FOLDSEEK_DIR}/aa_foldseek.fa" -v 3
