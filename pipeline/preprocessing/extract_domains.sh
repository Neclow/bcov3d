#!/bin/bash

set -e

infile=$1
outfile=$2
version=${3:-5.77-108.0}

mkdir -p "$(dirname "$outfile")" temp

sudo docker run --rm \
    -v "$(pwd)/data/interproscan-${version}/data:/opt/interproscan/data" \
    -v "$(pwd):/input" \
    -v "$(pwd)/temp:/temp" \
    "interpro/interproscan:${version}" \
    --input "/input/${infile}" \
    --output-file-base "/input/${outfile%.tsv}" \
    --formats tsv \
    --tempdir /temp \
    --goterms --pathways \
    --cpu 16

echo "Done. Output: $outfile"
