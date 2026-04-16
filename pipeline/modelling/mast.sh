#!/bin/bash

set -e

s=$1          # input sequence (subset.fa)
m=$2          # model string (from best_model.txt)
te=$3         # input trees file (cut.trees)
prefix=$4     # output prefix for iterative MAST results
seed=$5
threads=$6
metric=${7:-BIC}

# Step 1: Initial MAST run with all trees
echo "=== Initial MAST run ==="
iqtree -s "$s" -m "${m}+T" -te "$te" \
    -wspm -wslm \
    -seed "$seed" -nt "$threads" \
    --prefix "${prefix}_start" --quiet --redo

# Step 2: Sort trees by site log-likelihood
echo "=== Sorting trees ==="
python -m pipeline.modelling.mast_sort \
    -wslm "${prefix}_start.sitelh" \
    -te "$te"

# Step 3: Iterative MAST
echo "=== Iterative MAST ==="
mkdir -p "$prefix"
python -m pipeline.modelling.mast_iterative \
    -s "$s" -m "$m" \
    -te "${te}.sorted" \
    --prefix_dir "$prefix" \
    --seed "$seed" --threads "$threads" --metric "$metric"

echo "Done."
