#!/bin/bash

set -e

s=$1
prefix=$2
mset=$3
mfreq=$4
mrate=$5
seed=$6
threads=$7

iqtree -s "$s" -m MFP \
    -mset "$mset" -mfreq "$mfreq" -mrate "$mrate" \
    -seed "$seed" -nt "$threads" \
    --prefix "$prefix" --quiet --redo

best_model=$(grep "Best-fit model according to BIC:" "$prefix.iqtree" | sed 's/.*: //')
echo "$best_model" > "$(dirname "$prefix")/best_model.txt"

echo "Best model: $best_model"
