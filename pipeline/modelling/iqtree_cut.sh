#!/bin/bash

set -e

facutdir=$1
treecutdir=$2
m=$3
seed=$4
threads=$5

# Basic IQ-TREE runs
files="${facutdir}/*.fa"
for s in $files; do
    echo "Running file ${s}"
    iqtree -s "$s" -m "$m" -seed "$seed" \
        -nt "$threads" \
        --prefix "$treecutdir/$(basename "$s")" \
        --quiet --redo
done

cat "$treecutdir"/*.treefile > "$(dirname "$treecutdir")/$(basename "$treecutdir").trees"

echo "Done."
