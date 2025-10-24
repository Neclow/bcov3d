#!/bin/bash

# Install extra R packages if missing
Rscript -e 'if (!requireNamespace("envalysis", quietly = TRUE)) \
install.packages("envalysis", repos="https://cloud.r-project.org")'

Rscript -e 'if (!requireNamespace("brms", quietly = TRUE)) \
install.packages("brms")'

# Initialize and update submodule recursively
git submodule update --init --recursive extern/famsa3di

# Build famsa3di
cd extern/famsa3di
make
