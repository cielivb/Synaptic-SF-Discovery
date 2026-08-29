#!/usr/bin/env bash

set -e # Exit on error
source "$(conda info --base)/etc/profile.d/conda.sh" # Locate conda
conda activate Human-NT-Receptor-Geometry-by-Synapse-Type-Env

mkdir -p data_pdb data_csv results

# Download .pdb files for Homo Sapiens proteins
if [ -z "$(ls -A data_pdb)" ]; then
    proteomes2structs --pdb UP000005640 data_pdb
else
    echo ".pdb data already downloaded - skipping download"
fi

# Extract coordinates from each .pdb file
if [ -z "$(ls -A data_csv)" ]; then
    ./scripts/extract_coordinates.sh data_pdb data_csv
else
    echo ".csv data already exists - skipping coordinate extraction"
fi

# Execute R analysis
Rscript scripts/analysis.R

# Done
conda deactivate
