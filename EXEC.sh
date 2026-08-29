#!/usr/bin/env bash
conda activate Human-NT-Receptor-Geometry-by-Synapse-Type-Env

mkdir -p data_pdb data_csv results

# Download .pdb files for Homo Sapiens proteins
proteomes2structs --pdb UP000005640 data_pdb

# Extract coordinates from each .pdb file
extract_coordinates.sh data_pdb data_csv

# Execute R analysis
Rscript scripts/analysis.R

conda deactivate
