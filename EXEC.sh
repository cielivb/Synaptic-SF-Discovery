#!/usr/bin/env bash
conda activate Human-NT-Receptor-Geometry-by-Synapse-Type-Env

# Download .pdb files for Homo Sapiens proteins
mkdir -p data
proteomes2structs --pdb UP000005640 data

# Execute R analysis
Rscript scripts/analysis.R

conda deactivate
