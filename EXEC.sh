#!/usr/bin/env bash

set -e # Exit on error
source "$(conda info --base)/etc/profile.d/conda.sh" # Locate conda
conda activate Human-NT-Receptor-Geometry-by-Synapse-Type-Env

mkdir -p data results

# Execute R analysis
Rscript scripts/analysis.R

# Done
conda deactivate
