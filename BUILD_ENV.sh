#!/usr/bin/env bash
# Create and activate conda environment

set -e # Exit on error
source "$(conda info --base)/etc/profile.d/conda.sh" # Locate conda

mamba env create -f environment.yml # Note: replace mamba with conda if required
conda activate Human-NT-Receptor-Geometry-by-Synapse-Type-Env
R --version # Validate environment

echo "Environment built successfully"
