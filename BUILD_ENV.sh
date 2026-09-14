#!/usr/bin/env bash
# Create and activate conda environment

set -e # Exit on error
source "$(conda info --base)/etc/profile.d/conda.sh" # Locate conda

mamba env create -f environment.yml # Note: replace mamba with conda if required
conda activate Synaptic-SF-Env
R --version # Validate environment
python --version

echo "Environment built successfully"
