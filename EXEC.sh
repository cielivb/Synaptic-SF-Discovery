#!/usr/bin/env bash
# Run scripts 01 through to 10

set -e # Exit on error
source "$(conda info --base)/etc/profile.d/conda.sh" # Locate conda
conda activate Synaptic-SF-Env
R --version # Validate environment
python --version # Validate environment


echo "Environment built successfully"

Rscript scripts/01_filter_proteins.R
Rscript scripts/02_extract_shape_features.R
Rscript scripts/03_extract_go_module_scores.R
Rscript scripts/04_cluster_by_gms.R
Rscript scripts/05_describe_gms_cluster_func.R
Rscript scripts/06_make_taxa_plots.R
Rscript scripts/07_make_shape_plots.R
Rscript scripts/08_core_analysis.R
Rscript scripts/09_make_gms_cluster_summaries.R
Rscript scripts/10_make_misc.R

echo "Analysis complete!"
