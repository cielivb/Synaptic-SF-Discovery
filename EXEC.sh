#!usr/bin/bash
Rscript -e "rmarkdown::render('scripts/analysis.Rmd', output_dir='results')"
