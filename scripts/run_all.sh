#!/usr/bin/env bash
set -euo pipefail
conda env create -f environment.yml || true
conda activate tcga-methylation-targeted
Rscript src/synthetic_data.R
Rscript analysis.R
