$ErrorActionPreference = "Stop"
conda env create -f environment.yml
conda activate tcga-methylation-targeted
Rscript src/synthetic_data.R
Rscript analysis.R
