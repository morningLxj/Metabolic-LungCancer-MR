# Targeted Promoter Methylation Validation in TCGA LUAD/LUSC

## Overview
- Probe-level promoter methylation analysis for MFAP2, CDK11A, PRKCZ, WRAP73.
- Validation of negative correlation between promoter Beta and RNA TPM.
- Outputs include Table S7-style CSVs and publication-ready figures.

## Methodology
- HM450K annotation merged with EPIC when available; standardized fields: Probe, Gene, Region.
- Promoter regions only: TSS1500, TSS200, 5'UTR, 1stExon.
- Strict synonym whitelist to avoid off-target aliases:
  - MFAP2: MFAP2, MAGP, MAGP1
  - CDK11A: CDK11A, CDC2L2, CDK11
  - WRAP73: WRAP73, WRAP53
  - PRKCZ: PRKCZ, PKC-ZETA, PKC2
- Differential methylation per probe: Δβ = mean(Tumor) – mean(Normal), Wilcoxon rank-sum, BH FDR.
- Best probe per gene: prioritize P < 0.05 and |Δβ| > 0.1; otherwise min FDR and max |Δβ|.
- Spearman correlation between promoter Beta and TPM; combined and per-cohort analyses.
 - Full manuscript can be placed under `docs/sn-article.pdf`. If you prefer Markdown, convert with Pandoc: `pandoc docs/sn-article.pdf -o docs/sn-article.md`.

## Data Description
- HM450K methylation TSV: first column `Composite Element REF` (probe ID), subsequent columns TCGA barcodes.
- RNA TPM TSV: first column `Ensembl_ID`, subsequent columns TCGA barcodes; Ensembl version suffix removed.
- Annotation via `minfi`: 450k primary, EPIC fallback.

## Reproducibility
- Environment pinned in `environment.yml` (Conda + Bioconductor).
- Configuration in `config.yaml`: paths, promoter tags, thresholds, synonyms.
- Synthetic data generator in `src/synthetic_data.R` enables end-to-end execution without private data.

## How to Run
1. Install Conda.
2. `conda env create -f environment.yml && conda activate tcga-methylation-targeted`
3. Place raw TSVs under `data/raw/` using filenames specified in `config.yaml`.
4. If real data cannot be shared, run:
   - `Rscript src/synthetic_data.R` (writes demo TSVs to `data/raw/`)
5. Run analysis:
   - `Rscript analysis.R`
6. Results:
   - Tables in `results/tables/*.csv`
   - Figures in `results/figures/*.png`

## System Requirements
- R 4.3+
- Conda with channels: conda-forge, bioconda
- Bioconductor packages listed in `environment.yml`

## Citation
- Authors: Your Name et al.
- Title: Targeted Promoter Methylation Validation in TCGA LUAD/LUSC
- Year: 2025
- DOI: TBD

## License
- MIT License, see `LICENSE`.

## Privacy Notice
- Real patient data is private. Synthetic data with identical schema is provided to allow reviewers to run the pipeline.
