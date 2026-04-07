# Metabolic-LungCancer-MR

Integrated causal-to-translational multi-omics framework for lung squamous cell carcinoma, spanning Mendelian randomization, mQTL triangulation, multi-omics validation, M12 prognostic modeling, and submission-ready manuscript generation.

[![R](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://cran.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository has been synchronized to the current submission-ready version of the project. The finalized analytical narrative is:

- causal inference across metabolic and inflammatory liabilities
- epigenetic triangulation using mQTL and Bayesian colocalization
- multi-omics interpretation focused on MFAP2-linked biology
- construction and validation of the 12-CpG prognostic signature (M12)
- manuscript, tables, and figure assets aligned with the current submission package

The current paper package is centered on:

- **Model**: M12 / M12_Top12
- **Training cohort**: TCGA-LUSC
- **External validation cohorts**: GSE39279 and GSE30219
- **Key performance metrics**: C-index = 0.664; 3-year AUC = 0.685; calibration slope = 1.227; calibration R² = 0.843

## Repository Layout

```text
.
├── README.md
├── R_packages.R
├── requirements.txt
├── environment.yml
├── docs/
│   ├── manuscript/                # synchronized submission-ready manuscript assets
│   └── methods/
├── results/
│   └── final_submission/          # synchronized tables and selected submission outputs
├── scripts/
│   ├── submission/                # current manuscript/table synchronization scripts
│   ├── R/                         # analytical R scripts retained from the project workflow
│   ├── python/                    # analytical Python scripts retained from the project workflow
│   └── utils/
└── tests/
```

## Synchronized Final Assets

After synchronization, the most authoritative outputs in this repository are:

- `docs/manuscript/Manuscript_Final_Revised.docx`
- `docs/manuscript/Manuscript_Final_Revised.md`
- `docs/manuscript/Methods_Revised.docx`
- `docs/manuscript/Methods_Revised.md`
- `docs/manuscript/Figure_Names_and_Legends_Submission.docx`
- `results/final_submission/tables/Main_Tables_Final.docx`
- `results/final_submission/tables/Supplementary_Tables_Final.docx`
- `results/final_submission/tables/Supplementary_Table_S1.csv` to `Supplementary_Table_S9.csv`
- `results/final_submission/figures/Figure7_M12_Summary.csv`

## Core Submission Scripts

The current submission package is driven by:

- `scripts/submission/generate_full_manuscript_final.py`
- `scripts/submission/generate_methods_docx.py`
- `scripts/submission/rebuild_supplementary_tables_m12.py`
- `scripts/submission/build_and_format_submission_tables_docx.py`
- `scripts/submission/verify_submission_docx_accuracy.py`
- `scripts/submission/rebuild_main_figure5_m12.py`
- `scripts/submission/rebuild_main_figure6_m12.py`
- `scripts/submission/rebuild_main_figure7_m12.R`
- `scripts/submission/regen_s6_enhanced.py`

## Reproducibility Notes

- Public summary-statistics and public molecular resources are required for full reruns.
- The repository includes the submission-facing scripts and synchronized outputs needed to inspect the final analytical state.
- Large raw datasets are intentionally excluded from version control.

## Citation

If you use this repository, cite the corresponding manuscript and repository URL:

- Repository: https://github.com/morningLxj/Metabolic-LungCancer-MR

## License

This project is released under the MIT License. See [LICENSE](LICENSE).

## Maintainer

- GitHub: https://github.com/morningLxj/Metabolic-LungCancer-MR

---

Last synchronized to the submission-ready M12 version: 2026-04-08
