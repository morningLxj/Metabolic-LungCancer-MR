# Mendelian Randomization Analysis: Metabolic Traits and Lung Cancer Subtypes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)

This repository contains the complete R code for the Mendelian randomization analysis presented in the manuscript:

> **"Causal Effects of Metabolic Traits on Lung Cancer Subtypes: A Comprehensive Mendelian Randomization Study"**

## 📋 Project Overview

This study investigates the causal relationships between 17 metabolic traits and 3 lung cancer subtypes using two-sample Mendelian randomization.

### Key Findings
- **BMI** → Increased risk of squamous cell lung cancer (OR = 1.43, p = 1.2×10⁻⁶)
- **LDL cholesterol** → Increased risk of lung adenocarcinoma (OR = 1.12, p = 0.0026)  
- **Apolipoprotein A1** → Protective effect against squamous cell lung cancer (OR = 0.91, p = 4.8×10⁻⁴)

## 🏗️ Repository Structure

```
Metabolic-LungCancer-MR/
├── 01_data_preparation/     # Data loading and IV selection
├── 02_mr_analysis/          # Primary MR analysis
├── 03_sensitivity_analysis/ # Sensitivity tests
├── 04_visualization/        # Figure generation
├── 05_results/             # Result tables
├── config/                 # Analysis parameters
├── scripts/                # Utility functions
└── logs/                   # Analysis logs
```

## 🚀 Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/[您的用户名]/Metabolic-LungCancer-MR.git
cd Metabolic-LungCancer-MR
```

### 2. Install Dependencies
```r
# Install required R packages
install.packages(c("TwoSampleMR", "MRPRESSO", "ggplot2", "dplyr", 
                   "readr", "openxlsx", "forestplot"))
```

### 3. Run Complete Analysis
```r
# Execute the full analysis pipeline
source("main_analysis.R")
```

### 4. Run Specific Components
```r
# Data preparation only
source("01_data_preparation/01_load_gwas_data.R")

# MR analysis only  
source("02_mr_analysis/01_primary_mr.R")

# Generate figures only
source("04_visualization/01_generate_figures.R")
```

## 📊 Analysis Pipeline

### Step 1: Data Preparation
- Load GWAS summary statistics from OpenGWAS
- Select instrumental variables (p < 5×10⁻⁸, LD r² < 0.001)
- Harmonize exposure and outcome data

### Step 2: MR Analysis  
- Primary analysis using Inverse-Variance Weighted method
- Supplementary methods: MR-Egger, Weighted Median
- Calculate odds ratios and confidence intervals

### Step 3: Sensitivity Analysis
- Heterogeneity tests (Cochran's Q, I²)
- Pleiotropy assessment (MR-Egger intercept)
- Leave-one-out analysis
- MR-PRESSO for outlier detection

### Step 4: Multiple Testing Correction
- Bonferroni correction (p < 9.8×10⁻⁴)
- False Discovery Rate (FDR < 0.05)

## 🔧 Dependencies

- **R** ≥ 4.0.0
- **Key Packages**: 
  - `TwoSampleMR` - MR analysis
  - `MRPRESSO` - Pleiotropy testing  
  - `ggplot2` - Visualization
  - `dplyr` - Data manipulation

See `session_info.txt` for complete package versions.

## 📈 Outputs

The analysis generates:

- **Figures**: Forest plots, heterogeneity scatter plots
- **Tables**: Main results, sensitivity analyses, IV quality metrics
- **Results**: Complete MR estimates with confidence intervals

## 🎯 Data Sources

All GWAS data are sourced from public repositories:

| Trait | GWAS ID | Sample Size |
|-------|---------|-------------|
| Body Mass Index | ieu-b-40 | 681,275 |
| LDL Cholesterol | ieu-b-110 | 440,546 |
| Lung Adenocarcinoma | ieu-a-984 | 65,864 |
| ... | ... | ... |

See `config/analysis_parameters.R` for complete list.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{your2024metabolic,
  title={Causal Effects of Metabolic Traits on Lung Cancer Subtypes: A Comprehensive Mendelian Randomization Study},
  author={Your Name and Coauthors},
  journal={Journal Name},
  year={2024}
}
```

## 🙏 Acknowledgments

- GWAS consortia for making summary statistics publicly available
- Developers of the TwoSampleMR and MR-PRESSO packages
- Study participants and researchers involved in original GWAS