# Metabolic Lung Cancer MR Analysis - R Package Dependencies
# This file contains all R packages required for the analysis

# Core MR Analysis Packages
install.packages("TwoSampleMR")      # Two-sample Mendelian randomization
install.packages("MendelianRandomization")  # Alternative MR package
install.packages("MR-PRESSO")        # MR-PRESSO pleiotropy test
install.packages("causalacc")        # Causal inference methods

# Data Processing and Analysis
install.packages("data.table")       # Fast data manipulation
install.packages("dplyr")           # Data manipulation
install.packages("tidyr")           # Data tidying
install.packages("readr")           # Read rectangular data
install.packages("readxl")          # Read Excel files
install.packages("writexl")         # Write Excel files

# Statistical Analysis
install.packages("ggplot2")         # Data visualization
install.packages("forestplot")      # Forest plots
install.packages("meta")            # Meta-analysis
install.packages("corrplot")        # Correlation plots
install.packages("heatmap3")        # Enhanced heatmaps

# GWAS Specific Packages
install.packages("ieugwasr")        # Interface to IEU GWAS database
install.packages("genetics.binaR")  # Genetic data processing
install.packages("LDlinkR")         # Linkage disequilibrium
install.packages("gwasglue")        # GWAS analysis glue

# Functional Annotation
install.packages("clusterProfiler") # Functional enrichment
install.packages("ReactomePA")      # Reactome pathway analysis
install.packages("org.Hs.eg.db")    # Human gene annotations

# Reporting and Documentation
install.packages("rmarkdown")       # Dynamic documents
install.packages("knitR")           # Report generation
install.packages("stargazer")       # Regression tables
install.packages("kableExtra")      # Enhanced tables

# Development and Testing
install.packages("devtools")        # Development tools
install.packages("testthat")        # Testing framework
install.packages("roxygen2")        # Documentation

# Utility Packages
install.packages("tidyverse")       # Tidyverse collection
install.packages("magrittr")        # Pipe operators
install.packages("purrr")           # Functional programming
install.packages("stringr")         # String manipulation
install.packages("lubridate")       # Date/time handling
install.packages("here")            # File path management