#!/usr/bin/env Rscript
# =============================================================================
# Main Analysis Script for Metabolic Lung Cancer MR Analysis
# =============================================================================
# This script orchestrates the complete workflow for MR analysis
# Author: morningLxj
# Date: 2024-12-01
# =============================================================================

#' Main Analysis Pipeline
#'
#' This function runs the complete MR analysis pipeline including:
#' 1. Data preparation and QC
#' 2. Instrument extraction
#' 3. Two-sample MR analysis
#' 4. Sensitivity analysis
#' 5. Mediation analysis
#' 6. Colocalization analysis
#' 7. Results summarization
#'
#' @param exposure_gwas Path to exposure GWAS summary statistics
#' @param outcome_gwas Path to outcome GWAS summary statistics
#' @param output_dir Output directory for results
#' @param verbose Whether to print detailed progress
#' @return List containing all analysis results
#' @export
run_complete_analysis <- function(exposure_gwas = "data/raw/exposure_gwas.txt",
                                 outcome_gwas = "data/raw/outcome_gwas.txt",
                                 output_dir = "results/",
                                 verbose = TRUE) {
  
  # Initialize logging
  if (verbose) {
    cat("Starting Metabolic Lung Cancer MR Analysis Pipeline\n")
    cat("=================================================\n")
    cat("Timestamp:", as.character(Sys.time()), "\n\n")
  }
  
  # Load required packages
  check_and_load_packages()
  
  # Set output directories
  results_dir <- file.path(output_dir, "mr_results")
  figures_dir <- file.path(output_dir, "figures")
  tables_dir <- file.path(output_dir, "tables")
  reports_dir <- file.path(output_dir, "reports")
  
  # Create directories if they don't exist
  dirs_to_create <- c(results_dir, figures_dir, tables_dir, reports_dir)
  sapply(dirs_to_create, function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE))
  
  # Initialize results list
  all_results <- list()
  
  tryCatch({
    
    # Step 1: Data Preparation and Quality Control
    if (verbose) cat("Step 1: Data Preparation and QC...\n")
    data_prep_results <- step1_data_preparation(exposure_gwas, outcome_gwas, results_dir)
    all_results$data_prep <- data_prep_results
    
    # Step 2: Instrument Extraction
    if (verbose) cat("Step 2: Instrument Extraction...\n")
    instruments <- step2_extract_instruments(data_prep_results$exposure_data, results_dir)
    all_results$instruments <- instruments
    
    # Step 3: Two-Sample MR Analysis
    if (verbose) cat("Step 3: Two-Sample MR Analysis...\n")
    mr_results <- step3_run_mr_analysis(instruments, data_prep_results$outcome_data, results_dir)
    all_results$mr_results <- mr_results
    
    # Step 4: Sensitivity Analysis
    if (verbose) cat("Step 4: Sensitivity Analysis...\n")
    sensitivity_results <- step4_sensitivity_analysis(mr_results, results_dir)
    all_results$sensitivity <- sensitivity_results
    
    # Step 5: Mediation Analysis
    if (verbose) cat("Step 5: Mediation Analysis...\n")
    mediation_results <- step5_mediation_analysis(mr_results, results_dir)
    all_results$mediation <- mediation_results
    
    # Step 6: Colocalization Analysis
    if (verbose) cat("Step 6: Colocalization Analysis...\n")
    coloc_results <- step6_colocalization_analysis(instruments, results_dir)
    all_results$coloc <- coloc_results
    
    # Step 7: Results Summarization
    if (verbose) cat("Step 7: Results Summarization...\n")
    summary_report <- step7_summarize_results(all_results, output_dir)
    all_results$summary <- summary_report
    
    # Generate final report
    if (verbose) cat("Generating Final Report...\n")
    generate_final_report(all_results, reports_dir)
    
    if (verbose) {
      cat("\n=================================================\n")
      cat("Analysis completed successfully!\n")
      cat("Results saved to:", output_dir, "\n")
      cat("Timestamp:", as.character(Sys.time()), "\n")
    }
    
    return(all_results)
    
  }, error = function(e) {
    cat("ERROR in analysis pipeline:", conditionMessage(e), "\n")
    save(all_results, file = file.path(output_dir, "error_backup.RData"))
    stop(e)
  })
}

#' Check and Load Required Packages
#'
#' @keywords internal
check_and_load_packages <- function() {
  required_packages <- c(
    "TwoSampleMR", "MendelianRandomization", "data.table", 
    "dplyr", "ggplot2", "forestplot", "readr", "writexl"
  )
  
  missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
  
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages)
  }
  
  # Load all packages
  sapply(required_packages, library, character.only = TRUE)
}

#' Step 1: Data Preparation and QC
#'
#' @param exposure_gwas Path to exposure GWAS file
#' @param outcome_gwas Path to outcome GWAS file
#' @param output_dir Output directory
#' @keywords internal
step1_data_preparation <- function(exposure_gwas, outcome_gwas, output_dir) {
  cat("Loading exposure data from:", exposure_gwas, "\n")
  exposure_data <- read_exposure_data(exposure_gwas)
  
  cat("Loading outcome data from:", outcome_gwas, "\n")
  outcome_data <- read_outcome_data(outcome_gwas)
  
  # Quality control
  exposure_data <- exposure_data[!is.na(exposure_data$beta), ]
  exposure_data <- exposure_data[exposure_data$pval.exposure < 5e-8, ]
  
  outcome_data <- outcome_data[!is.na(outcome_data$beta), ]
  
  # Save processed data
  save(exposure_data, outcome_data, 
       file = file.path(output_dir, "step1_processed_data.RData"))
  
  return(list(exposure_data = exposure_data, outcome_data = outcome_data))
}

#' Step 2: Instrument Extraction
#'
#' @param exposure_data Processed exposure data
#' @param output_dir Output directory
#' @keywords internal
step2_extract_instruments <- function(exposure_data, output_dir) {
  instruments <- extract_instruments(exposure_data)
  cat("Extracted", nrow(instruments), "instruments\n")
  
  save(instruments, file = file.path(output_dir, "step2_instruments.RData"))
  return(instruments)
}

#' Step 3: Two-Sample MR Analysis
#'
#' @param instruments Instrumental variables
#' @param outcome_data Processed outcome data
#' @param output_dir Output directory
#' @keywords internal
step3_run_mr_analysis <- function(instruments, outcome_data, output_dir) {
  dat <- harmonise_data(instruments, outcome_data)
  mr_results <- mr(dat)
  
  # Additional methods
  mr_results_ivw <- mr(dat, method_list = "mr_ivw")
  mr_results_egger <- mr(dat, method_list = "mr_egger_regression")
  mr_results_weighted_median <- mr(dat, method_list = "mr_weighted_median")
  
  # Combine results
  all_mr_results <- mr_results %>%
    bind_rows(mr_results_ivw) %>%
    bind_rows(mr_results_egger) %>%
    bind_rows(mr_results_weighted_median)
  
  # Save results
  write.table(all_mr_results, 
              file = file.path(output_dir, "step3_mr_results.csv"),
              sep = ",", row.names = FALSE)
  
  save(all_mr_results, file = file.path(output_dir, "step3_mr_results.RData"))
  
  return(all_mr_results)
}

#' Step 4: Sensitivity Analysis
#'
#' @param mr_results MR analysis results
#' @param output_dir Output directory
#' @keywords internal
step4_sensitivity_analysis <- function(mr_results, output_dir) {
  # MR-PRESSO test
  presso_results <- run_mr_presso(mr_results)
  
  # Horizontal pleiotropy test
  pleiotropy_test <- mr_pleiotropy_test(mr_results)
  
  # Heterogeneity test
  het_test <- mr_heterogeneity(mr_results)
  
  # Save sensitivity results
  sensitivity_results <- list(
    presso = presso_results,
    pleiotropy = pleiotropy_test,
    heterogeneity = het_test
  )
  
  save(sensitivity_results, file = file.path(output_dir, "step4_sensitivity.RData"))
  
  return(sensitivity_results)
}

#' Step 5: Mediation Analysis
#'
#' @param mr_results MR analysis results
#' @param output_dir Output directory
#' @keywords internal
step5_mediation_analysis <- function(mr_results, output_dir) {
  # This would be implemented based on specific mediation analysis requirements
  # For now, return placeholder
  mediation_results <- list(message = "Mediation analysis results would be generated here")
  
  save(mediation_results, file = file.path(output_dir, "step5_mediation.RData"))
  return(mediation_results)
}

#' Step 6: Colocalization Analysis
#'
#' @param instruments Instrumental variables
#' @param output_dir Output directory
#' @keywords internal
step6_colocalization_analysis <- function(instruments, output_dir) {
  # This would be implemented using coloc package or similar
  # For now, return placeholder
  coloc_results <- list(message = "Colocalization analysis results would be generated here")
  
  save(coloc_results, file = file.path(output_dir, "step6_coloc.RData"))
  return(coloc_results)
}

#' Step 7: Results Summarization
#'
#' @param all_results All analysis results
#' @param output_dir Output directory
#' @keywords internal
step7_summarize_results <- function(all_results, output_dir) {
  # Generate summary statistics
  summary_stats <- list(
    n_instruments = nrow(all_results$instruments),
    n_mr_methods = nrow(all_results$mr_results),
    significant_associations = sum(all_results$mr_results$pval < 0.05)
  )
  
  # Save summary
  save(summary_stats, file = file.path(output_dir, "step7_summary.RData"))
  return(summary_stats)
}

#' Generate Final Report
#'
#' @param all_results All analysis results
#' @param output_dir Output directory
#' @keywords internal
generate_final_report <- function(all_results, output_dir) {
  report_content <- paste0(
    "# Metabolic Lung Cancer MR Analysis Report\n\n",
    "## Analysis Summary\n\n",
    "- Total instruments: ", nrow(all_results$instruments), "\n",
    "- MR methods tested: ", nrow(all_results$mr_results), "\n",
    "- Significant associations: ", sum(all_results$mr_results$pval < 0.05), "\n\n",
    "## Main Results\n\n",
    "Detailed results are available in the individual output files.\n\n",
    "---\n",
    "Generated on: ", as.character(Sys.time()), "\n"
  )
  
  writeLines(report_content, file.path(output_dir, "final_report.md"))
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript main_analysis.R [exposure_gwas] [outcome_gwas] [output_dir]\n")
    cat("Running with default parameters...\n\n")
    
    # Run with default parameters
    results <- run_complete_analysis()
    
  } else if (length(args) == 3) {
    results <- run_complete_analysis(args[1], args[2], args[3])
    
  } else {
    stop("Invalid arguments. Usage: Rscript main_analysis.R [exposure_gwas] [outcome_gwas] [output_dir]")
  }
}