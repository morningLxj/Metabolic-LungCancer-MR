############################################################################
# Step 7b: Mediation Analysis - Sensitivity Analysis (English Version)
# Mendelian Randomization Study - Metabolic Traits, Inflammatory Markers, and Lung Cancer Subtypes
# 
# Contents:
#   1. Load mediation analysis results from Step 7
#   2. Sensitivity analyses for significant mediation pathways
#   3. Generate English version figures for publication
#   4. Heterogeneity and pleiotropy tests
#   5. Leave-one-out analysis
#   6. MR method comparison
############################################################################

cat("Step 7b: Mediation Analysis - Sensitivity Analysis (English Version)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Note: This script performs comprehensive sensitivity analyses on mediation results\n\n")

# Load necessary packages
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

# Declare global variables
utils::globalVariables(c("exposure", "outcome", "method", "pathway_label",
                        "SNP", "b", "se", "leave_one_out_effect"))

# Create output directories
dirs <- c(
  "results/mediation_analysis/sensitivity",
  "results/figures/step07b_sensitivity",
  "results/figures/step07b_english",
  "results/tables/step07b_sensitivity"
)
for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Set timeout and options
options(timeout = 300)

# ============================================================================
# Step 1: Load Step 7 mediation results
# ============================================================================

cat("【Step 1】Load mediation analysis results from Step 7\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Check if Step 7 results exist
if (!file.exists("data/step07_all_mediation_results.RData")) {
  stop("Error: Step 7 mediation results not found. Please run Step 7 first.")
}

# Load results
load("data/step07_all_mediation_results.RData")
cat(sprintf("✓ Loaded mediation results: %d pathways\n", nrow(mediation_summary)))
cat(sprintf("✓ Loaded full results: %d pathway objects\n", length(mediation_results)))

# Check if mediation_summary has required columns
required_cols <- c("exposure", "mediator", "outcome", "indirect_effect", 
                   "indirect_effect_pval", "significant_fdr")
missing_cols <- setdiff(required_cols, names(mediation_summary))
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Missing required columns: %s", 
              paste(missing_cols, collapse = ", ")))
}

cat("\n")

# ============================================================================
# Step 2: Define English label mapping
# ============================================================================

cat("【Step 2】Define English label mapping\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Exposure name mapping (Chinese to English)
exposure_labels <- c(
  circulating_leptin = "Circulating Leptin",
  vitamin_D = "Vitamin D",
  HbA1c = "HbA1c",
  ApoB = "Apolipoprotein B",
  ApoA1 = "Apolipoprotein A-I",
  IGF1 = "IGF-1",
  ApoB_ApoA1_ratio = "ApoB/ApoA-I Ratio",
  HDL_diameter = "HDL Diameter",
  HDL_large = "Large HDL",
  remnant_cholesterol = "Remnant Cholesterol",
  LDL_small = "Small LDL",
  BCAA = "Branched-Chain Amino Acids",
  HDL_very_large = "Very Large HDL",
  BMI = "Body Mass Index",
  HDL_cholesterol = "HDL Cholesterol",
  LDL_cholesterol = "LDL Cholesterol",
  smoking_initiation = "Smoking Initiation",
  alcohol_drinks = "Alcohol Consumption",
  fasting_glucose = "Fasting Glucose",
  fasting_insulin = "Fasting Insulin",
  SBP = "Systolic Blood Pressure",
  DBP = "Diastolic Blood Pressure",
  hypertension = "Hypertension",
  triglycerides = "Triglycerides",
  GGT = "Gamma-Glutamyl Transferase"
)

# Mediator name mapping
mediator_labels <- c(
  CRP = "C-Reactive Protein",
  WBC = "White Blood Cell Count",
  IL6 = "Interleukin-6",
  IL6R = "IL-6 Receptor",
  TNFR1 = "TNF Receptor 1"
)

# Outcome name mapping
outcome_labels <- c(
  lung_cancer_overall = "Overall Lung Cancer",
  lung_adenocarcinoma = "Lung Adenocarcinoma",
  squamous_cell_lung = "Squamous Cell Lung Cancer"
)

# Strategy name mapping
strategy_labels <- c(
  smart_selection = "Evidence-Based Selection",
  exhaustive_screening = "Exhaustive Screening"
)

# Add English labels to mediation_summary
mediation_summary <- mediation_summary %>%
  mutate(
    exposure_label = ifelse(exposure %in% names(exposure_labels),
                           exposure_labels[exposure],
                           exposure),
    mediator_label = ifelse(mediator %in% names(mediator_labels),
                           mediator_labels[mediator],
                           mediator),
    outcome_label = ifelse(outcome %in% names(outcome_labels),
                          outcome_labels[outcome],
                          outcome),
    strategy_label = ifelse(source %in% names(strategy_labels),
                           strategy_labels[source],
                           source)
  )

cat("✓ English labels added to mediation summary\n\n")

# ============================================================================
# Step 3: Heterogeneity and pleiotropy tests
# ============================================================================

cat("【Step 3】Heterogeneity and pleiotropy tests\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Filter significant pathways for sensitivity analysis
sig_pathways <- mediation_summary %>%
  filter(significant_fdr == TRUE | all_paths_significant == TRUE) %>%
  arrange(fdr_pval_indirect)

cat(sprintf("Selected %d significant pathways for sensitivity analysis\n", 
           nrow(sig_pathways)))

if (nrow(sig_pathways) > 0) {
  
  # Initialize results data frame
  sensitivity_results <- data.frame()
  
  cat("\nPerforming heterogeneity and pleiotropy tests...\n")
  pb <- txtProgressBar(min = 0, max = nrow(sig_pathways), style = 3)
  
  for (i in seq_len(nrow(sig_pathways))) {
    pathway <- sig_pathways[i, ]
    
    setTxtProgressBar(pb, i)
    
    # Get the corresponding mediation result
    key <- paste(pathway$exposure, pathway$mediator, pathway$outcome, sep = "_")
    
    if (key %in% names(mediation_results)) {
      result <- mediation_results[[key]]
      
      # Extract heterogeneity and pleiotropy for each step
      # Step 1: Exposure -> Mediator
      if (!is.null(result$exp_to_med) && 
          !is.null(result$exp_to_med$heterogeneity)) {
        het_exp_med <- result$exp_to_med$heterogeneity %>%
          filter(method == "Inverse variance weighted")
        
        Q_exp_med <- if (nrow(het_exp_med) > 0) het_exp_med$Q[1] else NA
        Q_pval_exp_med <- if (nrow(het_exp_med) > 0) het_exp_med$Q_pval[1] else NA
      } else {
        Q_exp_med <- NA
        Q_pval_exp_med <- NA
      }
      
      if (!is.null(result$exp_to_med) && 
          !is.null(result$exp_to_med$pleiotropy)) {
        pleio_exp_med <- result$exp_to_med$pleiotropy
        egger_int_exp_med <- if (nrow(pleio_exp_med) > 0) {
          pleio_exp_med$egger_intercept[1]
        } else NA
        egger_pval_exp_med <- if (nrow(pleio_exp_med) > 0) {
          pleio_exp_med$pval[1]
        } else NA
      } else {
        egger_int_exp_med <- NA
        egger_pval_exp_med <- NA
      }
      
      # Step 2: Mediator -> Outcome
      if (!is.null(result$med_to_out) && 
          !is.null(result$med_to_out$heterogeneity)) {
        het_med_out <- result$med_to_out$heterogeneity %>%
          filter(method == "Inverse variance weighted")
        
        Q_med_out <- if (nrow(het_med_out) > 0) het_med_out$Q[1] else NA
        Q_pval_med_out <- if (nrow(het_med_out) > 0) het_med_out$Q_pval[1] else NA
      } else {
        Q_med_out <- NA
        Q_pval_med_out <- NA
      }
      
      if (!is.null(result$med_to_out) && 
          !is.null(result$med_to_out$pleiotropy)) {
        pleio_med_out <- result$med_to_out$pleiotropy
        egger_int_med_out <- if (nrow(pleio_med_out) > 0) {
          pleio_med_out$egger_intercept[1]
        } else NA
        egger_pval_med_out <- if (nrow(pleio_med_out) > 0) {
          pleio_med_out$pval[1]
        } else NA
      } else {
        egger_int_med_out <- NA
        egger_pval_med_out <- NA
      }
      
      # Step 3: Exposure -> Outcome
      if (!is.null(result$exp_to_out) && 
          !is.null(result$exp_to_out$heterogeneity)) {
        het_exp_out <- result$exp_to_out$heterogeneity %>%
          filter(method == "Inverse variance weighted")
        
        Q_exp_out <- if (nrow(het_exp_out) > 0) het_exp_out$Q[1] else NA
        Q_pval_exp_out <- if (nrow(het_exp_out) > 0) het_exp_out$Q_pval[1] else NA
      } else {
        Q_exp_out <- NA
        Q_pval_exp_out <- NA
      }
      
      if (!is.null(result$exp_to_out) && 
          !is.null(result$exp_to_out$pleiotropy)) {
        pleio_exp_out <- result$exp_to_out$pleiotropy
        egger_int_exp_out <- if (nrow(pleio_exp_out) > 0) {
          pleio_exp_out$egger_intercept[1]
        } else NA
        egger_pval_exp_out <- if (nrow(pleio_exp_out) > 0) {
          pleio_exp_out$pval[1]
        } else NA
      } else {
        egger_int_exp_out <- NA
        egger_pval_exp_out <- NA
      }
      
      # Create result row
      sens_row <- data.frame(
        pathway_id = i,
        exposure = pathway$exposure,
        exposure_label = pathway$exposure_label,
        mediator = pathway$mediator,
        mediator_label = pathway$mediator_label,
        outcome = pathway$outcome,
        outcome_label = pathway$outcome_label,
        
        # Exposure -> Mediator
        Q_exp_to_med = Q_exp_med,
        Q_pval_exp_to_med = Q_pval_exp_med,
        egger_intercept_exp_to_med = egger_int_exp_med,
        egger_pval_exp_to_med = egger_pval_exp_med,
        heterogeneity_exp_to_med = ifelse(is.na(Q_pval_exp_med), NA,
                                         ifelse(Q_pval_exp_med < 0.05, 
                                               "Significant", "Not Significant")),
        pleiotropy_exp_to_med = ifelse(is.na(egger_pval_exp_med), NA,
                                      ifelse(egger_pval_exp_med < 0.05, 
                                            "Detected", "Not Detected")),
        
        # Mediator -> Outcome
        Q_med_to_out = Q_med_out,
        Q_pval_med_to_out = Q_pval_med_out,
        egger_intercept_med_to_out = egger_int_med_out,
        egger_pval_med_to_out = egger_pval_med_out,
        heterogeneity_med_to_out = ifelse(is.na(Q_pval_med_out), NA,
                                         ifelse(Q_pval_med_out < 0.05, 
                                               "Significant", "Not Significant")),
        pleiotropy_med_to_out = ifelse(is.na(egger_pval_med_out), NA,
                                      ifelse(egger_pval_med_out < 0.05, 
                                            "Detected", "Not Detected")),
        
        # Exposure -> Outcome
        Q_exp_to_out = Q_exp_out,
        Q_pval_exp_to_out = Q_pval_exp_out,
        egger_intercept_exp_to_out = egger_int_exp_out,
        egger_pval_exp_to_out = egger_pval_exp_out,
        heterogeneity_exp_to_out = ifelse(is.na(Q_pval_exp_out), NA,
                                         ifelse(Q_pval_exp_out < 0.05, 
                                               "Significant", "Not Significant")),
        pleiotropy_exp_to_out = ifelse(is.na(egger_pval_exp_out), NA,
                                      ifelse(egger_pval_exp_out < 0.05, 
                                            "Detected", "Not Detected")),
        
        stringsAsFactors = FALSE
      )
      
      sensitivity_results <- rbind(sensitivity_results, sens_row)
    }
  }
  
  close(pb)
  
  # Save sensitivity results
  write.xlsx(
    sensitivity_results,
    "results/tables/step07b_sensitivity/heterogeneity_pleiotropy_tests.xlsx",
    rowNames = FALSE
  )
  
  write.csv(
    sensitivity_results,
    "results/tables/step07b_sensitivity/heterogeneity_pleiotropy_tests.csv",
    row.names = FALSE
  )
  
  cat("\n✓ Heterogeneity and pleiotropy tests completed\n")
  cat(sprintf("✓ Saved: results/tables/step07b_sensitivity/heterogeneity_pleiotropy_tests.xlsx\n\n"))
  
} else {
  cat("⚠ No significant pathways found for sensitivity analysis\n\n")
  sensitivity_results <- data.frame()
}

# ============================================================================
# Step 4: Leave-one-out sensitivity analysis
# ============================================================================

cat("【Step 4】Leave-one-out sensitivity analysis\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Select top pathways for leave-one-out analysis (computationally intensive)
if (nrow(sig_pathways) > 0) {
  
  # Limit to top 10 most significant pathways
  top_pathways <- head(sig_pathways, 10)
  
  cat(sprintf("Performing leave-one-out analysis for top %d pathways...\n", 
             nrow(top_pathways)))
  
  loo_results_list <- list()
  
  pb <- txtProgressBar(min = 0, max = nrow(top_pathways), style = 3)
  
  for (i in seq_len(nrow(top_pathways))) {
    pathway <- top_pathways[i, ]
    
    setTxtProgressBar(pb, i)
    
    key <- paste(pathway$exposure, pathway$mediator, pathway$outcome, sep = "_")
    
    if (key %in% names(mediation_results)) {
      result <- mediation_results[[key]]
      
      # Leave-one-out for the mediation pathway
      # Focus on Exposure -> Mediator (most critical step)
      if (!is.null(result$exp_to_med) && 
          !is.null(result$exp_to_med$harmonized) &&
          nrow(result$exp_to_med$harmonized) > 3) {
        
        harmonized_data <- result$exp_to_med$harmonized
        n_snps <- nrow(harmonized_data)
        
        loo_effects <- data.frame()
        
        for (j in seq_len(n_snps)) {
          # Remove one SNP
          harmonized_loo <- harmonized_data[-j, ]
          
          # Perform MR
          mr_loo <- tryCatch({
            mr(harmonized_loo, method_list = c("mr_ivw"))
          }, error = function(e) {
            NULL
          })
          
          if (!is.null(mr_loo) && nrow(mr_loo) > 0) {
            ivw_result <- mr_loo %>% filter(method == "Inverse variance weighted")
            
            if (nrow(ivw_result) > 0) {
              loo_effects <- rbind(loo_effects, data.frame(
                removed_snp = harmonized_data$SNP[j],
                loo_beta = ivw_result$b[1],
                loo_se = ivw_result$se[1],
                loo_pval = ivw_result$pval[1],
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if (nrow(loo_effects) > 0) {
          loo_effects$pathway_id <- i
          loo_effects$exposure <- pathway$exposure
          loo_effects$exposure_label <- pathway$exposure_label
          loo_effects$mediator <- pathway$mediator
          loo_effects$mediator_label <- pathway$mediator_label
          loo_effects$outcome <- pathway$outcome
          loo_effects$outcome_label <- pathway$outcome_label
          loo_effects$original_beta <- pathway$exp_to_med_beta
          loo_effects$original_se <- pathway$exp_to_med_se
          
          loo_results_list[[key]] <- loo_effects
        }
      }
      
      Sys.sleep(1)
    }
  }
  
  close(pb)
  
  # Combine all LOO results
  if (length(loo_results_list) > 0) {
    loo_combined <- do.call(rbind, loo_results_list)
    
    # Save LOO results
    write.xlsx(
      loo_combined,
      "results/tables/step07b_sensitivity/leave_one_out_results.xlsx",
      rowNames = FALSE
    )
    
    write.csv(
      loo_combined,
      "results/tables/step07b_sensitivity/leave_one_out_results.csv",
      row.names = FALSE
    )
    
    cat("\n✓ Leave-one-out analysis completed\n")
    cat(sprintf("✓ Analyzed %d pathways\n", length(unique(loo_combined$pathway_id))))
    cat(sprintf("✓ Saved: results/tables/step07b_sensitivity/leave_one_out_results.xlsx\n\n"))
    
    # Create leave-one-out plot for top pathway
    if (nrow(top_pathways) > 0) {
      top_pathway_key <- paste(top_pathways$exposure[1], 
                              top_pathways$mediator[1], 
                              top_pathways$outcome[1], sep = "_")
      
      if (top_pathway_key %in% names(loo_results_list)) {
        loo_data_plot <- loo_results_list[[top_pathway_key]]
        
        pathway_title <- sprintf("%s → %s → %s",
                                top_pathways$exposure_label[1],
                                top_pathways$mediator_label[1],
                                top_pathways$outcome_label[1])
        
        p_loo <- ggplot(loo_data_plot, 
                       aes(x = loo_beta, 
                           y = reorder(removed_snp, loo_beta))) +
          geom_point(size = 2, color = "#E69F00") +
          geom_errorbarh(aes(xmin = loo_beta - 1.96 * loo_se,
                            xmax = loo_beta + 1.96 * loo_se),
                        height = 0.3, alpha = 0.6) +
          geom_vline(xintercept = loo_data_plot$original_beta[1], 
                    linetype = "dashed", color = "red", size = 0.8) +
          labs(
            title = "Leave-One-Out Sensitivity Analysis",
            subtitle = pathway_title,
            x = "Effect Estimate (β) with 95% CI",
            y = "Removed SNP",
            caption = "Red dashed line: original estimate including all SNPs"
          ) +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 7),
            plot.caption = element_text(size = 8, hjust = 0)
          )
        
        ggsave(
          "results/figures/step07b_sensitivity/leave_one_out_example.png",
          p_loo,
          width = 10,
          height = min(45, max(6, nrow(loo_data_plot) * 0.25)),
          dpi = 300,
          limitsize = FALSE
        )
        
        cat("✓ Leave-one-out plot saved\n\n")
      }
    }
    
  } else {
    cat("\n⚠ No leave-one-out results generated\n\n")
    loo_combined <- data.frame()
  }
  
} else {
  cat("⚠ No significant pathways for leave-one-out analysis\n\n")
  loo_combined <- data.frame()
}

# ============================================================================
# Step 5: MR method comparison
# ============================================================================

cat("【Step 5】MR method comparison\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (nrow(sig_pathways) > 0) {
  
  # Extract MR results for different methods
  method_comparison <- data.frame()
  
  cat("Extracting MR estimates from different methods...\n")
  
  for (i in seq_len(nrow(sig_pathways))) {
    pathway <- sig_pathways[i, ]
    
    key <- paste(pathway$exposure, pathway$mediator, pathway$outcome, sep = "_")
    
    if (key %in% names(mediation_results)) {
      result <- mediation_results[[key]]
      
      # Extract all methods for Exposure -> Mediator
      if (!is.null(result$exp_to_med) && 
          !is.null(result$exp_to_med$mr_results)) {
        
        mr_exp_med <- result$exp_to_med$mr_results %>%
          mutate(
            pathway_id = i,
            exposure = pathway$exposure,
            exposure_label = pathway$exposure_label,
            mediator = pathway$mediator,
            mediator_label = pathway$mediator_label,
            outcome = pathway$outcome,
            outcome_label = pathway$outcome_label,
            step = "Exposure → Mediator"
          ) %>%
          select(pathway_id, exposure, exposure_label, mediator, mediator_label,
                outcome, outcome_label, step, method, b, se, pval, nsnp)
        
        method_comparison <- rbind(method_comparison, mr_exp_med)
      }
      
      # Extract all methods for Mediator -> Outcome
      if (!is.null(result$med_to_out) && 
          !is.null(result$med_to_out$mr_results)) {
        
        mr_med_out <- result$med_to_out$mr_results %>%
          mutate(
            pathway_id = i,
            exposure = pathway$exposure,
            exposure_label = pathway$exposure_label,
            mediator = pathway$mediator,
            mediator_label = pathway$mediator_label,
            outcome = pathway$outcome,
            outcome_label = pathway$outcome_label,
            step = "Mediator → Outcome"
          ) %>%
          select(pathway_id, exposure, exposure_label, mediator, mediator_label,
                outcome, outcome_label, step, method, b, se, pval, nsnp)
        
        method_comparison <- rbind(method_comparison, mr_med_out)
      }
    }
  }
  
  if (nrow(method_comparison) > 0) {
    
    # Clean method names
    method_comparison <- method_comparison %>%
      mutate(
        method_clean = case_when(
          method == "Inverse variance weighted" ~ "IVW",
          method == "MR Egger" ~ "MR-Egger",
          method == "Weighted median" ~ "Weighted Median",
          method == "Weighted mode" ~ "Weighted Mode",
          method == "Simple mode" ~ "Simple Mode",
          TRUE ~ method
        )
      )
    
    # Save method comparison
    write.xlsx(
      method_comparison,
      "results/tables/step07b_sensitivity/mr_method_comparison.xlsx",
      rowNames = FALSE
    )
    
    write.csv(
      method_comparison,
      "results/tables/step07b_sensitivity/mr_method_comparison.csv",
      row.names = FALSE
    )
    
    cat(sprintf("✓ MR method comparison completed: %d estimates\n", 
               nrow(method_comparison)))
    cat("✓ Saved: results/tables/step07b_sensitivity/mr_method_comparison.xlsx\n\n")
    
  } else {
    cat("⚠ No MR method comparison data available\n\n")
  }
  
} else {
  cat("⚠ No significant pathways for method comparison\n\n")
  method_comparison <- data.frame()
}

# ============================================================================
# Step 6: Generate English version figures
# ============================================================================

cat("【Step 6】Generate English version figures for publication\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (nrow(mediation_summary) > 0) {
  
  # 6.1 Forest plot of indirect effects (English version)
  sig_mediation_eng <- mediation_summary %>%
    filter(all_paths_significant == TRUE | indirect_significant == TRUE) %>%
    mutate(
      pathway_label = paste0(exposure_label, " → ", mediator_label, 
                            " → ", outcome_label),
      or_indirect = exp(indirect_effect),
      or_lci = exp(indirect_effect - 1.96 * indirect_effect_se),
      or_uci = exp(indirect_effect + 1.96 * indirect_effect_se),
      significance = ifelse(significant_fdr, "FDR < 0.05", "P < 0.05")
    )
  
  if (nrow(sig_mediation_eng) > 0) {
    
    p_forest_eng <- ggplot(sig_mediation_eng, 
                          aes(x = or_indirect, 
                              y = reorder(pathway_label, or_indirect),
                              fill = strategy_label)) +
      geom_point(aes(shape = significance), size = 3) +
      geom_errorbarh(aes(xmin = or_lci, xmax = or_uci),
                    height = 0.3, alpha = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      scale_x_continuous(trans = "log10") +
      scale_shape_manual(values = c("FDR < 0.05" = 16, "P < 0.05" = 1)) +
      scale_fill_manual(
        values = c("Evidence-Based Selection" = "#4E79A7", 
                  "Exhaustive Screening" = "#F28E2B")
      ) +
      labs(
        title = "Mediation Analysis: Indirect Effects via Inflammatory Markers",
        subtitle = sprintf("%d significant mediation pathways", 
                          nrow(sig_mediation_eng)),
        x = "Odds Ratio (95% CI, Indirect Effect)",
        y = "Mediation Pathway",
        shape = "Significance Level",
        fill = "Selection Strategy"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.box = "vertical"
      )
    
    ggsave(
      "results/figures/step07b_english/mediation_forest_plot_english.png",
      p_forest_eng, 
      width = 14, 
      height = min(45, max(8, nrow(sig_mediation_eng) * 0.3)), 
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      "results/figures/step07b_english/mediation_forest_plot_english.pdf",
      p_forest_eng, 
      width = 14, 
      height = min(45, max(8, nrow(sig_mediation_eng) * 0.3)),
      device = "pdf",
      limitsize = FALSE
    )
    
    cat("✓ Forest plot (English version) saved\n")
  }
  
  # 6.2 Mediation proportion bar chart (English version)
  sig_prop_eng <- mediation_summary %>%
    filter(all_paths_significant == TRUE,
           !is.na(mediation_proportion_percent),
           partial_success == FALSE) %>%
    arrange(desc(abs(mediation_proportion_percent))) %>%
    head(20) %>%
    mutate(
      pathway_label = paste0(exposure_label, " → ", mediator_label, 
                            " → ", outcome_label)
    )
  
  if (nrow(sig_prop_eng) > 0) {
    
    p_proportion_eng <- ggplot(sig_prop_eng, 
                              aes(x = reorder(pathway_label, 
                                             mediation_proportion_percent),
                                  y = mediation_proportion_percent,
                                  fill = strategy_label)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      scale_fill_manual(
        values = c("Evidence-Based Selection" = "#4E79A7", 
                  "Exhaustive Screening" = "#F28E2B")
      ) +
      labs(
        title = "Proportion of Effect Mediated by Inflammatory Markers",
        subtitle = "Top 20 pathways with significant indirect effects",
        x = "Mediation Pathway",
        y = "Mediation Proportion (%)",
        fill = "Selection Strategy"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom"
      )
    
    ggsave(
      "results/figures/step07b_english/mediation_proportion_english.png",
      p_proportion_eng, 
      width = 12, 
      height = min(45, max(8, nrow(sig_prop_eng) * 0.4)),
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      "results/figures/step07b_english/mediation_proportion_english.pdf",
      p_proportion_eng, 
      width = 12, 
      height = min(45, max(8, nrow(sig_prop_eng) * 0.4)),
      device = "pdf",
      limitsize = FALSE
    )
    
    cat("✓ Mediation proportion chart (English version) saved\n")
  }
  
  # 6.3 Method comparison plot (if available)
  if (exists("method_comparison") && nrow(method_comparison) > 0) {
    
    # Select top 5 pathways for method comparison visualization
    n_pathways <- length(unique(method_comparison$pathway_id))
    top5_pathways <- if (n_pathways > 0) {
      unique(method_comparison$pathway_id)[seq_len(min(5, n_pathways))]
    } else {
      integer(0)
    }
    
    method_comp_plot <- method_comparison %>%
      filter(pathway_id %in% top5_pathways,
             step == "Exposure → Mediator") %>%
      mutate(
        pathway_label = paste0(exposure_label, " → ", mediator_label),
        ci_lower = b - 1.96 * se,
        ci_upper = b + 1.96 * se
      )
    
    if (nrow(method_comp_plot) > 0) {
      
      p_method_comp <- ggplot(method_comp_plot,
                             aes(x = b, 
                                 y = method_clean,
                                 color = method_clean)) +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                      height = 0.3, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        facet_wrap(~ pathway_label, ncol = 1, scales = "free_y") +
        scale_color_brewer(palette = "Set2") +
        labs(
          title = "MR Method Comparison: Top 5 Mediation Pathways",
          subtitle = "Exposure → Mediator associations",
          x = "Effect Estimate (β) with 95% CI",
          y = "MR Method",
          color = "Method"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "bottom",
          strip.text = element_text(size = 9, face = "bold")
        )
      
      ggsave(
        "results/figures/step07b_sensitivity/mr_method_comparison_plot.png",
        p_method_comp,
        width = 12,
        height = min(45, max(8, length(top5_pathways) * 2)),
        dpi = 300,
        limitsize = FALSE
      )
      
      cat("✓ MR method comparison plot saved\n")
    }
  }
  
  # 6.4 Heterogeneity and pleiotropy summary plot
  if (exists("sensitivity_results") && nrow(sensitivity_results) > 0) {
    
    # Prepare data for visualization
    het_pleio_summary <- sensitivity_results %>%
      select(pathway_id, exposure_label, mediator_label, outcome_label,
             heterogeneity_exp_to_med, pleiotropy_exp_to_med,
             heterogeneity_med_to_out, pleiotropy_med_to_out) %>%
      mutate(
        pathway_label = paste0(exposure_label, " → ", mediator_label),
        het_exp_med_binary = ifelse(heterogeneity_exp_to_med == "Significant", 1, 0),
        pleio_exp_med_binary = ifelse(pleiotropy_exp_to_med == "Detected", 1, 0),
        het_med_out_binary = ifelse(heterogeneity_med_to_out == "Significant", 1, 0),
        pleio_med_out_binary = ifelse(pleiotropy_med_to_out == "Detected", 1, 0)
      ) %>%
      head(15)  # Top 15 pathways
    
    if (nrow(het_pleio_summary) > 0) {
      
      # Create heatmap-style plot
      het_pleio_long <- het_pleio_summary %>%
        select(pathway_label, het_exp_med_binary, pleio_exp_med_binary,
               het_med_out_binary, pleio_med_out_binary) %>%
        pivot_longer(
          cols = c(het_exp_med_binary, pleio_exp_med_binary,
                  het_med_out_binary, pleio_med_out_binary),
          names_to = "test",
          values_to = "result"
        ) %>%
        mutate(
          test_label = case_when(
            test == "het_exp_med_binary" ~ "Heterogeneity\n(Exp → Med)",
            test == "pleio_exp_med_binary" ~ "Pleiotropy\n(Exp → Med)",
            test == "het_med_out_binary" ~ "Heterogeneity\n(Med → Out)",
            test == "pleio_med_out_binary" ~ "Pleiotropy\n(Med → Out)",
            TRUE ~ test
          ),
          result_label = ifelse(result == 1, "Yes", "No")
        )
      
      p_het_pleio <- ggplot(het_pleio_long,
                           aes(x = test_label, 
                               y = reorder(pathway_label, result),
                               fill = result_label)) +
        geom_tile(color = "white", size = 0.5) +
        scale_fill_manual(
          values = c("Yes" = "#D55E00", "No" = "#009E73"),
          name = "Detected"
        ) +
        labs(
          title = "Heterogeneity and Pleiotropy Tests for Mediation Pathways",
          subtitle = "Top 15 significant pathways",
          x = "Sensitivity Test",
          y = "Mediation Pathway"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
          legend.position = "right"
        )
      
      ggsave(
        "results/figures/step07b_sensitivity/heterogeneity_pleiotropy_heatmap.png",
        p_het_pleio,
        width = 10,
        height = min(45, max(8, nrow(het_pleio_summary) * 0.4)),
        dpi = 300,
        limitsize = FALSE
      )
      
      cat("✓ Heterogeneity and pleiotropy heatmap saved\n")
    }
  }
  
  cat("\n")
}

# ============================================================================
# Step 7: Generate comprehensive sensitivity analysis report
# ============================================================================

cat("【Step 7】Generate comprehensive sensitivity analysis report\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# Generate report
capture.output({
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Step 7b: Mediation Analysis - Sensitivity Analysis Report\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  cat("【Analysis Overview】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat(sprintf("Total pathways analyzed:        %d\n", nrow(mediation_summary)))
  cat(sprintf("Significant pathways (FDR<0.05): %d\n", 
             sum(mediation_summary$significant_fdr, na.rm = TRUE)))
  cat(sprintf("Pathways in sensitivity tests:  %d\n\n", nrow(sig_pathways)))
  
  if (exists("sensitivity_results") && nrow(sensitivity_results) > 0) {
    cat("【Heterogeneity Tests】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    # Exposure -> Mediator
    n_het_exp_med <- sum(sensitivity_results$heterogeneity_exp_to_med == "Significant", 
                        na.rm = TRUE)
    cat(sprintf("Exposure → Mediator: %d/%d pathways show significant heterogeneity\n",
               n_het_exp_med, nrow(sensitivity_results)))
    
    # Mediator -> Outcome
    n_het_med_out <- sum(sensitivity_results$heterogeneity_med_to_out == "Significant", 
                        na.rm = TRUE)
    cat(sprintf("Mediator → Outcome:  %d/%d pathways show significant heterogeneity\n\n",
               n_het_med_out, nrow(sensitivity_results)))
    
    cat("【Pleiotropy Tests (MR-Egger)】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    # Exposure -> Mediator
    n_pleio_exp_med <- sum(sensitivity_results$pleiotropy_exp_to_med == "Detected", 
                          na.rm = TRUE)
    cat(sprintf("Exposure → Mediator: %d/%d pathways show evidence of pleiotropy\n",
               n_pleio_exp_med, nrow(sensitivity_results)))
    
    # Mediator -> Outcome
    n_pleio_med_out <- sum(sensitivity_results$pleiotropy_med_to_out == "Detected", 
                          na.rm = TRUE)
    cat(sprintf("Mediator → Outcome:  %d/%d pathways show evidence of pleiotropy\n\n",
               n_pleio_med_out, nrow(sensitivity_results)))
    
    # Robust pathways (no heterogeneity or pleiotropy)
    robust_pathways <- sensitivity_results %>%
      filter(
        (is.na(heterogeneity_exp_to_med) | 
         heterogeneity_exp_to_med == "Not Significant") &
        (is.na(pleiotropy_exp_to_med) | 
         pleiotropy_exp_to_med == "Not Detected") &
        (is.na(heterogeneity_med_to_out) | 
         heterogeneity_med_to_out == "Not Significant") &
        (is.na(pleiotropy_med_to_out) | 
         pleiotropy_med_to_out == "Not Detected")
      )
    
    cat("【Robust Pathways】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("Pathways passing all sensitivity tests: %d/%d\n\n",
               nrow(robust_pathways), nrow(sensitivity_results)))
    
    if (nrow(robust_pathways) > 0) {
      cat("Top Robust Pathways:\n")
      for (i in seq_len(min(10, nrow(robust_pathways)))) {
        pathway <- robust_pathways[i, ]
        cat(sprintf("%d. %s → %s → %s\n",
                   i,
                   pathway$exposure_label,
                   pathway$mediator_label,
                   pathway$outcome_label))
      }
      cat("\n")
    }
  }
  
  if (exists("loo_combined") && nrow(loo_combined) > 0) {
    cat("【Leave-One-Out Analysis】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("Pathways with LOO analysis:     %d\n",
               length(unique(loo_combined$pathway_id))))
    cat(sprintf("Total LOO iterations performed: %d\n\n",
               nrow(loo_combined)))
    
    # Check stability of estimates
    loo_stability <- loo_combined %>%
      group_by(pathway_id, exposure_label, mediator_label, outcome_label) %>%
      summarise(
        mean_loo_beta = mean(loo_beta, na.rm = TRUE),
        sd_loo_beta = sd(loo_beta, na.rm = TRUE),
        original_beta = first(original_beta),
        max_deviation = max(abs(loo_beta - original_beta), na.rm = TRUE),
        relative_deviation = max_deviation / abs(original_beta),
        .groups = "drop"
      ) %>%
      mutate(
        stable = ifelse(relative_deviation < 0.2, "Stable", "Variable")
      )
    
    n_stable <- sum(loo_stability$stable == "Stable", na.rm = TRUE)
    cat(sprintf("Stable pathways (< 20%% deviation): %d/%d\n\n",
               n_stable, nrow(loo_stability)))
  }
  
  if (exists("method_comparison") && nrow(method_comparison) > 0) {
    cat("【MR Method Comparison】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    method_counts <- method_comparison %>%
      group_by(method_clean) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))
    
    cat("Methods used:\n")
    for (i in seq_len(nrow(method_counts))) {
      cat(sprintf("  - %-20s: %d estimates\n",
                 method_counts$method_clean[i],
                 method_counts$n[i]))
    }
    cat("\n")
    
    # Check consistency across methods
    method_consistency <- method_comparison %>%
      filter(method_clean %in% c("IVW", "Weighted Median", "MR-Egger")) %>%
      group_by(pathway_id, exposure_label, mediator_label, step) %>%
      summarise(
        n_methods = n(),
        all_same_direction = length(unique(sign(b))) == 1,
        .groups = "drop"
      ) %>%
      filter(n_methods >= 2)
    
    if (nrow(method_consistency) > 0) {
      n_consistent <- sum(method_consistency$all_same_direction, na.rm = TRUE)
      cat(sprintf("Pathways with consistent direction across methods: %d/%d\n\n",
                 n_consistent, nrow(method_consistency)))
    }
  }
  
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Report generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
}, file = "results/step07b_sensitivity_analysis_report.txt")

cat("✓ Comprehensive sensitivity analysis report saved\n")
cat("  → results/step07b_sensitivity_analysis_report.txt\n\n")

# ============================================================================
# Final summary
# ============================================================================

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Step 7b Completed!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【Files Generated】\n\n")

cat("Sensitivity Analysis Results:\n")
if (file.exists("results/tables/step07b_sensitivity/heterogeneity_pleiotropy_tests.xlsx")) {
  cat("  ✓ results/tables/step07b_sensitivity/heterogeneity_pleiotropy_tests.xlsx\n")
}
if (file.exists("results/tables/step07b_sensitivity/leave_one_out_results.xlsx")) {
  cat("  ✓ results/tables/step07b_sensitivity/leave_one_out_results.xlsx\n")
}
if (file.exists("results/tables/step07b_sensitivity/mr_method_comparison.xlsx")) {
  cat("  ✓ results/tables/step07b_sensitivity/mr_method_comparison.xlsx\n")
}

cat("\nEnglish Version Figures:\n")
if (file.exists("results/figures/step07b_english/mediation_forest_plot_english.png")) {
  cat("  ✓ results/figures/step07b_english/mediation_forest_plot_english.png\n")
  cat("  ✓ results/figures/step07b_english/mediation_forest_plot_english.pdf\n")
}
if (file.exists("results/figures/step07b_english/mediation_proportion_english.png")) {
  cat("  ✓ results/figures/step07b_english/mediation_proportion_english.png\n")
  cat("  ✓ results/figures/step07b_english/mediation_proportion_english.pdf\n")
}

cat("\nSensitivity Analysis Figures:\n")
if (file.exists("results/figures/step07b_sensitivity/leave_one_out_example.png")) {
  cat("  ✓ results/figures/step07b_sensitivity/leave_one_out_example.png\n")
}
if (file.exists("results/figures/step07b_sensitivity/mr_method_comparison_plot.png")) {
  cat("  ✓ results/figures/step07b_sensitivity/mr_method_comparison_plot.png\n")
}
if (file.exists("results/figures/step07b_sensitivity/heterogeneity_pleiotropy_heatmap.png")) {
  cat("  ✓ results/figures/step07b_sensitivity/heterogeneity_pleiotropy_heatmap.png\n")
}

cat("\nReport:\n")
if (file.exists("results/step07b_sensitivity_analysis_report.txt")) {
  cat("  ✓ results/step07b_sensitivity_analysis_report.txt\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Sensitivity analysis completed successfully!\n")
cat("All figures are now available in English for publication.\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

