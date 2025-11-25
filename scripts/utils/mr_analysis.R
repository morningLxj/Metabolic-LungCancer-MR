# =============================================================================
# MR Analysis Utility Functions for Metabolic Lung Cancer Study
# =============================================================================
# Author: morningLxj
# Date: 2024-12-01
# =============================================================================

#' Extract Instruments from GWAS Data
#'
#' @param exposure_data Exposure GWAS data
#' @param pval_threshold P-value threshold for instrument selection
#' @param clump_r2 LD R² threshold for clumping
#' @param clump_kb Physical distance threshold for clumping (kb)
#' @param fstat_threshold F-statistic threshold
#' @param verbose Whether to print progress
#' @return Instrument data frame
#' @export
extract_instruments <- function(exposure_data, pval_threshold = 5e-8,
                               clump_r2 = 0.001, clump_kb = 10000,
                               fstat_threshold = 10, verbose = TRUE) {
  
  if (verbose) {
    cat("Extracting instruments from", nrow(exposure_data), "variants\n")
  }
  
  # Filter by P-value
  if ("P" %in% colnames(exposure_data)) {
    instruments <- exposure_data[exposure_data$P < pval_threshold, ]
  } else {
    stop("P-value column not found in exposure data")
  }
  
  if (verbose) {
    cat("After P-value threshold (<", pval_threshold, "):", nrow(instruments), "variants\n")
  }
  
  # Check for required columns
  required_cols <- c("SNP", "BETA", "SE", "EA", "OA")
  missing_cols <- required_cols[!required_cols %in% colnames(instruments)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Remove variants with invalid effect sizes
  valid_beta <- !is.na(instruments$BETA) & instruments$BETA != 0
  valid_se <- !is.na(instruments$SE) & instruments$SE > 0
  instruments <- instruments[valid_beta & valid_se, ]
  
  if (verbose) {
    cat("After filtering invalid effect sizes:", nrow(instruments), "variants\n")
  }
  
  # Clump variants based on LD (simplified implementation)
  if (nrow(instruments) > 1) {
    # Sort by P-value (lowest first)
    instruments <- instruments[order(instruments$P), ]
    
    # Simple clumping: keep only variants not in high LD
    # This is a placeholder - real implementation would use LD reference data
    clumped_indices <- c(1)  # Keep the best variant
    
    for (i in 2:nrow(instruments)) {
      # Simple distance-based clumping (if position available)
      if (all(c("CHR", "POS") %in% colnames(instruments))) {
        same_chr <- instruments$CHR[i] == instruments$CHR[clumped_indices]
        if (any(same_chr)) {
          same_chr_idx <- clumped_indices[same_chr]
          distances <- abs(instruments$POS[i] - instruments$POS[same_chr_idx])
          min_distance <- min(distances, na.rm = TRUE) / 1000  # Convert to kb
          
          if (min_distance < clump_kb) {
            next  # Skip this variant (too close)
          }
        }
      }
      
      clumped_indices <- c(clumped_indices, i)
    }
    
    instruments <- instruments[clumped_indices, ]
    
    if (verbose) {
      cat("After clumping:", nrow(instruments), "independent instruments\n")
    }
  }
  
  # Add instrument strength information
  if ("EAF" %in% colnames(instruments) && "N" %in% colnames(instruments)) {
    # Calculate F-statistic
    instruments$F_stat <- (instruments$BETA / instruments$SE)^2
    
    # Flag weak instruments
    instruments$weak_instrument <- instruments$F_stat < fstat_threshold
    
    weak_count <- sum(instruments$weak_instrument)
    if (verbose && weak_count > 0) {
      cat("Warning:", weak_count, "instruments have F-stat <", fstat_threshold, "\n")
    }
  }
  
  # Remove weak instruments if needed
  if ("F_stat" %in% colnames(instruments)) {
    before <- nrow(instruments)
    instruments <- instruments[!instruments$weak_instrument, ]
    removed <- before - nrow(instruments)
    if (removed > 0 && verbose) {
      cat("Removed", removed, "weak instruments\n")
    }
  }
  
  if (verbose) {
    cat("Final instruments:", nrow(instruments), "\n")
  }
  
  return(instruments)
}

#' Run Comprehensive Two-Sample MR Analysis
#'
#' @param instruments Instrument data frame
#' @param outcome_data Outcome GWAS data
#' @param methods MR methods to use
#' @param harmonize Whether to harmonize alleles
#' @param verbose Whether to print progress
#' @return Combined MR results
#' @export
run_mr_analysis <- function(instruments, outcome_data, 
                           methods = c("mr_ivw_fixed", "mr_egger_regression", 
                                     "mr_weighted_median", "mr_simple_mode"),
                           harmonize = TRUE, verbose = TRUE) {
  
  if (verbose) cat("Starting MR analysis...\n")
  
  # Harmonize data if requested
  if (harmonize) {
    mr_data <- harmonize_mr_data(instruments, outcome_data, verbose = verbose)
  } else {
    mr_data <- list(exposure = instruments, outcome = outcome_data)
  }
  
  # Prepare data for MR analysis
  # Convert to TwoSampleMR format
  exposure_mr <- instruments
  outcome_mr <- outcome_data
  
  # Ensure required columns
  if (!all(c("SNP", "BETA", "SE", "EA", "OA") %in% colnames(exposure_mr))) {
    stop("Missing required columns in exposure data")
  }
  if (!all(c("SNP", "BETA", "SE", "EA", "OA") %in% colnames(outcome_mr))) {
    stop("Missing required columns in outcome data")
  }
  
  # Combine data (simplified harmonization)
  harmonized_data <- merge(exposure_mr, outcome_mr, by = "SNP", suffixes = c(".exp", ".out"))
  
  # Flip outcome effect sizes for mismatched alleles
  flip_condition <- harmonized_data$EA.exp != harmonized_data$EA.out
  harmonized_data$BETA.out[flip_condition] <- -harmonized_data$BETA.out[flip_condition]
  
  if (verbose) {
    cat("Harmonized", nrow(harmonized_data), "SNPs\n")
  }
  
  # Run different MR methods
  all_results <- list()
  
  # Inverse Variance Weighted (Fixed Effects)
  if ("mr_ivw_fixed" %in% methods) {
    if (verbose) cat("Running IVW fixed effects...\n")
    ivw_result <- run_ivw_fixed(harmonized_data)
    all_results$ivw_fixed <- ivw_result
  }
  
  # MR-Egger
  if ("mr_egger_regression" %in% methods) {
    if (verbose) cat("Running MR-Egger regression...\n")
    egger_result <- run_egger_regression(harmonized_data)
    all_results$egger <- egger_result
  }
  
  # Weighted Median
  if ("mr_weighted_median" %in% methods) {
    if (verbose) cat("Running weighted median MR...\n")
    median_result <- run_weighted_median(harmonized_data)
    all_results$weighted_median <- median_result
  }
  
  # Simple Mode
  if ("mr_simple_mode" %in% methods) {
    if (verbose) cat("Running simple mode MR...\n")
    mode_result <- run_simple_mode(harmonized_data)
    all_results$simple_mode <- mode_result
  }
  
  # Combine all results
  combined_results <- do.call(rbind, lapply(names(all_results), function(method) {
    result <- all_results[[method]]
    result$method <- method
    return(result)
  }))
  
  if (verbose) {
    cat("MR analysis completed. Total results:", nrow(combined_results), "\n")
  }
  
  return(combined_results)
}

#' Inverse Variance Weighted Fixed Effects MR
#'
#' @param data Harmonized MR data
#' @return MR result
#' @keywords internal
run_ivw_fixed <- function(data) {
  # Calculate weights
  weights <- 1 / (data$SE.out^2)
  
  # Calculate weighted ratio estimates
  ratio_estimates <- data$BETA.out / data$BETA.exp
  
  # Weighted variance of ratio estimates
  ratio_var <- data$SE.out^2 / (data$BETA.exp^2)
  
  # IVW estimate
  ivw_beta <- sum(weights * ratio_estimates) / sum(weights)
  ivw_se <- sqrt(1 / sum(weights))
  ivw_pval <- 2 * pnorm(-abs(ivw_beta / ivw_se))
  
  # Calculate Q-statistic for heterogeneity
  q_stat <- sum(weights * (ratio_estimates - ivw_beta)^2)
  q_df <- length(ratio_estimates) - 1
  q_pval <- 1 - pchisq(q_stat, q_df)
  
  result <- data.frame(
    method = "IVW (fixed effects)",
    nsnp = nrow(data),
    exposure = "metabolic_trait",
    outcome = "lung_cancer",
    b = ivw_beta,
    se = ivw_se,
    pval = ivw_pval,
    or = exp(ivw_beta),
    or_lci95 = exp(ivw_beta - 1.96 * ivw_se),
    or_uci95 = exp(ivw_beta + 1.96 * ivw_se),
    q_statistic = q_stat,
    q_pval = q_pval,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' MR-Egger Regression
#'
#' @param data Harmonized MR data
#' @return MR result
#' @keywords internal
run_egger_regression <- function(data) {
  # Prepare data
  ratio_estimates <- data$BETA.out / data$BETA.exp
  ratio_se <- data$SE.out / abs(data$BETA.exp)
  
  # Add intercept term (reweighted by 1/(SE^2))
  weights <- 1 / (ratio_se^2)
  x_weighted <- data$BETA.exp * sqrt(weights)
  y_weighted <- ratio_estimates * sqrt(weights)
  
  # Eger regression: y = a + b*x + error
  x_mean <- sum(x_weighted) / sum(weights)
  y_mean <- sum(y_weighted) / sum(weights)
  
  cov_xy <- sum(weights * (data$BETA.exp - x_mean) * (ratio_estimates - y_mean))
  var_x <- sum(weights * (data$BETA.exp - x_mean)^2)
  
  # Eger estimates
  egger_b <- cov_xy / var_x
  intercept <- y_mean - egger_b * x_mean
  
  # Standard errors
  resid_var <- sum(weights * (ratio_estimates - intercept - egger_b * data$BETA.exp)^2) / 
               (length(ratio_estimates) - 2)
  se_b <- sqrt(resid_var / var_x)
  se_intercept <- sqrt(resid_var / sum(weights))
  
  # P-values
  b_pval <- 2 * pnorm(-abs(egger_b / se_b))
  intercept_pval <- 2 * pnorm(-abs(intercept / se_intercept))
  
  result <- data.frame(
    method = "MR-Egger",
    nsnp = nrow(data),
    exposure = "metabolic_trait", 
    outcome = "lung_cancer",
    b = egger_b,
    se = se_b,
    pval = b_pval,
    or = exp(egger_b),
    or_lci95 = exp(egger_b - 1.96 * se_b),
    or_uci95 = exp(egger_b + 1.96 * se_b),
    intercept = intercept,
    intercept_se = se_intercept,
    intercept_pval = intercept_pval,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Weighted Median MR
#'
#' @param data Harmonized MR data
#' @return MR result
#' @keywords internal
run_weighted_median <- function(data) {
  # Calculate ratio estimates
  ratio_estimates <- data$BETA.out / data$BETA.exp
  ratio_se <- data$SE.out / abs(data$BETA.exp)
  weights <- 1 / (ratio_se^2)
  
  # Weighted median using bootstrap
  set.seed(123)
  n_boot <- 10000
  boot_estimates <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    # Bootstrap sample
    boot_idx <- sample(length(ratio_estimates), replace = TRUE)
    boot_estimates[i] <- weighted.median(ratio_estimates[boot_idx], 
                                        weights[boot_idx])
  }
  
  # Calculate median and confidence intervals
  median_beta <- median(ratio_estimates, weights = weights)
  ci_lower <- quantile(boot_estimates, 0.025)
  ci_upper <- quantile(boot_estimates, 0.975)
  
  # Standard error from bootstrap
  se_beta <- sd(boot_estimates)
  pval <- 2 * pnorm(-abs(median_beta / se_beta))
  
  result <- data.frame(
    method = "Weighted Median",
    nsnp = nrow(data),
    exposure = "metabolic_trait",
    outcome = "lung_cancer", 
    b = median_beta,
    se = se_beta,
    pval = pval,
    or = exp(median_beta),
    or_lci95 = exp(ci_lower),
    or_uci95 = exp(ci_upper),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Simple Mode MR
#'
#' @param data Harmonized MR data
#' @return MR result
#' @keywords internal
run_simple_mode <- function(data) {
  # Calculate ratio estimates
  ratio_estimates <- data$BETA.out / data$BETA.exp
  ratio_se <- data$SE.out / abs(data$BETA.exp)
  weights <- 1 / (ratio_se^2)
  
  # Simple mode (unweighted)
  mode_result <- weighted.median(ratio_estimates, weights = weights)
  
  # Calculate standard error using asymptotic approximation
  se_beta <- sqrt(sum(weights^2 * ratio_se^2) / (sum(weights))^2)
  pval <- 2 * pnorm(-abs(mode_result / se_beta))
  
  result <- data.frame(
    method = "Simple Mode",
    nsnp = nrow(data),
    exposure = "metabolic_trait",
    outcome = "lung_cancer",
    b = mode_result,
    se = se_beta,
    pval = pval,
    or = exp(mode_result),
    or_lci95 = exp(mode_result - 1.96 * se_beta),
    or_uci95 = exp(mode_result + 1.96 * se_beta),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Run Sensitivity Analysis
#'
#' @param mr_results MR analysis results
#' @param mr_data Harmonized MR data
#' @param methods Sensitivity analysis methods
#' @return Sensitivity analysis results
#' @export
run_sensitivity_analysis <- function(mr_results, mr_data, 
                                   methods = c("leave_one_out", "mr_presso", 
                                             "pleiotropy_test", "heterogeneity_test")) {
  
  sensitivity_results <- list()
  
  # Leave-one-out analysis
  if ("leave_one_out" %in% methods) {
    sensitivity_results$leave_one_out <- run_leave_one_out(mr_data)
  }
  
  # MR-PRESSO (simplified)
  if ("mr_presso" %in% methods) {
    sensitivity_results$mr_presso <- run_mr_presso(mr_data)
  }
  
  # Pleiotropy test
  if ("pleiotropy_test" %in% methods) {
    sensitivity_results$pleiotropy <- test_horizontal_pleiotropy(mr_data)
  }
  
  # Heterogeneity test
  if ("heterogeneity_test" %in% methods) {
    sensitivity_results$heterogeneity <- test_heterogeneity(mr_data)
  }
  
  return(sensitivity_results)
}

#' Leave-One-Out Analysis
#'
#' @param data Harmonized MR data
#' @return Leave-one-out results
#' @keywords internal
run_leave_one_out <- function(data) {
  n_snps <- nrow(data)
  loo_results <- data.frame(
    excluded_snp = character(n_snps),
    beta = numeric(n_snps),
    se = numeric(n_snps),
    pval = numeric(n_snps),
    stringsAsFactors = FALSE
  )
  
  # Calculate total effect
  weights <- 1 / (data$SE.out^2)
  ratio_estimates <- data$BETA.out / data$BETA.exp
  total_effect <- sum(weights * ratio_estimates) / sum(weights)
  
  for (i in 1:n_snps) {
    loo_data <- data[-i, ]
    if (nrow(loo_data) > 1) {
      loo_weights <- 1 / (loo_data$SE.out^2)
      loo_ratio <- loo_data$BETA.out / loo_data$BETA.exp
      loo_beta <- sum(loo_weights * loo_ratio) / sum(loo_weights)
      loo_se <- sqrt(1 / sum(loo_weights))
      loo_pval <- 2 * pnorm(-abs(loo_beta / loo_se))
      
      loo_results$excluded_snp[i] <- as.character(data$SNP[i])
      loo_results$beta[i] <- loo_beta
      loo_results$se[i] <- loo_se
      loo_results$pval[i] <- loo_pval
    }
  }
  
  return(loo_results)
}

#' MR-PRESSO Analysis (Simplified)
#'
#' @param data Harmonized MR data
#' @return MR-PRESSO results
#' @keywords internal
run_mr_presso <- function(data) {
  # Simplified MR-PRESSO implementation
  # In practice, this would use the MRPRESSO package
  
  ratio_estimates <- data$BETA.out / data$BETA.exp
  weights <- 1 / (data$SE.out^2)
  
  # Global test
  ivw_beta <- sum(weights * ratio_estimates) / sum(weights)
  q_stat <- sum(weights * (ratio_estimates - ivw_beta)^2)
  q_df <- length(ratio_estimates) - 1
  q_pval <- 1 - pchisq(q_stat, q_df)
  
 presso_result <- data.frame(
    test = "Global",
    statistic = q_stat,
    pval = q_pval,
    distortion = NA,
    stringsAsFactors = FALSE
  )
  
  return(presso_result)
}

#' Test Horizontal Pleiotropy
#'
#' @param data Harmonized MR data
#' @return Pleiotropy test results
#' @keywords internal
test_horizontal_pleiotropy <- function(data) {
  # Use MR-Egger intercept as test for pleiotropy
  ratio_estimates <- data$BETA.out / data$BETA.exp
  weights <- 1 / (data$SE.out / abs(data$BETA.exp))^2
  
  # Eger regression for intercept
  x_weighted <- data$BETA.exp * sqrt(weights)
  y_weighted <- ratio_estimates * sqrt(weights)
  
  weights_total <- sum(weights)
  x_mean <- sum(x_weighted) / weights_total
  y_mean <- sum(y_weighted) / weights_total
  
  cov_xy <- sum(weights * (data$BETA.exp - x_mean) * (ratio_estimates - y_mean))
  var_x <- sum(weights * (data$BETA.exp - x_mean)^2)
  
  intercept <- y_mean - (cov_xy / var_x) * x_mean
  
  # Standard error of intercept
  resid_var <- sum(weights * (ratio_estimates - intercept - (cov_xy / var_x) * data$BETA.exp)^2) / 
               (length(ratio_estimates) - 2)
  se_intercept <- sqrt(resid_var / sum(weights))
  
  # P-value for intercept
  intercept_pval <- 2 * pnorm(-abs(intercept / se_intercept))
  
  pleiotropy_result <- data.frame(
    method = "MR-Egger intercept",
    intercept = intercept,
    se = se_intercept,
    pval = intercept_pval,
    interpretation = ifelse(intercept_pval < 0.05, 
                           "Evidence of horizontal pleiotropy", 
                           "No evidence of horizontal pleiotropy"),
    stringsAsFactors = FALSE
  )
  
  return(pleiotropy_result)
}

#' Test Heterogeneity
#'
#' @param data Harmonized MR data
#' @return Heterogeneity test results
#' @keywords internal
test_heterogeneity <- function(data) {
  ratio_estimates <- data$BETA.out / data$BETA.exp
  weights <- 1 / (data$SE.out^2)
  
  # IVW estimate
  ivw_beta <- sum(weights * ratio_estimates) / sum(weights)
  
  # Q-statistic
  q_stat <- sum(weights * (ratio_estimates - ivw_beta)^2)
  q_df <- length(ratio_estimates) - 1
  q_pval <- 1 - pchisq(q_stat, q_df)
  
  heterogeneity_result <- data.frame(
    method = "Cochran's Q",
    statistic = q_stat,
    df = q_df,
    pval = q_pval,
    interpretation = ifelse(q_pval < 0.05, 
                           "Evidence of heterogeneity", 
                           "No evidence of heterogeneity"),
    stringsAsFactors = FALSE
  )
  
  return(heterogeneity_result)
}

#' Harmonize MR Data
#'
#' @param exposure_data Exposure data
#' @param outcome_data Outcome data
#' @param verbose Whether to print progress
#' @return Harmonized MR data
#' @export
harmonize_mr_data <- function(exposure_data, outcome_data, verbose = TRUE) {
  
  if (verbose) cat("Harmonizing MR data...\n")
  
  # Check required columns
  exp_cols <- c("SNP", "BETA", "SE", "EA", "OA")
  out_cols <- c("SNP", "BETA", "SE", "EA", "OA")
  
  if (!all(exp_cols %in% colnames(exposure_data))) {
    stop("Missing required columns in exposure data")
  }
  if (!all(out_cols %in% colnames(outcome_data))) {
    stop("Missing required columns in outcome data")
  }
  
  # Merge on SNP
  harmonized <- merge(exposure_data, outcome_data, by = "SNP", suffixes = c(".exp", ".out"))
  
  if (verbose) cat("Merged", nrow(harmonized), "SNPs\n")
  
  # Check allele compatibility
  if (verbose) {
    same_direction <- sum(harmonized$EA.exp == harmonized$EA.out)
    cat("SNPs with same effect allele:", same_direction, "out of", nrow(harmonized), "\n")
  }
  
  # Remove palindromic SNPs
  palindromic <- harmonized$EA.exp == harmonized$OA.out | 
                (toupper(harmonized$EA.exp) == toupper(harmonized$OA.out) & 
                 toupper(harmonized$OA.exp) == toupper(harmonized$EA.out))
  
  if (sum(palindromic) > 0) {
    if (verbose) cat("Removing", sum(palindromic), "palindromic SNPs\n")
    harmonized <- harmonized[!palindromic, ]
  }
  
  # Flip outcome effect sizes for mismatched alleles
  mismatched <- harmonized$EA.exp != harmonized$EA.out
  n_flipped <- sum(mismatched)
  
  if (n_flipped > 0) {
    if (verbose) cat("Flipping", n_flipped, "outcome effect sizes\n")
    harmonized$BETA.out[mismatched] <- -harmonized$BETA.out[mismatched]
    
    # Also swap outcome alleles
    temp_ea <- harmonized$EA.out[mismatched]
    harmonized$EA.out[mismatched] <- harmonized$OA.out[mismatched]
    harmonized$OA.out[mismatched] <- temp_ea
  }
  
  if (verbose) cat("Final harmonized data:", nrow(harmonized), "SNPs\n")
  
  return(harmonized)
}

#' Generate MR Results Summary
#'
#' @param mr_results MR analysis results
#' @param sensitivity_results Sensitivity analysis results
#' @param output_dir Output directory
#' @return NULL
#' @export
generate_mr_summary <- function(mr_results, sensitivity_results = NULL, 
                               output_dir = "results/mr_results/") {
  
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save main MR results
  write.csv(mr_results, file.path(output_dir, "mr_results.csv"), row.names = FALSE)
  
  # Generate summary text
  summary_text <- paste0(
    "# MR Analysis Summary\n\n",
    "## Main Results\n\n",
    "| Method | Beta | SE | P-value | OR | 95% CI |\n",
    "|--------|------|----|---------|----|---------|\n"
  )
  
  for (i in 1:nrow(mr_results)) {
    row <- mr_results[i, ]
    summary_text <- paste0(summary_text,
                          "| ", row$method, " | ",
                          sprintf("%.4f", row$b), " | ",
                          sprintf("%.4f", row$se), " | ",
                          sprintf("%.2e", row$pval), " | ",
                          sprintf("%.3f", row$or), " | ",
                          sprintf("%.3f-%.3f", row$or_lci95, row$or_uci95), " |\n")
  }
  
  # Add sensitivity results if available
  if (!is.null(sensitivity_results)) {
    summary_text <- paste0(summary_text, "\n## Sensitivity Analysis\n\n")
    
    if (!is.null(sensitivity_results$pleiotropy)) {
      summary_text <- paste0(summary_text, "### Pleiotropy Test\n\n")
      pleio <- sensitivity_results$pleiotropy
      summary_text <- paste0(summary_text,
                            "- Intercept: ", sprintf("%.4f", pleio$intercept), "\n",
                            "- SE: ", sprintf("%.4f", pleio$se), "\n",
                            "- P-value: ", sprintf("%.2e", pleio$pval), "\n",
                            "- Interpretation: ", pleio$interpretation, "\n\n")
    }
    
    if (!is.null(sensitivity_results$heterogeneity)) {
      summary_text <- paste0(summary_text, "### Heterogeneity Test\n\n")
      het <- sensitivity_results$heterogeneity
      summary_text <- paste0(summary_text,
                            "- Q-statistic: ", sprintf("%.4f", het$statistic), "\n",
                            "- P-value: ", sprintf("%.2e", het$pval), "\n",
                            "- Interpretation: ", het$interpretation, "\n\n")
    }
  }
  
  summary_text <- paste0(summary_text, "---\nGenerated on: ", as.character(Sys.time()), "\n")
  
  # Save summary
  writeLines(summary_text, file.path(output_dir, "mr_summary.md"))
  
  cat("MR summary generated and saved to", output_dir, "\n")
}