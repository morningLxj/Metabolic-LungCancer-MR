# =============================================================================
# Data Processing Utility Functions for MR Analysis
# =============================================================================
# Author: morningLxj
# Date: 2024-12-01
# =============================================================================

#' Load and Process GWAS Summary Statistics
#'
#' @param file_path Path to GWAS summary statistics file
#' @param format File format (auto-detected from extension)
#' @param required_cols Required column names
#' @param verbose Whether to print progress messages
#' @return Processed GWAS data frame
#' @export
load_gwas_data <- function(file_path, format = "auto", 
                          required_cols = NULL, verbose = TRUE) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Auto-detect format
  if (format == "auto") {
    ext <- tolower(tools::file_ext(file_path))
    format <- switch(ext,
                    "txt" = "tab",
                    "tsv" = "tab", 
                    "csv" = "csv",
                    "gz" = "gzip",
                    "bgz" = "gzip",
                    stop("Unsupported file format: ", ext))
  }
  
  if (verbose) cat("Loading GWAS data from:", file_path, "\n")
  
  # Load data based on format
  if (format == "csv") {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
  } else if (format == "tab") {
    data <- read.table(file_path, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, quote = "")
  } else if (format == "gzip") {
    data <- read.table(gzfile(file_path), header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")
  } else {
    stop("Unsupported format: ", format)
  }
  
  if (verbose) cat("Loaded", nrow(data), "rows and", ncol(data), "columns\n")
  
  # Validate required columns if specified
  if (!is.null(required_cols)) {
    missing_cols <- required_cols[!required_cols %in% colnames(data)]
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  
  return(data)
}

#' Standardize Column Names for GWAS Data
#'
#' @param data GWAS data frame
#' @param mapping Column name mapping (optional)
#' @return Data frame with standardized column names
#' @export
standardize_columns <- function(data, mapping = NULL) {
  
  # Default column name mapping
  default_mapping <- list(
    "SNP" = c("snp", "rsid", "rs_id", "rs_number"),
    "CHR" = c("chr", "chromosome", "chrom"),
    "POS" = c("pos", "bp", "position"),
    "EA" = c("ea", "effect_allele", "allele1", "a1"),
    "OA" = c("oa", "other_allele", "allele2", "a2"),
    "EAF" = c("eaf", "effect_allele_freq", "freq", "maf"),
    "BETA" = c("beta", "effect", "or", "logor"),
    "SE" = c("se", "standard_error", "std_error"),
    "P" = c("p", "pval", "p_value", "pvalue"),
    "N" = c("n", "sample_size", "n_samples", "effective_n"),
    "NCASE" = c("ncase", "cases", "n_cases"),
    "NCTRL" = c("nctrl", "controls", "n_controls")
  )
  
  # Use custom mapping if provided
  if (!is.null(mapping)) {
    mapping <- c(default_mapping, mapping)
  } else {
    mapping <- default_mapping
  }
  
  # Find and rename columns
  new_names <- colnames(data)
  for (new_name in names(mapping)) {
    possible_names <- mapping[[new_name]]
    match_idx <- match(tolower(possible_names), tolower(colnames(data)))
    valid_idx <- !is.na(match_idx)
    
    if (any(valid_idx)) {
      old_name <- colnames(data)[match_idx[valid_idx][1]]
      new_names[new_names == old_name] <- new_name
      if (old_name != new_name) {
        cat("Renaming column:", old_name, "->", new_name, "\n")
      }
    }
  }
  
  colnames(data) <- new_names
  return(data)
}

#' Quality Control for GWAS Data
#'
#' @param data GWAS data frame
#' @param params QC parameters list
#' @return QC'd data frame
#' @export
quality_control <- function(data, params = NULL) {
  
  # Default QC parameters
  default_params <- list(
    pval_min = 0,
    pval_max = 1,
    maf_min = 0.01,
    info_min = 0.8,
    hwe_min = 1e-10,
    miss_rate_max = 0.05
  )
  
  # Merge with provided parameters
  if (!is.null(params)) {
    default_params[names(params)] <- params
  }
  
  cat("Starting quality control...\n")
  cat("Original data:", nrow(data), "variants\n")
  
  # Filter by P-value
  if ("P" %in% colnames(data)) {
    before <- nrow(data)
    data <- data[data$P >= default_params$pval_min & data$P <= default_params$pval_max, ]
    after <- nrow(data)
    cat("P-value filter: removed", before - after, "variants\n")
  }
  
  # Filter by MAF (if EAF available)
  if ("EAF" %in% colnames(data) && default_params$maf_min > 0) {
    before <- nrow(data)
    data <- data[data$EAF >= default_params$maf_min & data$EAF <= (1 - default_params$maf_min), ]
    after <- nrow(data)
    cat("MAF filter: removed", before - after, "variants\n")
  }
  
  # Filter by INFO score (if available)
  if ("INFO" %in% colnames(data)) {
    before <- nrow(data)
    data <- data[data$INFO >= default_params$info_min, ]
    after <- nrow(data)
    cat("INFO filter: removed", before - after, "variants\n")
  }
  
  # Filter by HWE (if available)
  if ("HWE" %in% colnames(data)) {
    before <- nrow(data)
    data <- data[data$HWE >= default_params$hwe_min, ]
    after <- nrow(data)
    cat("HWE filter: removed", before - after, "variants\n")
  }
  
  cat("Final data after QC:", nrow(data), "variants\n")
  return(data)
}

#' Harmonize Effect Alleles Between Exposure and Outcome
#'
#' @param exposure_data Exposure GWAS data
#' @param outcome_data Outcome GWAS data
#' @param action Action when alleles don't match ("flip", "drop", or "keep")
#' @return Harmonized data list
#' @export
harmonize_alleles <- function(exposure_data, outcome_data, action = "flip") {
  
  cat("Harmonizing effect alleles...\n")
  
  # Ensure SNP columns exist
  if (!"SNP" %in% colnames(exposure_data) || !"SNP" %in% colnames(outcome_data)) {
    stop("SNP column not found in both datasets")
  }
  
  # Merge on SNP
  harmonized <- merge(exposure_data, outcome_data, by = "SNP", suffixes = c(".exp", ".out"))
  cat("Merged on", nrow(harmonized), "SNPs\n")
  
  # Check and harmonize alleles
  if (all(c("EA.exp", "OA.exp", "EA.out", "OA.out") %in% colnames(harmonized))) {
    
    before <- nrow(harmonized)
    
    # Identify palindromic SNPs
    palindromic <- harmonized$EA.exp == harmonized$OA.out | 
                  (harmonized$EA.exp == toupper(harmonized$OA.out) & 
                   harmonized$OA.exp == toupper(harmonized$EA.out))
    
    # Handle mismatches
    mismatch <- harmonized$EA.exp != harmonized$EA.out
    
    if (action == "flip") {
      # Flip mismatched SNPs
      flip_idx <- mismatch & !palindromic
      harmonized$BETA.out[flip_idx] <- -harmonized$BETA.out[flip_idx]
      
      # Also flip EA/OA for outcome if needed
      temp_ea <- harmonized$EA.out[flip_idx]
      harmonized$EA.out[flip_idx] <- harmonized$OA.out[flip_idx]
      harmonized$OA.out[flip_idx] <- temp_ea
      
      if (sum(flip_idx) > 0) {
        cat("Flipped", sum(flip_idx), "SNPs\n")
      }
      
    } else if (action == "drop") {
      # Drop mismatched SNPs
      harmonized <- harmonized[!mismatch, ]
      cat("Dropped", before - nrow(harmonized), "mismatched SNPs\n")
    }
    
    # Handle palindromic SNPs
    pal_count <- sum(palindromic)
    if (pal_count > 0) {
      cat("Found", pal_count, "palindromic SNPs\n")
      
      if (action %in% c("flip", "drop")) {
        # Remove palindromic SNPs to avoid ambiguity
        harmonized <- harmonized[!palindromic, ]
        cat("Removed", pal_count, "palindromic SNPs\n")
      }
    }
  }
  
  cat("Final harmonized data:", nrow(harmonized), "SNPs\n")
  return(harmonized)
}

#' Remove Duplicate Variants
#'
#' @param data Data frame
#' @param by Columns to identify duplicates (default: SNP)
#' @param keep Keep "first", "last", or "none"
#' @return Data frame without duplicates
#' @export
remove_duplicates <- function(data, by = "SNP", keep = "first") {
  
  if (!by %in% colnames(data)) {
    stop("Duplicate identification column not found: ", by)
  }
  
  before <- nrow(data)
  
  if (keep == "first") {
    data <- data[!duplicated(data[, by]), ]
  } else if (keep == "last") {
    data <- data[!duplicated(data[, by], fromLast = TRUE), ]
  } else if (keep == "none") {
    # Remove all duplicates
    dup_idx <- duplicated(data[, by])
    data <- data[!dup_idx, ]
  }
  
  after <- nrow(data)
  cat("Removed", before - after, "duplicates\n")
  
  return(data)
}

#' Calculate Linkage Disequilibrium (LD) Matrix
#'
#' @param snp_list List of SNPs
#' @param population Population code (e.g., "EUR", "EAS")
#' @param r2_threshold R² threshold for clumping
#' @param kb_threshold Physical distance threshold in kb
#' @return Data frame with LD clumped results
#' @export
calculate_ld_matrix <- function(snp_list, population = "EUR", 
                               r2_threshold = 0.001, kb_threshold = 10000) {
  
  # This is a placeholder function
  # In practice, this would use external LD reference data
  # or call APIs like LDlink
  
  cat("Calculating LD for", length(snp_list), "SNPs in", population, "population\n")
  
  # Placeholder LD matrix
  n_snps <- length(snp_list)
  ld_matrix <- matrix(1, nrow = n_snps, ncol = n_snps,
                     dimnames = list(snp_list, snp_list))
  
  # Add some realistic LD structure (for demonstration)
  set.seed(123)
  for (i in 1:(n_snps-1)) {
    for (j in (i+1):n_snps) {
      # Simulate LD based on distance and random factors
      ld_value <- runif(1) * exp(-abs(i-j)/100)  # Decreases with distance
      ld_value <- min(ld_value, r2_threshold * 0.8)  # Keep below threshold
      ld_matrix[i,j] <- ld_matrix[j,i] <- ld_value
    }
  }
  
  return(ld_matrix)
}

#' Clump SNPs Based on LD
#'
#' @param data Input data frame
#' @param snp_col SNP column name
#' @param pval_col P-value column name
#' @param chr_col Chromosome column name (optional)
#' @param pos_col Position column name (optional)
#' @param r2_threshold R² threshold for clumping
#' @param kb_threshold Physical distance threshold in kb
#' @return Clumped data frame
#' @export
clump_snps <- function(data, snp_col = "SNP", pval_col = "P",
                      chr_col = NULL, pos_col = NULL,
                      r2_threshold = 0.001, kb_threshold = 10000) {
  
  cat("Clumping SNPs with R² <=", r2_threshold, "and distance >", kb_threshold, "kb\n")
  
  if (!all(c(snp_col, pval_col) %in% colnames(data))) {
    stop("Required columns not found: ", 
         paste(setdiff(c(snp_col, pval_col), colnames(data)), collapse = ", "))
  }
  
  # Sort by P-value
  data <- data[order(data[[pval_col]]), ]
  
  # Initialize clumped results
  clumped_snps <- character(0)
  remaining_snps <- data[[snp_col]]
  
  # Simple clumping algorithm (placeholder)
  while (length(remaining_snps) > 0) {
    lead_snp <- remaining_snps[1]  # Best P-value SNP
    clumped_snps <- c(clumped_snps, lead_snp)
    
    # Remove SNPs in LD with lead SNP
    # This is simplified - in practice would use actual LD data
    remaining_snps <- remaining_snps[-1]
  }
  
  # Filter to clumped SNPs
  clumped_data <- data[data[[snp_col]] %in% clumped_snps, ]
  
  cat("Clumped to", nrow(clumped_data), "independent SNPs\n")
  return(clumped_data)
}

#' Save Processed Data
#'
#' @param data Data to save
#' @param file_path Output file path
#' @param format Output format ("csv", "txt", "rds")
#' @param compress Whether to compress output
#' @return NULL
#' @export
save_processed_data <- function(data, file_path, format = "csv", compress = FALSE) {
  
  # Ensure output directory exists
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save based on format
  if (format == "csv") {
    write.csv(data, file_path, row.names = FALSE)
  } else if (format == "txt") {
    write.table(data, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (format == "rds") {
    saveRDS(data, file_path, compress = compress)
  } else if (format == "gz") {
    write.table(data, gzfile(file_path), sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    stop("Unsupported output format: ", format)
  }
  
  cat("Saved", nrow(data), "rows to", file_path, "\n")
}

#' Generate Data Summary Report
#'
#' @param data Input data
#' @param title Report title
#' @param output_file Output file path
#' @return NULL
#' @export
generate_data_summary <- function(data, title = "Data Summary", 
                                 output_file = "data_summary.txt") {
  
  summary_text <- paste0(
    title, "\n",
    paste(rep("=", nchar(title)), collapse = ""), "\n\n",
    "Date: ", as.character(Sys.time()), "\n\n",
    "Dataset Information:\n",
    "- Rows: ", nrow(data), "\n",
    "- Columns: ", ncol(data), "\n",
    "- Column names: ", paste(colnames(data), collapse = ", "), "\n\n",
    
    "Summary Statistics:\n"
  )
  
  # Add column-specific summaries
  for (col in colnames(data)) {
    if (is.numeric(data[[col]])) {
      summary_text <- paste0(summary_text, "\n", col, ":\n")
      summary_text <- paste0(summary_text, 
                            "  Min: ", min(data[[col]], na.rm = TRUE), "\n",
                            "  Max: ", max(data[[col]], na.rm = TRUE), "\n",
                            "  Mean: ", round(mean(data[[col]], na.rm = TRUE), 4), "\n",
                            "  Missing: ", sum(is.na(data[[col]])), "\n")
    }
  }
  
  # Write to file
  writeLines(summary_text, output_file)
  cat("Data summary written to", output_file, "\n")
}