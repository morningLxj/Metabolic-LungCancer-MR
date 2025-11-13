# =============================================================================
# GWAS数据质量控制和预处理工具函数
# Data Quality Control and Preprocessing Utilities for GWAS Analysis
# =============================================================================

#' GWAS数据质量控制综合分析
#' Comprehensive Quality Control Analysis for GWAS Data
#'
#' @param gwas_data data.frame GWAS数据框，包含SNP信息
#' @param exposure_name character 暴露变量名称
#' @param outcome_name character 结果变量名称  
#' @param sample_size integer 样本量
#' @param maf_threshold numeric MAF阈值 (默认0.01)
#' @param hwe_threshold numeric Hardy-Weinberg平衡阈值 (默认1e-6)
#' @param info_threshold numeric 信息量阈值 (默认0.8)
#' @param miss_threshold numeric 缺失率阈值 (默认0.05)
#'
#' @return list 包含QC结果和过滤统计的信息
#' @export
#'
#' @examples
#' \dontrun{
#' qc_results <- gwas_qc_analysis(gwas_data, "BMI", "Lung_Cancer", 500000)
#' }

gwas_qc_analysis <- function(gwas_data, 
                           exposure_name, 
                           outcome_name, 
                           sample_size,
                           maf_threshold = 0.01,
                           hwe_threshold = 1e-6,
                           info_threshold = 0.8,
                           miss_threshold = 0.05) {
  
  cat("开始GWAS数据质量控制分析...\n")
  cat("暴露变量:", exposure_name, "\n")
  cat("结果变量:", outcome_name, "\n")
  cat("样本量:", sample_size, "\n\n")
  
  # 初始化QC统计
  qc_stats <- list(
    original_snps = nrow(gwas_data),
    after_maf_filter = 0,
    after_hwe_filter = 0,
    after_info_filter = 0,
    after_miss_filter = 0,
    final_snps = 0,
    filtered_snps = 0
  )
  
  # 原始数据摘要
  cat("=== 原始数据摘要 ===\n")
  print_summary_stats(gwas_data)
  
  # 1. MAF过滤
  if ("maf" %in% names(gwas_data)) {
    before_maf <- nrow(gwas_data)
    gwas_data <- gwas_data[gwas_data$maf >= maf_threshold, ]
    after_maf <- nrow(gwas_data)
    qc_stats$after_maf_filter <- after_maf
    cat(sprintf("MAF过滤 (>=%.3f): %d -> %d SNPs (移除 %d, %.2f%%)\n", 
                maf_threshold, before_maf, after_maf, 
                before_maf - after_maf, 
                (before_maf - after_maf) / before_maf * 100))
  }
  
  # 2. Hardy-Weinberg平衡过滤
  if ("hwe_p" %in% names(gwas_data)) {
    before_hwe <- nrow(gwas_data)
    gwas_data <- gwas_data[gwas_data$hwe_p >= hwe_threshold, ]
    after_hwe <- nrow(gwas_data)
    qc_stats$after_hwe_filter <- after_hwe
    cat(sprintf("HWE过滤 (>=%.2e): %d -> %d SNPs (移除 %d, %.2f%%)\n", 
                hwe_threshold, before_hwe, after_hwe, 
                before_hwe - after_hwe, 
                (before_hwe - after_hwe) / before_hwe * 100))
  }
  
  # 3. 信息量过滤
  if ("info" %in% names(gwas_data)) {
    before_info <- nrow(gwas_data)
    gwas_data <- gwas_data[gwas_data$info >= info_threshold, ]
    after_info <- nrow(gwas_data)
    qc_stats$after_info_filter <- after_info
    cat(sprintf("信息量过滤 (>=%.2f): %d -> %d SNPs (移除 %d, %.2f%%)\n", 
                info_threshold, before_info, after_info, 
                before_info - after_info, 
                (before_info - after_info) / before_info * 100))
  }
  
  # 4. 缺失率过滤
  if ("missing_rate" %in% names(gwas_data)) {
    before_miss <- nrow(gwas_data)
    gwas_data <- gwas_data[gwas_data$missing_rate <= miss_threshold, ]
    after_miss <- nrow(gwas_data)
    qc_stats$after_miss_filter <- after_miss
    cat(sprintf("缺失率过滤 (<=%.3f): %d -> %d SNPs (移除 %d, %.2f%%)\n", 
                miss_threshold, before_miss, after_miss, 
                before_miss - after_miss, 
                (before_miss - after_miss) / before_miss * 100))
  }
  
  # 最终统计
  final_snps <- nrow(gwas_data)
  qc_stats$final_snps <- final_snps
  qc_stats$filtered_snps <- qc_stats$original_snps - final_snps
  
  cat("\n=== QC结果摘要 ===\n")
  cat(sprintf("原始SNPs: %d\n", qc_stats$original_snps))
  cat(sprintf("最终SNPs: %d\n", final_snps))
  cat(sprintf("过滤掉SNPs: %d (%.2f%%)\n", 
              qc_stats$filtered_snps, 
              qc_stats$filtered_snps / qc_stats$original_snps * 100))
  
  # 返回结果
  return(list(
    filtered_data = gwas_data,
    qc_stats = qc_stats,
    parameters = list(
      maf_threshold = maf_threshold,
      hwe_threshold = hwe_threshold,
      info_threshold = info_threshold,
      miss_threshold = miss_threshold
    )
  ))
}

#' 工具变量质量控制
#' Quality Control for Instrumental Variables
#'
#' @param iv_data data.frame 工具变量数据
#' @param r2_threshold numeric R²阈值 (默认0.001)
#' @param f_stat_threshold numeric F统计量阈值 (默认10)
#' @param pval_threshold numeric p值阈值 (默认5e-8)
#'
#' @return list 过滤后的工具变量数据
#' @export
#'
iv_quality_control <- function(iv_data, 
                              r2_threshold = 0.001,
                              f_stat_threshold = 10,
                              pval_threshold = 5e-8) {
  
  cat("开始工具变量质量控制...\n")
  
  original_ivs <- nrow(iv_data)
  
  # F统计量过滤
  if ("f_stat" %in% names(iv_data)) {
    before_f <- nrow(iv_data)
    iv_data <- iv_data[iv_data$f_stat >= f_stat_threshold, ]
    after_f <- nrow(iv_data)
    cat(sprintf("F统计量过滤 (>=%.1f): %d -> %d IVs (移除 %d)\n", 
                f_stat_threshold, before_f, after_f, before_f - after_f))
  }
  
  # R²过滤
  if ("r2" %in% names(iv_data)) {
    before_r2 <- nrow(iv_data)
    iv_data <- iv_data[iv_data$r2 >= r2_threshold, ]
    after_r2 <- nrow(iv_data)
    cat(sprintf("R²过滤 (>=%.4f): %d -> %d IVs (移除 %d)\n", 
                r2_threshold, before_r2, after_r2, before_r2 - after_r2))
  }
  
  # p值过滤
  if ("pval" %in% names(iv_data)) {
    before_p <- nrow(iv_data)
    iv_data <- iv_data[iv_data$pval <= pval_threshold, ]
    after_p <- nrow(iv_data)
    cat(sprintf("p值过滤 (<=%.2e): %d -> %d IVs (移除 %d)\n", 
                pval_threshold, before_p, after_p, before_p - after_p))
  }
  
  final_ivs <- nrow(iv_data)
  cat(sprintf("最终工具变量数量: %d (保留了 %.2f%%)\n", 
              final_ivs, final_ivs / original_ivs * 100))
  
  return(iv_data)
}

#' 数据标准化和协调
#' Data Standardization and Harmonization
#'
#' @param exp_data data.frame 暴露数据
#' @param out_data data.frame 结果数据
#' @param snp_col character SNP列名
#' @param effect_col character 效应列名
#' @param se_col character 标准误列名
#' @param pval_col character p值列名
#'
#' @return list 标准化后的数据列表
#' @export
#'
standardize_gwas_data <- function(exp_data, 
                                out_data,
                                snp_col = "SNP",
                                effect_col = "beta",
                                se_col = "se",
                                pval_col = "pval") {
  
  cat("开始GWAS数据标准化...\n")
  
  # 确保必要列存在
  required_cols <- c(snp_col, effect_col, se_col, pval_col)
  
  # 检查暴露数据
  missing_exp <- required_cols[!required_cols %in% names(exp_data)]
  if (length(missing_exp) > 0) {
    stop(sprintf("暴露数据缺少列: %s", paste(missing_exp, collapse = ", ")))
  }
  
  # 检查结果数据  
  missing_out <- required_cols[!required_cols %in% names(out_data)]
  if (length(missing_out) > 0) {
    stop(sprintf("结果数据缺少列: %s", paste(missing_out, collapse = ", ")))
  }
  
  # 重命名列以保持一致
  col_mapping <- setNames(required_cols, c("SNP", "beta", "se", "pval"))
  
  exp_std <- exp_data %>%
    rename(!!!col_mapping) %>%
    mutate(
      # 确保数据类型正确
      beta = as.numeric(beta),
      se = as.numeric(se),
      pval = as.numeric(pval),
      # 计算Z分数
      z_score = beta / se,
      # 标准化
      beta_std = (beta - mean(beta, na.rm = TRUE)) / sd(beta, na.rm = TRUE),
      se_std = se / sd(beta, na.rm = TRUE)
    )
  
  out_std <- out_data %>%
    rename(!!!col_mapping) %>%
    mutate(
      beta = as.numeric(beta),
      se = as.numeric(se),
      pval = as.numeric(pval),
      z_score = beta / se,
      beta_std = (beta - mean(beta, na.rm = TRUE)) / sd(beta, na.rm = TRUE),
      se_std = se / sd(beta, na.rm = TRUE)
    )
  
  cat(sprintf("暴露数据: %d SNPs\n", nrow(exp_std)))
  cat(sprintf("结果数据: %d SNPs\n", nrow(out_std)))
  
  return(list(
    exposure = exp_std,
    outcome = out_std
  ))
}

#' 生成数据质量报告
#' Generate Data Quality Report
#'
#' @param qc_results list QC结果列表
#' @param output_file character 输出文件名
#' @export
generate_qc_report <- function(qc_results, output_file = "data_qc_report.html") {
  
  # 创建报告内容
  report_content <- sprintf("
  <html>
  <head>
    <title>GWAS Data Quality Control Report</title>
    <style>
      body { font-family: Arial, sans-serif; margin: 40px; }
      .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
      .stats { background-color: #e8f4fd; padding: 15px; margin: 10px 0; border-radius: 5px; }
      .parameters { background-color: #fff2e8; padding: 15px; margin: 10px 0; border-radius: 5px; }
      table { border-collapse: collapse; width: 100%%; }
      th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
      th { background-color: #f2f2f2; }
    </style>
  </head>
  <body>
    <div class='header'>
      <h1>GWAS Data Quality Control Report</h1>
      <p>Generated on: %s</p>
    </div>
    
    <div class='stats'>
      <h2>Quality Control Statistics</h2>
      <table>
        <tr><th>Stage</th><th>SNPs Remaining</th><th>SNPs Removed</th><th>Removal Rate</th></tr>
        <tr><td>Original Data</td><td>%d</td><td>-</td><td>-</td></tr>
        <tr><td>After MAF Filter</td><td>%d</td><td>%d</td><td>%.2f%%</td></tr>
        <tr><td>After HWE Filter</td><td>%d</td><td>%d</td><td>%.2f%%</td></tr>
        <tr><td>After Info Filter</td><td>%d</td><td>%d</td><td>%.2f%%</td></tr>
        <tr><td>After Missing Filter</td><td>%d</td><td>%d</td><td>%.2f%%</td></tr>
        <tr><td><strong>Final Dataset</strong></td><td><strong>%d</strong></td><td><strong>%d</strong></td><td><strong>%.2f%%</strong></td></tr>
      </table>
    </div>
    
    <div class='parameters'>
      <h2>Quality Control Parameters</h2>
      <ul>
        <li>MAF Threshold: %.3f</li>
        <li>HWE p-value Threshold: %.2e</li>
        <li>INFO Score Threshold: %.2f</li>
        <li>Missing Rate Threshold: %.3f</li>
      </ul>
    </div>
  </body>
  </html>
  ", 
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  qc_results$qc_stats$original_snps,
  qc_results$qc_stats$after_maf_filter,
  qc_results$qc_stats$original_snps - qc_results$qc_stats$after_maf_filter,
  ifelse(qc_results$qc_stats$original_snps > 0, 
         (qc_results$qc_stats$original_snps - qc_results$qc_stats$after_maf_filter) / qc_results$qc_stats$original_snps * 100, 0),
  qc_results$qc_stats$after_hwe_filter,
  qc_results$qc_stats$after_maf_filter - qc_results$qc_stats$after_hwe_filter,
  ifelse(qc_results$qc_stats$after_maf_filter > 0,
         (qc_results$qc_stats$after_maf_filter - qc_results$qc_stats$after_hwe_filter) / qc_results$qc_stats$after_maf_filter * 100, 0),
  qc_results$qc_stats$after_info_filter,
  qc_results$qc_stats$after_hwe_filter - qc_results$qc_stats$after_info_filter,
  ifelse(qc_results$qc_stats$after_hwe_filter > 0,
         (qc_results$qc_stats$after_hwe_filter - qc_results$qc_stats$after_info_filter) / qc_results$qc_stats$after_hwe_filter * 100, 0),
  qc_results$qc_stats$after_miss_filter,
  qc_results$qc_stats$after_info_filter - qc_results$qc_stats$after_miss_filter,
  ifelse(qc_results$qc_stats$after_info_filter > 0,
         (qc_results$qc_stats$after_info_filter - qc_results$qc_stats$after_miss_filter) / qc_results$qc_stats$after_info_filter * 100, 0),
  qc_results$qc_stats$final_snps,
  qc_results$qc_stats$filtered_snps,
  ifelse(qc_results$qc_stats$original_snps > 0,
         qc_results$qc_stats$filtered_snps / qc_results$qc_stats$original_snps * 100, 0),
  qc_results$parameters$maf_threshold,
  qc_results$parameters$hwe_threshold,
  qc_results$parameters$info_threshold,
  qc_results$parameters$miss_threshold)
  
  # 保存报告
  writeLines(report_content, output_file)
  cat(sprintf("数据质量报告已保存至: %s\n", output_file))
  
  return(output_file)
}

# 辅助函数：打印数据摘要统计
print_summary_stats <- function(data) {
  if (nrow(data) == 0) {
    cat("数据为空\n")
    return()
  }
  
  cat(sprintf("总行数: %d\n", nrow(data)))
  cat(sprintf("总列数: %d\n", ncol(data)))
  
  # 如果有MAF列
  if ("maf" %in% names(data)) {
    maf_stats <- summary(data$maf, na.rm = TRUE)
    cat(sprintf("MAF统计: Min=%.4f, Q1=%.4f, Median=%.4f, Q3=%.4f, Max=%.4f\n",
                maf_stats[1], maf_stats[2], maf_stats[3], maf_stats[5], maf_stats[6]))
  }
  
  # 如果有p值列
  if ("pval" %in% names(data)) {
    sig_snps <- sum(data$pval < 5e-8, na.rm = TRUE)
    cat(sprintf("显著性SNPs (p<5e-8): %d (%.2f%%)\n", 
                sig_snps, sig_snps / nrow(data) * 100))
  }
  
  cat("\n")
}