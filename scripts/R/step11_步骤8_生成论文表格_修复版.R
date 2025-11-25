############################################################################
# 步骤8：生成增强版论文表格（修复版）
# 用途：从已保存的分析结果生成所有论文表格
# 修复：解决"更换参数长度为零"的错误
############################################################################

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤8：生成增强版论文表格（修复版）\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ============================================================================
# 步骤1：环境准备与包加载
# ============================================================================

cat("【步骤1】环境准备与包加载\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 加载必要的R包
required_packages <- c("TwoSampleMR", "dplyr", "readr", "ggplot2", 
                     "openxlsx", "gridExtra", "tidyr", "stringr", 
                     "RColorBrewer", "pheatmap")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装包: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(openxlsx)
  library(gridExtra)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
})

cat("✓ 所有必需的R包已加载\n\n")

# 创建输出目录
dir.create("results/tables/paper_tables", showWarnings = FALSE, recursive = TRUE)
cat("✓ 输出目录已创建\n\n")

# ============================================================================
# 步骤2：加载已保存的分析结果
# ============================================================================

cat("【步骤2】加载已保存的分析结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

results_file <- "data/step11_reverse_mr_complete_results.RData"

if (!file.exists(results_file)) {
  stop(sprintf("\n错误：未找到保存的分析结果文件！\n请先运行 step11_反向MR分析_优化版.R 进行完整分析。\n文件路径: %s\n", results_file))
}

cat(sprintf("正在加载: %s\n", results_file))
load(results_file)

# 检查必需的对象
if (!exists("all_results")) {
  stop("错误：all_results 对象不存在！")
}

if (!exists("analysis_log")) {
  stop("错误：analysis_log 对象不存在！")
}

if (!exists("outcome_category")) {
  stop("错误：outcome_category 对象不存在！")
}

# 提取成功的结果
successful_results <- Filter(function(x) x$status == "Success", all_results)

if (length(successful_results) == 0) {
  stop("错误：没有成功的分析结果，无法生成表格。")
}

cat(sprintf("✓ 成功加载 %d 个分析结果\n\n", length(successful_results)))

# ============================================================================
# 步骤3：定义安全函数（修复ifelse错误）
# ============================================================================

cat("【步骤3】定义安全函数\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# 安全提取函数（避免"更换参数长度为零"错误）
safe_extract <- function(x, idx = 1, default = NA) {
  tryCatch({
    if (is.null(x)) {
      return(default)
    }
    # 处理长度为0的向量
    if (length(x) == 0) {
      return(default)
    }
    # 处理索引超出范围
    if (is.null(idx) || idx < 1 || idx > length(x)) {
      return(default)
    }
    # 提取值
    val <- x[idx]
    # 处理NULL
    if (is.null(val)) {
      return(default)
    }
    # 处理长度为0
    if (length(val) == 0) {
      return(default)
    }
    # 如果值是向量，取第一个元素
    if (length(val) > 1) {
      val <- val[1]
    }
    return(val)
  }, error = function(e) {
    return(default)
  })
}

# 安全的ifelse函数（修复"更换参数长度为零"错误）
safe_ifelse <- function(test, yes, no) {
  # 确保test是长度为1的逻辑值
  if (length(test) == 0) {
    return(no)
  }
  if (length(test) > 1) {
    test <- test[1]
  }
  if (is.na(test)) {
    return(no)
  }
  if (test) {
    # 确保yes是标量
    if (length(yes) == 0 || is.null(yes)) {
      return(no)
    }
    if (length(yes) > 1) {
      return(yes[1])
    }
    return(yes)
  } else {
    # 确保no是标量
    if (length(no) == 0 || is.null(no)) {
      return(NA)
    }
    if (length(no) > 1) {
      return(no[1])
    }
    return(no)
  }
}

# 安全获取策略字符串
safe_get_strategy <- function(extraction_info) {
  if (is.null(extraction_info) || nrow(extraction_info) == 0) {
    return("Unknown")
  }
  # analysis_log 使用 extraction_strategy 列
  strategy_val <- extraction_info$extraction_strategy[1]
  if (is.null(strategy_val) || is.na(strategy_val) || length(strategy_val) == 0) {
    return("Unknown")
  }
  if (is.na(strategy_val)) {
    return("Unknown")
  }
  return(as.character(strategy_val[1]))
}

# 安全获取n_snps
safe_get_nsnp <- function(nsnp_val, default_val) {
  if (is.null(nsnp_val) || length(nsnp_val) == 0 || is.na(nsnp_val)) {
    if (is.null(default_val) || length(default_val) == 0 || is.na(default_val)) {
      return(NA)
    }
    return(default_val)
  }
  return(nsnp_val)
}

cat("✓ 安全函数已定义\n\n")

# ============================================================================
# 步骤4：生成论文表格
# ============================================================================

cat("【步骤4】生成论文表格\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# Table S1: 完整MR结果（所有方法）
table_s1_data <- list()
# Table 2增强版: 主MR结果（包含所有方法列）
table2_enhanced_data <- list()
# Table 5: 敏感性分析汇总
table5_data <- list()
# Table S2: 敏感性分析详情
table_s2_data <- list()
# Table S3: Leave-one-out详情
table_s3_data <- list()
# Table S6: 工具变量强度
table_s6_data <- list()

cat("正在处理分析结果...\n")
processed_count <- 0

for (key in names(successful_results)) {
  result <- successful_results[[key]]
  if (is.null(result) || is.null(result$mr_results) || nrow(result$mr_results) == 0) {
    next
  }
  
  tryCatch({
    exposure <- safe_extract(result$exposure_name, 1, "Unknown")
    outcome <- safe_extract(result$outcome_name, 1, "Unknown")
    
    # 安全获取category
    category <- outcome_category[[outcome]]
    if (is.null(category) || length(category) == 0 || is.na(category)) {
      category <- "Unknown"
    }
    
    # 安全获取策略
    extraction_info <- analysis_log[analysis_log$exposure == exposure & 
                                   analysis_log$outcome == outcome, ]
    strategy <- safe_get_strategy(extraction_info)
    
    # 提取所有MR方法结果
    if (nrow(result$mr_results) > 0) {
      for (i in seq_len(nrow(result$mr_results))) {
        method_row <- result$mr_results[i, ]
        method_name <- safe_extract(method_row$method, 1, "Unknown")
        
        beta <- safe_extract(method_row$b, 1, NA)
        se <- safe_extract(method_row$se, 1, NA)
        pval <- safe_extract(method_row$pval, 1, NA)
        nsnp_val <- safe_extract(method_row$nsnp, 1, NA)
        
        # 安全获取n_snps_harmonized
        default_nsnp <- safe_extract(result$n_snps_harmonized, 1, NA)
        nsnp <- safe_get_nsnp(nsnp_val, default_nsnp)
        
        # 跳过无效结果
        if (is.na(beta) || is.na(se)) next
        
        or <- exp(beta)
        or_lci <- exp(beta - 1.96 * se)
        or_uci <- exp(beta + 1.96 * se)
        
        table_s1_data[[length(table_s1_data) + 1]] <- data.frame(
          category = category,
          exposure = exposure,
          outcome = outcome,
          method = method_name,
          n_snps = nsnp,
          beta = beta,
          se = se,
          pval = pval,
          or = or,
          or_lci = or_lci,
          or_uci = or_uci,
          or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
          extraction_strategy = strategy,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # IVW结果作为主结果
    ivw_result <- result$mr_results[result$mr_results$method == "Inverse variance weighted", ]
    if (nrow(ivw_result) == 0) {
      ivw_result <- result$mr_results[1, ]
    }
    
    if (nrow(ivw_result) > 0) {
      ivw_beta <- safe_extract(ivw_result$b, 1, NA)
      ivw_se <- safe_extract(ivw_result$se, 1, NA)
      ivw_pval <- safe_extract(ivw_result$pval, 1, NA)
      ivw_nsnp_val <- safe_extract(ivw_result$nsnp, 1, NA)
      default_nsnp <- safe_extract(result$n_snps_harmonized, 1, NA)
      ivw_nsnp <- safe_get_nsnp(ivw_nsnp_val, default_nsnp)
      
      # 检查是否有效值
      if (is.na(ivw_beta) || is.na(ivw_se)) {
        next  # 跳过无效结果
      }
      
      ivw_or <- exp(ivw_beta)
      ivw_or_lci <- exp(ivw_beta - 1.96 * ivw_se)
      ivw_or_uci <- exp(ivw_beta + 1.96 * ivw_se)
      
      # 提取其他方法
      egger_result <- result$mr_results[result$mr_results$method == "MR Egger", ]
      wm_result <- result$mr_results[result$mr_results$method == "Weighted median", ]
      wmode_result <- result$mr_results[result$mr_results$method == "Weighted mode", ]
      smode_result <- result$mr_results[result$mr_results$method == "Simple mode", ]
      
      # 异质性检验
      het_p <- NA
      if (!is.null(result$heterogeneity) && nrow(result$heterogeneity) > 0) {
        ivw_het <- result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", ]
        if (nrow(ivw_het) > 0) {
          het_p <- safe_extract(ivw_het$Q_pval, 1, NA)
        }
      }
      
      # 多效性检验
      pleo_p <- NA
      pleo_intercept <- NA
      if (!is.null(result$pleiotropy) && nrow(result$pleiotropy) > 0) {
        pleo_p <- safe_extract(result$pleiotropy$pval, 1, NA)
        pleo_intercept <- safe_extract(result$pleiotropy$egger_intercept, 1, NA)
      }
      
      # MR-PRESSO结果
      presso_global_p <- NA
      presso_outlier_corrected_or <- NA
      presso_outlier_corrected_or_lci <- NA
      presso_outlier_corrected_or_uci <- NA
      presso_n_outliers <- 0
      
      if (!is.null(result$presso)) {
        tryCatch({
          presso_global_p <- result$presso$`MR-PRESSO results`$`Global Test`$Pvalue
          if (!is.null(result$presso$`MR-PRESSO results`$`Outlier Test`)) {
            presso_n_outliers <- nrow(result$presso$`MR-PRESSO results`$`Outlier Test`)
          }
          if (!is.null(result$presso$`Main MR results`)) {
            presso_beta <- result$presso$`Main MR results`[1, "Causal Estimate"]
            presso_se <- result$presso$`Main MR results`[1, "Sd"]
            presso_outlier_corrected_or <- exp(presso_beta)
            presso_outlier_corrected_or_lci <- exp(presso_beta - 1.96 * presso_se)
            presso_outlier_corrected_or_uci <- exp(presso_beta + 1.96 * presso_se)
          }
        }, error = function(e) {
          # 静默处理
        })
      }
      
      # Leave-one-out稳定性
      loo_stable <- "NA"
      loo_all_significant <- FALSE
      if (!is.null(result$loo) && nrow(result$loo) > 1) {
        loo_pvals <- result$loo$p[!is.na(result$loo$p)]
        if (length(loo_pvals) > 0) {
          all_significant <- all(loo_pvals < 0.05, na.rm = TRUE)
          all_non_significant <- all(loo_pvals >= 0.05, na.rm = TRUE)
          
          if (all_significant && !is.na(ivw_pval) && ivw_pval < 0.05) {
            loo_stable <- "Stable (all significant)"
            loo_all_significant <- TRUE
          } else if (all_non_significant && (is.na(ivw_pval) || ivw_pval >= 0.05)) {
            loo_stable <- "Stable (all non-significant)"
          } else {
            loo_stable <- "Unstable"
          }
        }
      }
      
      # 工具变量强度（使用safe_ifelse）
      mean_f <- NA
      r_squared <- NA
      if (!is.null(result$iv_strength)) {
        mean_f <- safe_extract(result$iv_strength$mean_f, 1, NA)
        r_squared <- safe_extract(result$iv_strength$r_squared, 1, NA)
      }
      
      # 安全提取Egger结果
      egger_beta_val <- if (nrow(egger_result) > 0) safe_extract(egger_result$b, 1, NA) else NA
      egger_se_val <- if (nrow(egger_result) > 0) safe_extract(egger_result$se, 1, NA) else NA
      egger_pval_val <- if (nrow(egger_result) > 0) safe_extract(egger_result$pval, 1, NA) else NA
      egger_or_95ci_val <- if (nrow(egger_result) > 0 && !is.na(egger_beta_val) && !is.na(egger_se_val)) {
        egger_or <- exp(egger_beta_val)
        sprintf("%.3f (%.3f-%.3f)", 
               egger_or, 
               exp(egger_beta_val - 1.96 * egger_se_val),
               exp(egger_beta_val + 1.96 * egger_se_val))
      } else {
        NA
      }
      
      # 安全提取Weighted Median结果
      wm_beta_val <- if (nrow(wm_result) > 0) safe_extract(wm_result$b, 1, NA) else NA
      wm_se_val <- if (nrow(wm_result) > 0) safe_extract(wm_result$se, 1, NA) else NA
      wm_pval_val <- if (nrow(wm_result) > 0) safe_extract(wm_result$pval, 1, NA) else NA
      wm_or_95ci_val <- if (nrow(wm_result) > 0 && !is.na(wm_beta_val) && !is.na(wm_se_val)) {
        wm_or <- exp(wm_beta_val)
        sprintf("%.3f (%.3f-%.3f)", 
               wm_or, 
               exp(wm_beta_val - 1.96 * wm_se_val),
               exp(wm_beta_val + 1.96 * wm_se_val))
      } else {
        NA
      }
      
      # Table 2增强版
      table2_enhanced_data[[length(table2_enhanced_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        method = "IVW",
        n_snps = ivw_nsnp,
        or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
        pval = ivw_pval,
        egger_or_95ci = egger_or_95ci_val,
        egger_pval = egger_pval_val,
        wm_or_95ci = wm_or_95ci_val,
        wm_pval = wm_pval_val,
        pleiotropy_p = pleo_p,
        heterogeneity_p = het_p,
        stringsAsFactors = FALSE
      )
      
      # Table 5: 敏感性分析汇总（使用safe_ifelse）
      presso_or_95ci_str <- NA
      if (!is.na(presso_outlier_corrected_or) && !is.na(presso_outlier_corrected_or_lci) && !is.na(presso_outlier_corrected_or_uci)) {
        presso_or_95ci_str <- sprintf("%.3f (%.3f-%.3f)", 
                                     presso_outlier_corrected_or,
                                     presso_outlier_corrected_or_lci,
                                     presso_outlier_corrected_or_uci)
      }
      
      table5_data[[length(table5_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        n_snps = ivw_nsnp,
        ivw_or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
        ivw_pval = ivw_pval,
        egger_intercept_p = pleo_p,
        heterogeneity_q_p = het_p,
        presso_global_p = presso_global_p,
        presso_n_outliers = presso_n_outliers,
        presso_outlier_corrected_or_95ci = presso_or_95ci_str,
        loo_stability = loo_stable,
        stringsAsFactors = FALSE
      )
      
      # 安全提取异质性Q值
      het_q <- NA
      if (!is.null(result$heterogeneity) && nrow(result$heterogeneity) > 0) {
        ivw_het <- result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", ]
        if (nrow(ivw_het) > 0) {
          het_q <- safe_extract(ivw_het$Q, 1, NA)
        }
      }
      
      # Table S2: 敏感性分析详情
      table_s2_data[[length(table_s2_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        method = "IVW",
        beta = ivw_beta,
        se = ivw_se,
        pval = ivw_pval,
        or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
        egger_beta = egger_beta_val,
        egger_se = egger_se_val,
        egger_pval = egger_pval_val,
        wm_beta = wm_beta_val,
        wm_se = wm_se_val,
        wm_pval = wm_pval_val,
        pleiotropy_intercept = pleo_intercept,
        pleiotropy_p = pleo_p,
        heterogeneity_q = het_q,
        heterogeneity_q_p = het_p,
        presso_global_p = presso_global_p,
        presso_outlier_corrected_or_95ci = presso_or_95ci_str,
        stringsAsFactors = FALSE
      )
      
      # Table S6: 工具变量强度
      table_s6_data[[length(table_s6_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        n_snps = ivw_nsnp,
        mean_f_statistic = mean_f,
        r_squared = r_squared,
        extraction_strategy = strategy,
        stringsAsFactors = FALSE
      )
      
      # Table S3: Leave-one-out详情
      if (!is.null(result$loo) && nrow(result$loo) > 0) {
        for (j in seq_len(nrow(result$loo))) {
          loo_row <- result$loo[j, ]
          loo_beta <- safe_extract(loo_row$b, 1, NA)
          loo_se <- safe_extract(loo_row$se, 1, NA)
          
          if (is.na(loo_beta) || is.na(loo_se)) next
          
          loo_or <- exp(loo_beta)
          loo_or_lci <- exp(loo_beta - 1.96 * loo_se)
          loo_or_uci <- exp(loo_beta + 1.96 * loo_se)
          
          loo_snp_name <- safe_extract(loo_row$SNP, 1, paste("SNP", j))
          if (is.null(loo_snp_name) || is.na(loo_snp_name) || length(loo_snp_name) == 0) {
            loo_snp_name <- paste("SNP", j)
          }
          
          loo_nsnp <- safe_extract(loo_row$nsnp, 1, NA)
          loo_pval <- safe_extract(loo_row$p, 1, NA)
          
          loo_significant <- "No"
          if (!is.na(loo_pval) && loo_pval < 0.05) {
            loo_significant <- "Yes"
          }
          
          table_s3_data[[length(table_s3_data) + 1]] <- data.frame(
            category = category,
            exposure = exposure,
            outcome = outcome,
            snp_removed = as.character(loo_snp_name[1]),
            n_snps_remaining = loo_nsnp,
            or_95ci = sprintf("%.3f (%.3f-%.3f)", loo_or, loo_or_lci, loo_or_uci),
            pval = loo_pval,
            significant = loo_significant,
            stringsAsFactors = FALSE
          )
        }
      }
      
      processed_count <- processed_count + 1
      if (processed_count %% 10 == 0) {
        cat(sprintf("  已处理: %d/%d\n", processed_count, length(successful_results)))
      }
    }
  }, error = function(e) {
    cat(sprintf("  警告：处理 %s 时出错: %s\n", key, e$message))
  })
}

cat(sprintf("✓ 处理完成: %d 个结果\n\n", processed_count))

# ============================================================================
# 步骤5：合并和保存表格
# ============================================================================

cat("【步骤5】合并和保存表格\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# 合并所有表格
table_s1_full <- if (length(table_s1_data) > 0) {
  do.call(rbind, table_s1_data)
} else {
  data.frame()
}

table2_enhanced_full <- if (length(table2_enhanced_data) > 0) {
  do.call(rbind, table2_enhanced_data)
} else {
  data.frame()
}

table5_full <- if (length(table5_data) > 0) {
  do.call(rbind, table5_data)
} else {
  data.frame()
}

table_s2_full <- if (length(table_s2_data) > 0) {
  do.call(rbind, table_s2_data)
} else {
  data.frame()
}

table_s3_full <- if (length(table_s3_data) > 0) {
  do.call(rbind, table_s3_data)
} else {
  data.frame()
}

table_s6_full <- if (length(table_s6_data) > 0) {
  do.call(rbind, table_s6_data)
} else {
  data.frame()
}

# 添加显著性标记并排序
if (nrow(table_s1_full) > 0) {
  table_s1_full$significance <- ifelse(table_s1_full$pval < 0.001, "***",
                                       ifelse(table_s1_full$pval < 0.01, "**",
                                              ifelse(table_s1_full$pval < 0.05, "*", "")))
  table_s1_full <- table_s1_full[order(table_s1_full$category, 
                                      table_s1_full$exposure, 
                                      table_s1_full$outcome,
                                      table_s1_full$method), ]
}

if (nrow(table2_enhanced_full) > 0) {
  table2_enhanced_full$ivw_significant <- ifelse(table2_enhanced_full$pval < 0.05, "Yes", "No")
  table2_enhanced_full <- table2_enhanced_full[order(table2_enhanced_full$category, 
                                                    table2_enhanced_full$exposure, 
                                                    table2_enhanced_full$outcome), ]
}

# 保存为CSV
cat("正在保存CSV文件...\n")
if (nrow(table_s1_full) > 0) {
  write.csv(table_s1_full, "results/tables/paper_tables/Table_S1_Reverse_MR_Full_Results.csv", row.names = FALSE)
  cat("  ✓ Table_S1_Reverse_MR_Full_Results.csv\n")
}
if (nrow(table2_enhanced_full) > 0) {
  write.csv(table2_enhanced_full, "results/tables/paper_tables/Table_2_Reverse_MR_Enhanced.csv", row.names = FALSE)
  cat("  ✓ Table_2_Reverse_MR_Enhanced.csv\n")
}
if (nrow(table5_full) > 0) {
  write.csv(table5_full, "results/tables/paper_tables/Table_5_Reverse_MR_Sensitivity.csv", row.names = FALSE)
  cat("  ✓ Table_5_Reverse_MR_Sensitivity.csv\n")
}
if (nrow(table_s2_full) > 0) {
  write.csv(table_s2_full, "results/tables/paper_tables/Table_S2_Reverse_MR_Sensitivity_Details.csv", row.names = FALSE)
  cat("  ✓ Table_S2_Reverse_MR_Sensitivity_Details.csv\n")
}
if (nrow(table_s3_full) > 0) {
  write.csv(table_s3_full, "results/tables/paper_tables/Table_S3_Reverse_MR_Leave_One_Out.csv", row.names = FALSE)
  cat("  ✓ Table_S3_Reverse_MR_Leave_One_Out.csv\n")
}
if (nrow(table_s6_full) > 0) {
  write.csv(table_s6_full, "results/tables/paper_tables/Table_S6_Reverse_MR_IV_Strength.csv", row.names = FALSE)
  cat("  ✓ Table_S6_Reverse_MR_IV_Strength.csv\n")
}

# 创建多工作表Excel文件
cat("\n正在创建Excel文件...\n")
wb <- createWorkbook()
if (nrow(table_s1_full) > 0) {
  addWorksheet(wb, "Table_S1_Full_MR")
  writeData(wb, "Table_S1_Full_MR", table_s1_full)
  cat("  ✓ 已添加工作表: Table_S1_Full_MR\n")
}
if (nrow(table2_enhanced_full) > 0) {
  addWorksheet(wb, "Table_2_Enhanced")
  writeData(wb, "Table_2_Enhanced", table2_enhanced_full)
  cat("  ✓ 已添加工作表: Table_2_Enhanced\n")
}
if (nrow(table5_full) > 0) {
  addWorksheet(wb, "Table_5_Sensitivity")
  writeData(wb, "Table_5_Sensitivity", table5_full)
  cat("  ✓ 已添加工作表: Table_5_Sensitivity\n")
}
if (nrow(table_s2_full) > 0) {
  addWorksheet(wb, "Table_S2_Details")
  writeData(wb, "Table_S2_Details", table_s2_full)
  cat("  ✓ 已添加工作表: Table_S2_Details\n")
}
if (nrow(table_s3_full) > 0) {
  addWorksheet(wb, "Table_S3_LOO")
  writeData(wb, "Table_S3_LOO", table_s3_full)
  cat("  ✓ 已添加工作表: Table_S3_LOO\n")
}
if (nrow(table_s6_full) > 0) {
  addWorksheet(wb, "Table_S6_IV_Strength")
  writeData(wb, "Table_S6_IV_Strength", table_s6_full)
  cat("  ✓ 已添加工作表: Table_S6_IV_Strength\n")
}

saveWorkbook(wb, "results/tables/paper_tables/step11_all_paper_tables.xlsx", overwrite = TRUE)
cat("  ✓ step11_all_paper_tables.xlsx\n\n")

cat("✓ 所有论文表格已保存到 results/tables/paper_tables/\n")
cat("  - step11_all_paper_tables.xlsx (多工作表)\n")
cat("  - Table_S1_Reverse_MR_Full_Results.csv\n")
cat("  - Table_2_Reverse_MR_Enhanced.csv\n")
cat("  - Table_5_Reverse_MR_Sensitivity.csv\n")
cat("  - Table_S2_Reverse_MR_Sensitivity_Details.csv\n")
cat("  - Table_S3_Reverse_MR_Leave_One_Out.csv\n")
cat("  - Table_S6_Reverse_MR_IV_Strength.csv\n\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤8完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

