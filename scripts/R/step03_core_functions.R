
############################################################################
# 步骤3：定义核心分析函数库（精简版）
# 专注于基础MR和中介分析逻辑，不包含数据提取
# 包含防御性检查，提高代码健壮性
############################################################################

cat("步骤3：定义核心分析函数库...\n")

# 加载必要的包
required_packages <- c("TwoSampleMR", "dplyr", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# 确保包已正确加载
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  stop("TwoSampleMR包未安装，请先安装该包")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("dplyr包未安装，请先安装该包")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2包未安装，请先安装该包")
}

# 3.1 基础MR分析函数（简化版，供后续步骤调用）
perform_basic_mr <- function(exposure_data, outcome_id, outcome_name) {
  tryCatch({
    # 检查输入数据
    if (is.null(exposure_data) || nrow(exposure_data) == 0) {
      warning("MR分析：暴露数据为空")
      return(NULL)
    }
    
    if (!"SNP" %in% names(exposure_data)) {
      warning("MR分析：暴露数据缺少SNP列")
      return(NULL)
    }
    
    outcome_data <- TwoSampleMR::extract_outcome_data(
      snps = exposure_data$SNP,
      outcomes = outcome_id
    )
    
    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      return(NULL)
    }
    
    harmonized <- TwoSampleMR::harmonise_data(exposure_data, outcome_data)
    
    if (is.null(harmonized) || nrow(harmonized) == 0) {
      return(NULL)
    }
    
    mr_results <- TwoSampleMR::mr(harmonized)
    
    # 检查mr_results是否有效
    if (is.null(mr_results) || nrow(mr_results) == 0) {
      return(NULL)
    }
    
    # 安全计算OR值（检查必需的列是否存在）
    if ("b" %in% names(mr_results) && "se" %in% names(mr_results)) {
      mr_results$OR <- exp(mr_results$b)
      mr_results$OR_lci95 <- exp(mr_results$b - 1.96 * mr_results$se)
      mr_results$OR_uci95 <- exp(mr_results$b + 1.96 * mr_results$se)
    } else {
      warning("MR分析结果缺少必需的列 (b, se)，无法计算OR值")
    }
    
    # 敏感性分析（可能失败，使用tryCatch包裹）
    heterogeneity <- tryCatch({
      TwoSampleMR::mr_heterogeneity(harmonized)
    }, error = function(e) {
      warning("异质性检验失败:", e$message)
      return(NULL)
    })
    
    pleiotropy <- tryCatch({
      TwoSampleMR::mr_pleiotropy_test(harmonized)
    }, error = function(e) {
      warning("多效性检验失败:", e$message)
      return(NULL)
    })
    
    return(list(
      harmonized_data = harmonized,
      mr_results = mr_results,
      heterogeneity = heterogeneity,
      pleiotropy = pleiotropy,
      n_snps = nrow(harmonized)
    ))
    
  }, error = function(e) {
    cat("MR分析错误:", e$message, "\n")
    return(NULL)
  })
}

# 3.2 中介效应计算函数
calculate_mediation_effect <- function(exp_to_med_result, 
                                      med_to_out_result, 
                                      exp_to_out_result = NULL) {
  
  # 检查输入数据的有效性
  if (is.null(exp_to_med_result) || is.null(exp_to_med_result$mr_results) ||
      is.null(med_to_out_result) || is.null(med_to_out_result$mr_results)) {
    return(NULL)
  }
  
  # 检查必需的列是否存在
  required_cols <- c("method", "b", "se")
  if (!all(required_cols %in% names(exp_to_med_result$mr_results)) ||
      !all(required_cols %in% names(med_to_out_result$mr_results))) {
    warning("中介效应计算：缺少必需的列 (method, b, se)")
    return(NULL)
  }
  
  # 提取IVW结果
  exp_to_med_ivw <- exp_to_med_result$mr_results %>% 
    dplyr::filter(method == "Inverse variance weighted")
  med_to_out_ivw <- med_to_out_result$mr_results %>% 
    dplyr::filter(method == "Inverse variance weighted")
  
  if (nrow(exp_to_med_ivw) == 0 || nrow(med_to_out_ivw) == 0) {
    return(NULL)
  }
  
  # 安全提取系数和标准误
  alpha <- if ("b" %in% names(exp_to_med_ivw)) exp_to_med_ivw$b[1] else NA_real_
  beta <- if ("b" %in% names(med_to_out_ivw)) med_to_out_ivw$b[1] else NA_real_
  se_alpha <- if ("se" %in% names(exp_to_med_ivw)) exp_to_med_ivw$se[1] else NA_real_
  se_beta <- if ("se" %in% names(med_to_out_ivw)) med_to_out_ivw$se[1] else NA_real_
  
  if (is.na(alpha) || is.na(beta) || is.na(se_alpha) || is.na(se_beta)) {
    return(NULL)
  }
  
  # 计算间接效应
  indirect_effect <- alpha * beta
  se_indirect <- sqrt((alpha^2 * se_beta^2) + (beta^2 * se_alpha^2))
  z_indirect <- indirect_effect / se_indirect
  p_indirect <- 2 * (1 - pnorm(abs(z_indirect)))
  
  # 计算总效应和直接效应
  total_effect <- NULL
  direct_effect <- NULL
  mediation_proportion <- NULL
  
  if (!is.null(exp_to_out_result) && !is.null(exp_to_out_result$mr_results)) {
    if (all(required_cols %in% names(exp_to_out_result$mr_results))) {
      exp_to_out_ivw <- exp_to_out_result$mr_results %>% 
        dplyr::filter(method == "Inverse variance weighted")
      
      if (nrow(exp_to_out_ivw) > 0 && "b" %in% names(exp_to_out_ivw)) {
        total_effect <- exp_to_out_ivw$b[1]
        direct_effect <- total_effect - indirect_effect
        
        if (!is.na(total_effect) && abs(total_effect) > 1e-10) {
          mediation_proportion <- indirect_effect / total_effect
        }
      }
    }
  }
  
  # 安全提取p值
  exp_to_med_pval <- if ("pval" %in% names(exp_to_med_ivw)) exp_to_med_ivw$pval[1] else NA_real_
  med_to_out_pval <- if ("pval" %in% names(med_to_out_ivw)) med_to_out_ivw$pval[1] else NA_real_
  
  return(list(
    indirect_effect = indirect_effect,
    indirect_se = se_indirect,
    indirect_pval = p_indirect,
    total_effect = total_effect,
    direct_effect = direct_effect,
    mediation_proportion = mediation_proportion,
    alpha = alpha,
    beta = beta,
    exp_to_med_pval = exp_to_med_pval,
    med_to_out_pval = med_to_out_pval
  ))
}

# 3.3 结果汇总函数
summarize_mr_result <- function(mr_result, exposure_name, outcome_name) {
  if (is.null(mr_result) || is.null(mr_result$mr_results)) {
    return(NULL)
  }
  
  # 检查必需的列是否存在
  if (!"method" %in% names(mr_result$mr_results)) {
    warning("结果汇总：缺少method列")
    return(NULL)
  }
  
  ivw_result <- mr_result$mr_results[mr_result$mr_results$method == "Inverse variance weighted", ]
  
  if (nrow(ivw_result) > 0) {
    # 安全提取列
    beta <- if ("b" %in% names(ivw_result)) ivw_result$b[1] else NA_real_
    se <- if ("se" %in% names(ivw_result)) ivw_result$se[1] else NA_real_
    pval <- if ("pval" %in% names(ivw_result)) ivw_result$pval[1] else NA_real_
    n_snps <- if (!is.null(mr_result$n_snps)) mr_result$n_snps else NA_integer_
    
    if (is.na(beta) || is.na(se)) {
      return(NULL)
    }
    
    return(data.frame(
      exposure = exposure_name,
      outcome = outcome_name,
      n_snps = n_snps,
      beta = beta,
      se = se,
      pval = pval,
      or = exp(beta),
      or_lci = exp(beta - 1.96 * se),
      or_uci = exp(beta + 1.96 * se),
      stringsAsFactors = FALSE
    ))
  }
  
  return(NULL)
}

# 3.4 可视化函数定义
create_forest_plot <- function(data, outcome_filter, title) {
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  # 检查必需的列是否存在
  required_cols <- c("outcome", "or", "exposure")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    warning(paste("创建森林图：缺少必需的列:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  plot_data <- data[data$outcome == outcome_filter & !is.na(data$or), ]
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  plot_data <- plot_data[order(plot_data$or), ]
  
  # 检查是否有category列，如果没有则创建一个默认值
  if (!"category" %in% names(plot_data)) {
    plot_data$category <- "Unknown"
  }
  
  # 检查是否有置信区间列
  has_ci <- "or_lci" %in% names(plot_data) && "or_uci" %in% names(plot_data)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = or, y = stats::reorder(exposure, or))) +
    ggplot2::geom_point(ggplot2::aes(color = category), size = 3) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8)
  
  # 如果有置信区间，添加误差条
  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = or_lci, xmax = or_uci, color = category),
                            height = 0.2, alpha = 0.7)
  }
  
  p <- p +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::labs(title = title, x = "Odds Ratio (95% CI)", y = "Exposure") +
    ggplot2::theme_bw()
  
  return(p)
}

# 3.5 质量检查函数
check_instrument_quality <- function(instruments) {
  if (is.null(instruments) || nrow(instruments) == 0) {
    return(list(valid = FALSE, message = "No instruments available"))
  }
  
  # 检查必需的列是否存在，尝试多个可能的列名
  beta_col <- NULL
  se_col <- NULL
  
  # 优先查找 beta.exposure 和 se.exposure
  if ("beta.exposure" %in% names(instruments)) {
    beta_col <- "beta.exposure"
  } else if ("beta" %in% names(instruments)) {
    beta_col <- "beta"
  }
  
  if ("se.exposure" %in% names(instruments)) {
    se_col <- "se.exposure"
  } else if ("se" %in% names(instruments)) {
    se_col <- "se"
  }
  
  if (is.null(beta_col) || is.null(se_col)) {
    return(list(
      valid = FALSE,
      n_snps = nrow(instruments),
      mean_f_stat = NA_real_,
      weak_instruments = NA_integer_,
      weak_proportion = NA_real_,
      message = "Missing required columns (beta.exposure/beta, se.exposure/se)"
    ))
  }
  
  # 计算F统计量
  beta_values <- instruments[[beta_col]]
  se_values <- instruments[[se_col]]
  
  # 安全计算F统计量（避免除零错误）
  valid_mask <- !is.na(beta_values) & !is.na(se_values) & se_values != 0
  if (sum(valid_mask) == 0) {
    return(list(
      valid = FALSE,
      n_snps = nrow(instruments),
      mean_f_stat = NA_real_,
      weak_instruments = NA_integer_,
      weak_proportion = NA_real_,
      message = "No valid data for F-statistic calculation"
    ))
  }
  
  instruments$f_stat <- NA_real_
  instruments$f_stat[valid_mask] <- (beta_values[valid_mask]^2) / (se_values[valid_mask]^2)
  mean_f <- mean(instruments$f_stat, na.rm = TRUE)
  
  # 检查弱工具变量
  weak_instruments <- sum(instruments$f_stat < 10, na.rm = TRUE)
  weak_proportion <- weak_instruments / sum(valid_mask)
  
  return(list(
    valid = nrow(instruments) >= 3 && !is.na(mean_f) && mean_f >= 10,
    n_snps = nrow(instruments),
    mean_f_stat = mean_f,
    weak_instruments = weak_instruments,
    weak_proportion = weak_proportion,
    message = ifelse(nrow(instruments) < 3, "Too few SNPs",
                    ifelse(is.na(mean_f) || mean_f < 10, "Weak instruments", "Good quality"))
  ))
}

# 保存核心函数
dir.create("data", showWarnings = FALSE, recursive = TRUE)
save(perform_basic_mr, 
     calculate_mediation_effect,
     summarize_mr_result,
     create_forest_plot,
     check_instrument_quality,
     file = "data/step03_core_functions.RData")

cat("\n步骤3完成！\n")
cat("- 已定义5个核心函数\n")
cat("- perform_basic_mr: 基础MR分析\n")
cat("- calculate_mediation_effect: 中介效应计算\n")
cat("- summarize_mr_result: 结果汇总\n")
cat("- create_forest_plot: 森林图绘制\n")
cat("- check_instrument_quality: 工具变量质量检查\n")
cat("- 函数已保存至: data/step03_core_functions.RData\n")
cat("\n下一步：运行 step04_extract_instruments.R\n")

