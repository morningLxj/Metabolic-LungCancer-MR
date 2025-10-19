# 数据加载和预处理函数

#' 从OpenGWAS加载GWAS数据
load_gwas_data <- function(gwas_id, trait_name, type = "exposure") {
  cat("加载数据:", trait_name, "(", gwas_id, ")\n")
  
  tryCatch({
    if (type == "exposure") {
      # 对于暴露，我们需要提取工具变量
      dat <- TwoSampleMR::extract_instruments(
        outcomes = gwas_id,
        p1 = ANALYSIS_PARAMS$iv_selection$pval_threshold,
        clump = TRUE,
        r2 = ANALYSIS_PARAMS$iv_selection$clump_r2,
        kb = ANALYSIS_PARAMS$iv_selection$clump_kb
      )
    } else {
      # 对于结局，我们加载完整数据
      dat <- TwoSampleMR::extract_outcome_data(
        snps = NULL,  # 将在协调步骤中指定
        outcomes = gwas_id
      )
    }
    
    if (is.null(dat) || nrow(dat) == 0) {
      warning("未找到数据: ", trait_name)
      return(NULL)
    }
    
    # 添加性状名称
    dat$trait <- trait_name
    dat$gwas_id <- gwas_id
    
    return(dat)
    
  }, error = function(e) {
    warning("加载数据失败: ", trait_name, " - ", e$message)
    return(NULL)
  })
}

#' 计算F统计量
calculate_f_statistic <- function(beta, se, samplesize) {
  # F = (beta/se)^2
  f_stat <- (beta / se)^2
  return(f_stat)
}

#' 评估工具变量强度
assess_instrument_strength <- function(exposure_data) {
  if (is.null(exposure_data) || nrow(exposure_data) == 0) {
    return(data.frame(
      trait = NA,
      n_snps = 0,
      mean_fstat = NA,
      min_fstat = NA,
      max_fstat = NA
    ))
  }
  
  # 计算每个SNP的F统计量
  exposure_data$f_stat <- calculate_f_statistic(
    exposure_data$beta.exposure,
    exposure_data$se.exposure,
    exposure_data$samplesize.exposure
  )
  
  summary_stats <- data.frame(
    trait = unique(exposure_data$trait),
    n_snps = nrow(exposure_data),
    mean_fstat = mean(exposure_data$f_stat, na.rm = TRUE),
    min_fstat = min(exposure_data$f_stat, na.rm = TRUE),
    max_fstat = max(exposure_data$f_stat, na.rm = TRUE)
  )
  
  return(summary_stats)
}