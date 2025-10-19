# 执行孟德尔随机化分析
cat("开始孟德尔随机化分析...\n")

mr_results_list <- list()

for (exposure_name in names(exposure_data_list)) {
  for (outcome_name in names(outcome_data_list)) {
    
    cat("分析:", exposure_name, "->", outcome_name, "\n")
    
    exposure_data <- exposure_data_list[[exposure_name]]
    outcome_data <- outcome_data_list[[outcome_name]]
    
    if (is.null(exposure_data) || is.null(outcome_data)) {
      cat("  跳过: 数据缺失\n")
      next
    }
    
    # 协调数据
    harmonised_data <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_data,
      outcome_dat = outcome_data
    )
    
    if (nrow(harmonised_data) == 0) {
      cat("  跳过: 协调后无SNP\n")
      next
    }
    
    # 执行MR分析
    mr_results <- TwoSampleMR::mr(harmonised_data, method_list = ANALYSIS_PARAMS$mr_methods)
    
    # 计算OR和置信区间
    mr_results$or <- exp(mr_results$b)
    mr_results$or_lci95 <- exp(mr_results$b - 1.96 * mr_results$se)
    mr_results$or_uci95 <- exp(mr_results$b + 1.96 * mr_results$se)
    
    # 添加元数据
    mr_results$exposure <- exposure_name
    mr_results$outcome <- outcome_name
    mr_results$nsnps <- nrow(harmonised_data)
    
    # 存储结果
    result_id <- paste(exposure_name, outcome_name, sep = "_")
    mr_results_list[[result_id]] <- list(
      harmonised_data = harmonised_data,
      mr_results = mr_results
    )
    
    cat("  完成:", nrow(harmonised_data), "个SNP\n")
  }
}

# 保存MR结果
saveRDS(mr_results_list, "02_mr_analysis/mr_results_list.rds")

# 合并所有结果
all_mr_results <- do.call(rbind, lapply(mr_results_list, function(x) x$mr_results))
write.csv(all_mr_results, "02_mr_analysis/all_mr_results.csv", row.names = FALSE)

cat("孟德尔随机化分析完成!\n")