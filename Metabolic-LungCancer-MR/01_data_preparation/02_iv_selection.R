# 工具变量质量评估
cat("评估工具变量质量...\n")

# 计算每个暴露的工具变量强度
iv_strength_results <- list()

for (trait_name in names(exposure_data_list)) {
  exposure_data <- exposure_data_list[[trait_name]]
  strength_summary <- assess_instrument_strength(exposure_data)
  iv_strength_results[[trait_name]] <- strength_summary
}

# 合并结果
iv_strength_df <- do.call(rbind, iv_strength_results)
iv_strength_df$trait <- rownames(iv_strength_df)

# 输出工具变量质量报告
cat("\n工具变量质量报告:\n")
print(iv_strength_df)

# 保存质量评估结果
write.csv(iv_strength_df, "01_data_preparation/instrument_strength_assessment.csv", row.names = FALSE)

# 检查是否有工具变量过弱
weak_instruments <- iv_strength_df[iv_strength_df$mean_fstat < ANALYSIS_PARAMS$f_statistic_threshold, ]
if (nrow(weak_instruments) > 0) {
  cat("警告: 以下性状的工具变量可能偏弱 (F < 10):\n")
  print(weak_instruments)
} else {
  cat("✓ 所有工具变量强度均满足要求 (F > 10)\n")
}