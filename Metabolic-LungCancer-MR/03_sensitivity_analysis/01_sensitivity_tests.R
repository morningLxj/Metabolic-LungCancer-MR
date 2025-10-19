# 敏感性分析
cat("执行敏感性分析...\n")

sensitivity_results <- list()

for (result_id in names(mr_results_list)) {
  cat("敏感性分析:", result_id, "\n")
  
  harmonised_data <- mr_results_list[[result_id]]$harmonised_data
  
  if (nrow(harmonised_data) < 3) {
    cat("  跳过: SNP数量不足\n")
    next
  }
  
  # 异质性检验
  heterogeneity <- TwoSampleMR::mr_heterogeneity(harmonised_data)
  
  # 多效性检验
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonised_data)
  
  # 留一法分析
  leaveoneout <- TwoSampleMR::mr_leaveoneout(harmonised_data)
  
  # MR-PRESSO分析
  tryCatch({
    presso <- MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", 
      BetaExposure = "beta.exposure", 
      SdOutcome = "se.outcome", 
      SdExposure = "se.exposure", 
      OUTLIERtest = TRUE, 
      DISTORTIONtest = TRUE, 
      data = harmonised_data,
      NbDistribution = 1000, 
      SignifThreshold = 0.05
    )
  }, error = function(e) {
    presso <- NULL
  })
  
  sensitivity_results[[result_id]] <- list(
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy,
    leaveoneout = leaveoneout,
    presso = presso
  )
}

# 保存敏感性分析结果
saveRDS(sensitivity_results, "03_sensitivity_analysis/sensitivity_results.rds")
cat("敏感性分析完成!\n")