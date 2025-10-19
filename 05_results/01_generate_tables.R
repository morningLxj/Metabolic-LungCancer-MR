# 生成结果表格
cat("生成结果表格...\n")

library(openxlsx)

# 创建Excel工作簿
wb <- createWorkbook()

# 表1: 主要MR结果
addWorksheet(wb, "Table1_Main_MR_Results")
main_results_table <- all_mr_results[all_mr_results$method == "Inverse variance weighted", 
                                    c("exposure", "outcome", "nsnps", "pval", "or", "or_lci95", "or_uci95")]
writeData(wb, "Table1_Main_MR_Results", main_results_table)

# 表2: 工具变量质量
addWorksheet(wb, "Table2_IV_Quality")
writeData(wb, "Table2_IV_Quality", iv_strength_df)

# 表3: 敏感性分析结果
addWorksheet(wb, "Table3_Sensitivity_Analysis")
sensitivity_summary <- do.call(rbind, lapply(names(sensitivity_results), function(x) {
  data.frame(
    analysis = x,
    heterogeneity_pval = ifelse(!is.null(sensitivity_results[[x]]$heterogeneity), 
                               sensitivity_results[[x]]$heterogeneity$Q_pval[1], NA),
    pleiotropy_pval = ifelse(!is.null(sensitivity_results[[x]]$pleiotropy),
                            sensitivity_results[[x]]$pleiotropy$pval, NA)
  )
}))
writeData(wb, "Table3_Sensitivity_Analysis", sensitivity_summary)

# 保存Excel文件
saveWorkbook(wb, "05_results/mr_analysis_results.xlsx", overwrite = TRUE)

# 生成README文件
generate_readme_file()

cat("结果表格生成完成!\n")