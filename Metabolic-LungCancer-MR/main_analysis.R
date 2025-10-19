#!/usr/bin/env Rscript
# 主分析脚本 - 孟德尔随机化分析流程
# 基于项目数据集表格生成

cat("=====================================\n")
cat("孟德尔随机化分析：代谢性状与肺癌亚型\n")
cat("基于项目数据集表格\n") 
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=====================================\n\n")

# 加载配置和函数
source("config/analysis_parameters.R")
source("scripts/utility_functions.R")
source("scripts/data_loading_functions.R")

# 设置日志
setup_logging()

tryCatch({
  # 步骤1: 数据准备
  cat("步骤1: 数据准备...\n")
  source("01_data_preparation/01_load_gwas_data.R")
  
  # 步骤2: 工具变量筛选
  cat("步骤2: 工具变量筛选...\n")
  source("01_data_preparation/02_iv_selection.R")
  
  # 步骤3: 主MR分析
  cat("步骤3: 孟德尔随机化分析...\n")
  source("02_mr_analysis/01_primary_mr.R")
  
  # 步骤4: 敏感性分析
  cat("步骤4: 敏感性分析...\n")
  source("03_sensitivity_analysis/01_sensitivity_tests.R")
  
  # 步骤5: 可视化
  cat("步骤5: 生成图表...\n")
  source("04_visualization/01_generate_figures.R")
  
  # 步骤6: 结果汇总
  cat("步骤6: 生成结果表格...\n")
  source("05_results/01_generate_tables.R")
  
  cat("\n✅ 分析完成！\n")
  cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
}, error = function(e) {
  cat("❌ 分析出错:", e$message, "\n")
  save.image("debug_environment.RData")
  stop("分析失败")
})