# 工具函数

#' 设置日志
setup_logging <- function() {
  log_file <- paste0("logs/analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  if (!dir.exists("logs")) dir.create("logs")
  sink(file = log_file, split = TRUE)
}

#' 生成README文件
generate_readme_file <- function() {
  readme_content <- paste0(
    "# Mendelian Randomization Analysis: Metabolic Traits and Lung Cancer Subtypes\n\n",
    "## 分析概述\n",
    "本项目包含基于孟德尔随机化方法分析17种代谢性状与3种肺癌亚型因果关系的完整代码。\n\n",
    "## 数据来源\n",
    "所有GWAS数据来源于IEU OpenGWAS数据库和GWAS Catalog，具体ID见config/analysis_parameters.R\n\n",
    "## 分析步骤\n",
    "1. 数据加载和工具变量筛选\n", 
    "2. 孟德尔随机化主分析 (IVW方法)\n",
    "3. 敏感性分析 (MR-Egger, 加权中位数, MR-PRESSO)\n",
    "4. 异质性和多效性检验\n",
    "5. 结果可视化和表格生成\n\n",
    "## 运行说明\n",
    "运行 `main_analysis.R` 执行完整分析流程\n\n",
    "生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"
  )
  
  writeLines(readme_content, "README.md")
}

#' 安装所需包
install_required_packages <- function() {
  required_packages <- c(
    "TwoSampleMR", "MRPRESSO", "ggplot2", "dplyr", 
    "readr", "openxlsx", "forestplot", "knitr"
  )
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}