# 环境设置脚本
cat("设置分析环境...\n")

# 创建目录结构
dirs_to_create <- c(
  "01_data_preparation", "02_mr_analysis", "03_sensitivity_analysis",
  "04_visualization", "05_results", "scripts", "config", "logs"
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("创建目录:", dir, "\n")
  }
}

# 安装所需包
source("scripts/utility_functions.R")
install_required_packages()

# 保存会话信息
sink("session_info.txt")
print(sessionInfo())
sink()

cat("环境设置完成!\n")
cat("现在可以运行: source('main_analysis.R')\n")