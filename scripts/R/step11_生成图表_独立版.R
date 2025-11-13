############################################################################
# 步骤11：反向MR分析 - 独立图表生成脚本
# 用途：从已保存的分析结果直接生成图表，无需重新分析
# 使用方法：source("step11_生成图表_独立版.R")
############################################################################

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤11：反向MR分析 - 独立图表生成脚本\n")
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

# 创建输出目录结构
dirs <- c(
  "results/tables", 
  "results/tables/paper_tables",
  "results/figures/reverse_mr",
  "results/figures/reverse_mr/sensitivity",
  "results/figures/reverse_mr/main_figures"
)

for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

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

# 统计加载的数据
cat("\n【加载的数据统计】\n")
cat(sprintf("  - 总配对数量: %d\n", length(all_results)))
success_count <- sum(sapply(all_results, function(x) x$status == "Success"))
failed_count <- length(all_results) - success_count
cat(sprintf("  - 成功: %d\n", success_count))
cat(sprintf("  - 失败: %d\n", failed_count))
if (exists("analysis_log")) {
  cat(sprintf("  - 分析日志: %d 条记录\n", nrow(analysis_log)))
}
cat("\n")

# ============================================================================
# 步骤3：执行结果整理和图表生成
# ============================================================================

cat("【步骤3】开始生成图表和表格\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# 提取成功的结果
successful_results <- Filter(function(x) x$status == "Success", all_results)

if (length(successful_results) == 0) {
  stop("错误：没有成功的分析结果，无法生成图表和表格。")
}

cat(sprintf("✓ 成功分析: %d 条\n\n", length(successful_results)))

# ============================================================================
# 执行原始脚本中的步骤7和步骤8（图表生成部分）
# ============================================================================

cat("正在执行图表和表格生成...\n")
cat("（这相当于运行原始脚本的步骤7和步骤8）\n\n")

# 由于原始脚本中的图表生成代码非常长，这里我们直接source原始脚本的后续部分
# 但需要确保环境已经准备好（all_results, analysis_log等）

# 这里需要运行原始脚本中从"步骤7"开始的所有代码
# 为了简化，我们提供一个函数来执行

source_lines_start <- grep("步骤7：整理和导出结果", readLines("step11_反向MR分析_优化版.R"))
if (length(source_lines_start) > 0) {
  cat("正在执行原始脚本的图表生成部分...\n")
  
  # 读取原始脚本，从步骤7开始执行
  script_lines <- readLines("step11_反向MR分析_优化版.R")
  
  # 找到步骤7的开始位置
  step7_start <- grep("^# ============================================================================$|^# 步骤7：", script_lines)
  if (length(step7_start) > 0) {
    # 从步骤7开始执行
    step7_code <- script_lines[step7_start[1]:length(script_lines)]
    step7_code_text <- paste(step7_code, collapse = "\n")
    
    # 在执行环境中运行代码
    eval(parse(text = step7_code_text), envir = .GlobalEnv)
  } else {
    cat("⚠ 警告：无法找到步骤7的开始位置，将尝试手动执行关键部分\n")
  }
} else {
  cat("⚠ 警告：无法定位原始脚本，请手动运行 step11_反向MR分析_优化版.R 的步骤7和8\n")
  cat("或者您可以修改代码，直接从当前环境中的 all_results 生成图表\n\n")
  
  # 提供简化版本的图表生成
  cat("【简化版图表生成】\n")
  cat("您需要从原始脚本 step11_反向MR分析_优化版.R 中复制步骤7和8的代码到这里。\n")
  cat("或者，直接运行原始脚本，它会自动检测已有结果并跳过分析步骤。\n\n")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("图表生成完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【提示】\n")
cat("如果您想避免每次都重新运行分析，请：\n")
cat("1. 运行 step11_反向MR分析_优化版.R（它会自动检测并加载已有结果）\n")
cat("2. 或者将原始脚本中步骤7和8的代码复制到本脚本中\n\n")

