#!/usr/bin/env Rscript
# ============================================================================
# 甲基化重启分析 - 优化启动脚本
# 专门为已验证的39个候选探针设计的执行入口
# ============================================================================

cat("🚀 甲基化重启分析 - 启动脚本\n")
cat("=======================================\n\n")

# 设置工作目录
setwd("D:\\GWAS")

# 基础环境检查
cat("🔍 环境检查...\n")

# 检查关键文件
key_files <- c(
  "methylation_analysis/tcga_valid_probes.RData",
  "scripts/master_methylation_analysis.R"
)

all_files_exist <- TRUE
for(file in key_files) {
  if(file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  } else {
    cat(sprintf("✗ %s (缺失)\n", file))
    all_files_exist <- FALSE
  }
}

if(!all_files_exist) {
  cat("\n❌ 关键文件缺失，请检查目录结构\n")
  quit(status = 1)
}

# 加载候选探针数据
cat("\n📊 加载候选探针数据...\n")
tryCatch({
  load("methylation_analysis/tcga_valid_probes.RData")
  cat("✓ 候选探针数据加载成功\n")
  cat(sprintf("  - MFAP2探针: %d个 (总: %d, LUAD: %d, LUSC: %d)\n", 
              length(valid_mfap2), length(valid_mfap2), length(valid_mfap2_luad), length(valid_mfap2_lusc)))
  cat(sprintf("  - CDK11A探针: %d个 (总: %d, LUAD: %d, LUSC: %d)\n", 
              length(valid_cdk11a), length(valid_cdk11a), length(valid_cdk11a_luad), length(valid_cdk11a_lusc)))
}, error = function(e) {
  cat("❌ 加载探针数据失败:", e$message, "\n")
  quit(status = 1)
})

# 检查TCGA数据文件...
cat("\n📁 检查TCGA数据文件...\n")

tcga_dir <- "PDC/TCGA"
tcga_files <- c(
  paste0(tcga_dir, "/TCGA-LUAD.methylation450.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUSC.methylation450.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUAD.star_tpm.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUSC.star_tpm.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUAD.clinical.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUSC.clinical.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUAD.survival.tsv.gz"),
  paste0(tcga_dir, "/TCGA-LUSC.survival.tsv.gz")
)

missing_files <- tcga_files[!file.exists(tcga_files)]
if(length(missing_files) > 0) {
  cat("❌ 缺少TCGA数据文件:\n")
  for(file in missing_files) {
    cat("  -", file, "\n")
  }
  return(NULL)
} else {
  cat("✓ 所有TCGA数据文件检查完成\n")
}

# 创建输出目录
cat("\n📁 创建输出目录...\n")
output_dirs <- c(
  "methylation_analysis/dmp_results",
  "methylation_analysis/correlation_results", 
  "methylation_analysis/survival_results",
  "methylation_analysis/colocalization"
)

for(dir in output_dirs) {
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(sprintf("✓ 创建目录: %s\n", dir))
  }
}

# 执行分析脚本
cat("\n🎯 开始执行分析步骤...\n")

# 步骤1: DMP分析
cat("\n📊 步骤1: 差异甲基化探针 (DMP) 分析\n")
cat("-------------------------------------------\n")
if(file.exists("scripts/01_methylation_dmp_analysis.R")) {
  tryCatch({
    source("scripts/01_methylation_dmp_analysis.R")
    cat("✓ DMP分析脚本加载成功\n")
    cat("⚠ 请手动执行 run_complete_dmp_analysis() 函数\n")
  }, error = function(e) {
    cat("❌ DMP分析脚本加载失败:", e$message, "\n")
  })
} else {
  cat("❌ DMP分析脚本不存在\n")
}

# 步骤2: 相关性分析  
cat("\n🔗 步骤2: 甲基化-表达相关性分析\n")
cat("-------------------------------------------\n")
if(file.exists("scripts/02_methylation_expression_correlation.R")) {
  tryCatch({
    source("scripts/02_methylation_expression_correlation.R")
    cat("✓ 相关性分析脚本加载成功\n")
    cat("⚠ 请手动执行 run_correlation_analysis() 函数\n")
  }, error = function(e) {
    cat("❌ 相关性分析脚本加载失败:", e$message, "\n")
  })
} else {
  cat("❌ 相关性分析脚本不存在\n")
}

# 步骤3: 生存分析
cat("\n🏥 步骤3: 甲基化生存分析\n")
cat("-------------------------------------------\n")
if(file.exists("scripts/03_methylation_survival_analysis.R")) {
  tryCatch({
    source("scripts/03_methylation_survival_analysis.R")
    cat("✓ 生存分析脚本加载成功\n")
    cat("⚠ 请手动执行 run_survival_analysis() 函数\n")
  }, error = function(e) {
    cat("❌ 生存分析脚本加载失败:", e$message, "\n")
  })
} else {
  cat("❌ 生存分析脚本不存在\n")
}

# 步骤4: 共定位分析
cat("\n🎯 步骤4: mQTL与GWAS SNP共定位分析\n")
cat("-------------------------------------------\n")
if(file.exists("scripts/04_methylation_gwas_colocalization.R")) {
  tryCatch({
    source("scripts/04_methylation_gwas_colocalization.R")
    cat("✓ 共定位分析脚本加载成功\n")
    cat("⚠ 请手动执行 run_colocalization_analysis() 函数\n")
  }, error = function(e) {
    cat("❌ 共定位分析脚本加载失败:", e$message, "\n")
  })
} else {
  cat("❌ 共定位分析脚本不存在\n")
}

cat("\n🎉 甲基化重启分析启动完成！\n")
cat("==============================\n")
cat("📋 使用说明:\n")
cat("- 所有分析脚本已加载到当前环境\n")
cat("- 可以使用 ls() 查看可用函数\n")
cat("- 按照以下顺序执行分析:\n")
cat("  1. run_complete_dmp_analysis()\n")
cat("  2. run_correlation_analysis()\n") 
cat("  3. run_survival_analysis()\n")
cat("  4. run_colocalization_analysis()\n")
cat("\n✅ 环境准备完成，等待手动执行具体分析步骤\n")