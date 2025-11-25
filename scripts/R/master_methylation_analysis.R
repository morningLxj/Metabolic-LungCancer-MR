#!/usr/bin/env Rscript
# ============================================================================
# 甲基化重启分析 - 主控制脚本
# 基于已验证的39个候选探针（MFAP2: 18个, CDK11A: 21个）
# 执行四步完整分析流程
# ============================================================================

# 加载必要的包
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(survival)
  library(survminer)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(DT)
})

cat("🚀 甲基化重启分析 - 主控制脚本\n")
cat("=====================================\n\n")

# ============================================================================
# 1. 环境设置与数据加载
# ============================================================================

setup_analysis_environment <- function() {
  cat("🔧 设置分析环境...\n")
  
  # 创建输出目录
  dirs_to_create <- c(
    "methylation_analysis/dmp_results",
    "methylation_analysis/correlation_results", 
    "methylation_analysis/survival_results",
    "methylation_analysis/colocalization"
  )
  
  for(dir in dirs_to_create) {
    if(!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("  ✓ 创建目录:", dir, "\n")
    }
  }
  
  # 加载验证的候选探针
  load("methylation_analysis/tcga_probe_summary.RData")
  load("methylation_analysis/comprehensive_debug.RData")
  
  cat("📊 加载候选探针数据:\n")
  cat(sprintf("  - MFAP2: %d个探针\n", tcga_summary$mfap2_count))
  cat(sprintf("  - CDK11A: %d个探针\n", tcga_summary$cdk11a_count))
  cat(sprintf("  - 总计: %d个探针\n\n", tcga_summary$total_candidates))
  
  return(list(
    mfap2_probes = mfap2_candidates,
    cdk11a_probes = cdk11a_candidates,
    target_genes = c("MFAP2", "CDK11A")
  ))
}

# ============================================================================
# 2. TCGA数据检查与准备
# ============================================================================

check_tcga_data_availability <- function() {
  cat("📁 检查TCGA数据可用性...\n")
  
  # 检查TCGA文件
  tcga_files <- c(
    "PDC/TCGA/TCGA-LUAD.methylation450.tsv.gz",
    "PDC/TCGA/TCGA-LUSC.methylation450.tsv.gz",
    "PDC/TCGA/TCGA-LUAD.expression.tsv.gz",
    "PDC/TCGA/TCGA-LUSC.expression.tsv.gz"
  )
  
  available_files <- sapply(tcga_files, file.exists)
  cat("可用数据文件:\n")
  for(i in seq_along(tcga_files)) {
    status <- ifelse(available_files[i], "✓", "✗")
    cat(sprintf("  %s %s\n", status, tcga_files[i]))
  }
  
  return(available_files)
}

# ============================================================================
# 3. 分析流程控制
# ============================================================================

run_analysis_pipeline <- function(target_genes, available_data) {
  cat("🎯 开始执行四步分析流程...\n\n")
  
  results <- list()
  
  # 步骤1: DMP分析
  if(file.exists("scripts/01_methylation_dmp_analysis.R")) {
    cat("📊 步骤1: DMP (差异甲基化探针) 分析\n")
    cat("-------------------------------------------\n")
    source("scripts/01_methylation_dmp_analysis.R")
    
    dmp_results <- run_complete_dmp_analysis(target_genes, available_data)
    results$dmp <- dmp_results
    
    cat("✓ DMP分析完成\n\n")
  }
  
  # 步骤2: 甲基化-表达相关性分析
  if(file.exists("scripts/02_methylation_expression_correlation.R")) {
    cat("🔗 步骤2: 甲基化-表达相关性分析\n")
    cat("-------------------------------------------\n")
    source("scripts/02_methylation_expression_correlation.R")
    
    correlation_results <- run_correlation_analysis(target_genes, available_data)
    results$correlation <- correlation_results
    
    cat("✓ 相关性分析完成\n\n")
  }
  
  # 步骤3: 甲基化生存分析
  if(file.exists("scripts/03_methylation_survival_analysis.R")) {
    cat("🏥 步骤3: 甲基化生存分析\n")
    cat("-------------------------------------------\n")
    source("scripts/03_methylation_survival_analysis.R")
    
    survival_results <- run_survival_analysis(target_genes, available_data)
    results$survival <- survival_results
    
    cat("✓ 生存分析完成\n\n")
  }
  
  # 步骤4: mQTL与GWAS SNP共定位分析
  if(file.exists("scripts/04_methylation_gwas_colocalization.R")) {
    cat("🎯 步骤4: mQTL与GWAS SNP共定位分析\n")
    cat("-------------------------------------------\n")
    source("scripts/04_methylation_gwas_colocalization.R")
    
    colocalization_results <- run_colocalization_analysis(target_genes, available_data)
    results$colocalization <- colocalization_results
    
    cat("✓ 共定位分析完成\n\n")
  }
  
  return(results)
}

# ============================================================================
# 4. 结果汇总与报告生成
# ============================================================================

generate_comprehensive_report <- function(analysis_results) {
  cat("📝 生成综合分析报告...\n")
  
  report_content <- paste0(
    "# 甲基化重启分析 - 综合报告\n",
    "生成时间: ", Sys.time(), "\n\n",
    
    "## 📊 分析概览\n",
    "- 目标基因: MFAP2, CDK11A\n",
    "- 候选探针: 39个 (MFAP2: 18个, CDK11A: 21个)\n",
    "- 分析类型: DMP → 相关性 → 生存 → 共定位\n\n",
    
    "## 🔍 主要发现\n"
  )
  
  # 添加各步骤结果摘要
  if(!is.null(analysis_results$dmp)) {
    report_content <- paste0(report_content, 
      "### DMP分析结果\n",
      "- 显著差异甲基化探针数量: ", analysis_results$dmp$significant_probes, "\n",
      "- 最显著探针: ", analysis_results$dmp$top_probe, "\n\n"
    )
  }
  
  if(!is.null(analysis_results$correlation)) {
    report_content <- paste0(report_content,
      "### 相关性分析结果\n", 
      "- 显著负相关探针: ", analysis_results$correlation$negative_correlations, "\n",
      "- 最佳相关性: ", analysis_results$correlation$best_correlation, "\n\n"
    )
  }
  
  if(!is.null(analysis_results$survival)) {
    report_content <- paste0(report_content,
      "### 生存分析结果\n",
      "- 显著预后探针: ", analysis_results$survival$prognostic_probes, "\n", 
      "- 最佳风险比: ", analysis_results$survival$best_hr, "\n\n"
    )
  }
  
  if(!is.null(analysis_results$colocalization)) {
    report_content <- paste0(report_content,
      "### 共定位分析结果\n",
      "- mQTL信号数量: ", analysis_results$colocalization$mqtl_signals, "\n",
      "- GWAS重叠: ", analysis_results$colocalization$gwas_overlap, "\n\n"
    )
  }
  
  # 保存报告
  report_file <- "methylation_analysis/comprehensive_methylation_report.md"
  writeLines(report_content, report_file)
  
  cat(sprintf("✓ 报告已保存至: %s\n", report_file))
  
  return(report_file)
}

# ============================================================================
# 5. 主执行函数
# ============================================================================

main_analysis <- function(run_steps = "all") {
  cat("🎬 甲基化重启分析主程序启动\n\n")
  
  tryCatch({
    # 环境设置
    probe_data <- setup_analysis_environment()
    
    # 数据检查
    available_data <- check_tcga_data_availability()
    
    # 运行分析流程
    results <- run_analysis_pipeline(
      target_genes = probe_data$target_genes,
      available_data = available_data
    )
    
    # 生成报告
    report_file <- generate_comprehensive_report(results)
    
    cat("\n🎉 甲基化重启分析完成！\n")
    cat("=======================\n")
    cat("📁 输出文件:\n")
    cat("- 综合报告: methylation_analysis/comprehensive_methylation_report.md\n")
    cat("- DMP结果: methylation_analysis/dmp_results/\n")
    cat("- 相关性结果: methylation_analysis/correlation_results/\n")
    cat("- 生存分析结果: methylation_analysis/survival_results/\n") 
    cat("- 共定位结果: methylation_analysis/colocalization/\n")
    cat("\n✅ 所有分析任务已成功完成！\n")
    
    return(list(results = results, report = report_file))
    
  }, error = function(e) {
    cat("❌ 分析过程中出现错误:\n")
    cat(paste0("错误信息: ", e$message, "\n"))
    cat("请检查数据路径和脚本文件\n")
    return(NULL)
  })
}

# ============================================================================
# 脚本直接执行入口
# ============================================================================

if(!interactive()) {
  # 命令行执行
  args <- commandArgs(trailingOnly = TRUE)
  run_steps <- if(length(args) > 0) args[1] else "all"
  
  cat("命令行参数:", run_steps, "\n")
  main_analysis(run_steps)
} else {
  # 交互式执行
  cat("💡 提示: 可以使用 main_analysis(\"all\") 开始完整分析\n")
  cat("         或指定步骤: main_analysis(c(\"dmp\", \"correlation\"))\n")
}