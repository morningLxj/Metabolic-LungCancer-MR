#!/usr/bin/env Rscript

# =============================================================================
# 示例运行脚本 - 代谢肺癌MR分析
# Example Run Script - Metabolic Lung Cancer MR Analysis
# =============================================================================

# 加载必要的库
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
  library(forestplot)
  library(optparse)
})

# 加载工具函数
source("scripts/utils/data_processing.R")
source("scripts/utils/mr_analysis.R") 
source("scripts/utils/data_qc.R")
source("scripts/utils/visualization.R")

# 命令行参数设置
option_list <- list(
  make_option(c("--exposure", "-e"), type="character", default="BMI",
              help="暴露变量名称 [default: %default]"),
  make_option(c("--outcome", "-o"), type="character", default="Lung_Cancer",
              help="结果变量名称 [default: %default]"),
  make_option(c("--exposure_file", "-x"), type="character", default="data/raw/exposure_gwas.txt",
              help="暴露数据文件路径 [default: %default]"),
  make_option(c("--outcome_file", "-y"), type="character", default="data/raw/outcome_gwas.txt", 
              help="结果数据文件路径 [default: %default]"),
  make_option(c("--output", "-t"), type="character", default="results/",
              help="输出目录 [default: %default]"),
  make_option(c("--threads", "-p"), type="integer", default=1,
              help="并行线程数 [default: %default]"),
  make_option(c("--method", "-m"), type="character", default="all",
              help="MR方法选择 [default: %default]"),
  make_option(c("--sensitivity", "-s"), action="store_true", default=TRUE,
              help="是否进行敏感性分析 [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查参数
cat("=== 代谢肺癌MR分析参数 ===\n")
cat("暴露变量:", opt$exposure, "\n")
cat("结果变量:", opt$outcome, "\n")
cat("暴露数据文件:", opt$exposure_file, "\n")
cat("结果数据文件:", opt$outcome_file, "\n")
cat("输出目录:", opt$output, "\n")
cat("并行线程数:", opt$threads, "\n")
cat("MR方法:", opt$method, "\n")
cat("敏感性分析:", opt$sensitivity, "\n")
cat("================================\n\n")

# 创建输出目录
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output, "tables"), recursive = TRUE, showWarnings = FALSE)

# 示例数据创建函数（用于演示）
create_demo_data <- function(n_snps = 10000) {
  set.seed(123)
  
  # 创建示例暴露数据
  exposure_data <- data.frame(
    SNP = paste0("rs", 1:n_snps),
    chr = rep(1:22, each = ceiling(n_snps/22))[1:n_snps],
    pos = sample(1:250000000, n_snps, replace = TRUE),
    beta = rnorm(n_snps, 0, 0.1),
    se = runif(n_snps, 0.01, 0.15),
    pval = 10^(-runif(n_snps, 2, 8)),
    eaf = runif(n_snps, 0.05, 0.95),
    sample_size = rep(500000, n_snps),
    ncase = rep(0, n_snps),
    ncontrol = rep(500000, n_snps)
  )
  
  # 创建示例结果数据（与暴露有一定相关性）
  outcome_data <- exposure_data %>%
    mutate(
      beta = beta * 0.6 + rnorm(n_snps, 0, 0.08),
      se = se * runif(n_snps, 0.8, 1.2),
      pval = 10^(-runif(n_snps, 1.5, 7.5))
    )
  
  # 添加一些强相关SNPs
  strong_snps <- sample(1:n_snps, 100)
  exposure_data$pval[strong_snps] <- 10^(-runif(100, 5, 8))
  exposure_data$beta[strong_snps] <- rnorm(100, 0, 0.15)
  
  outcome_data$beta[strong_snps] <- exposure_data$beta[strong_snps] * 0.7 + rnorm(100, 0, 0.05)
  
  return(list(exposure = exposure_data, outcome = outcome_data))
}

# 数据加载函数
load_gwas_data <- function(exposure_file, outcome_file) {
  cat("正在加载GWAS数据...\n")
  
  # 检查文件是否存在，如果不存在则创建演示数据
  if (!file.exists(exposure_file) || !file.exists(outcome_file)) {
    cat("未找到指定的数据文件，创建演示数据...\n")
    demo_data <- create_demo_data()
    
    # 保存演示数据
    write.table(demo_data$exposure, exposure_file, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(demo_data$outcome, outcome_file, row.names = FALSE, sep = "\t", quote = FALSE)
    
    return(demo_data)
  }
  
  # 加载实际数据
  exposure_data <- read.table(exposure_file, header = TRUE, stringsAsFactors = FALSE)
  outcome_data <- read.table(outcome_file, header = TRUE, stringsAsFactors = FALSE)
  
  cat(sprintf("暴露数据: %d SNPs\n", nrow(exposure_data)))
  cat(sprintf("结果数据: %d SNPs\n", nrow(outcome_data)))
  
  return(list(exposure = exposure_data, outcome = outcome_data))
}

# 主要分析流程
main_analysis_pipeline <- function(exposure_data, outcome_data, exposure_name, outcome_name, 
                                 output_dir, method = "all", run_sensitivity = TRUE) {
  
  cat(sprintf("\n=== 开始主分析流程: %s -> %s ===\n", exposure_name, outcome_name))
  
  # 1. 数据质量控制
  cat("\n步骤1: 数据质量控制...\n")
  
  # 添加QC需要的列
  exposure_qc <- exposure_data %>%
    mutate(
      maf = pmin(eaf, 1 - eaf),
      hwe_p = runif(n(), 1e-10, 1),  # 演示数据用随机值
      info = runif(n(), 0.8, 1),
      missing_rate = runif(n(), 0, 0.05)
    )
  
  outcome_qc <- outcome_data %>%
    mutate(
      maf = pmin(eaf, 1 - eaf),
      hwe_p = runif(n(), 1e-10, 1),
      info = runif(n(), 0.8, 1),
      missing_rate = runif(n(), 0, 0.05)
    )
  
  # 执行QC
  qc_results <- list(
    exposure = gwas_qc_analysis(exposure_qc, exposure_name, outcome_name, 
                               exposure_data$sample_size[1]),
    outcome = gwas_qc_analysis(outcome_qc, exposure_name, outcome_name, 
                              outcome_data$sample_size[1])
  )
  
  # 2. 工具变量提取
  cat("\n步骤2: 提取工具变量...\n")
  
  strong_exposure <- exposure_data %>%
    filter(pval < 5e-8) %>%
    filter(maf >= 0.01) %>%
    filter(info >= 0.8)
  
  cat(sprintf("提取到 %d 个强工具变量\n", nrow(strong_exposure)))
  
  if (nrow(strong_exposure) < 3) {
    cat("警告: 工具变量数量不足，放宽阈值...\n")
    strong_exposure <- exposure_data %>%
      filter(pval < 1e-5) %>%
      filter(maf >= 0.01)
    cat(sprintf("放宽阈值后提取到 %d 个工具变量\n", nrow(strong_exposure)))
  }
  
  # 3. 数据协调
  cat("\n步骤3: 数据协调...\n")
  
  harmonized_data <- harmonise_data(
    strong_exposure,
    outcome_data
  )
  
  cat(sprintf("协调后保留 %d 个SNPs用于MR分析\n", nrow(harmonized_data)))
  
  if (nrow(harmonized_data) < 3) {
    stop("协调后的SNP数量不足，无法进行MR分析")
  }
  
  # 4. MR分析
  cat("\n步骤4: MR分析...\n")
  
  mr_results <- run_mr_analysis(harmonized_data, methods = method)
  
  # 5. 敏感性分析
  sensitivity_results <- NULL
  if (run_sensitivity && nrow(harmonized_data) >= 10) {
    cat("\n步骤5: 敏感性分析...\n")
    sensitivity_results <- run_sensitivity_analysis(harmonized_data)
  }
  
  # 6. 结果可视化
  cat("\n步骤6: 结果可视化...\n")
  
  # 创建森林图
  if (!is.null(mr_results) && nrow(mr_results) > 0) {
    forest_plot <- create_mr_forest_plot(
      mr_results, "IVW", exposure_name, outcome_name,
      save_path = file.path(output_dir, "plots", "forest_plot.png")
    )
    
    # 创建漏斗图
    funnel_plot <- create_mr_funnel_plot(
      harmonized_data, "IVW", exposure_name, outcome_name,
      save_path = file.path(output_dir, "plots", "funnel_plot.png")
    )
    
    # 创建方法比较图
    comparison_plot <- create_mr_comparison_plot(
      list("IVW" = mr_results), exposure_name, outcome_name,
      save_path = file.path(output_dir, "plots", "comparison_plot.png")
    )
  }
  
  # 7. 生成报告
  cat("\n步骤7: 生成分析报告...\n")
  
  # 保存结果表格
  if (!is.null(mr_results)) {
    write.csv(mr_results, file.path(output_dir, "tables", "mr_results.csv"), 
              row.names = FALSE)
  }
  
  if (!is.null(sensitivity_results)) {
    write.csv(sensitivity_results, file.path(output_dir, "tables", "sensitivity_results.csv"),
              row.names = FALSE)
  }
  
  # 生成QC报告
  qc_report_file <- file.path(output_dir, "data_qc_report.html")
  generate_qc_report(qc_results$exposure, qc_report_file)
  
  # 生成综合结果
  analysis_summary <- list(
    exposure_name = exposure_name,
    outcome_name = outcome_name,
    n_exposure_snps = nrow(exposure_data),
    n_outcome_snps = nrow(outcome_data),
    n_iv_snps = nrow(strong_exposure),
    n_harmonized_snps = nrow(harmonized_data),
    n_mr_methods = ifelse(is.null(mr_results), 0, nrow(mr_results)),
    significant_methods = ifelse(is.null(mr_results), 0, sum(mr_results$pval < 0.05)),
    main_ivw_result = ifelse(is.null(mr_results) || nrow(mr_results) == 0, 
                            NA, mr_results$or[1])
  )
  
  # 保存分析摘要
  saveRDS(analysis_summary, file.path(output_dir, "analysis_summary.rds"))
  
  cat("\n=== 分析完成 ===\n")
  cat(sprintf("结果已保存到: %s\n", output_dir))
  
  return(list(
    mr_results = mr_results,
    sensitivity = sensitivity_results,
    summary = analysis_summary
  ))
}

# 主函数
main <- function() {
  tryCatch({
    # 加载数据
    gwas_data <- load_gwas_data(opt$exposure_file, opt$outcome_file)
    
    # 运行主分析
    results <- main_analysis_pipeline(
      gwas_data$exposure, 
      gwas_data$outcome,
      opt$exposure,
      opt$outcome,
      opt$output,
      opt$method,
      opt$sensitivity
    )
    
    # 打印结果摘要
    cat("\n=== 分析结果摘要 ===\n")
    cat(sprintf("暴露变量: %s\n", opt$exposure))
    cat(sprintf("结果变量: %s\n", opt$outcome))
    cat(sprintf("分析SNPs: %d\n", results$summary$n_harmonized_snps))
    cat(sprintf("MR方法数: %d\n", results$summary$n_mr_methods))
    cat(sprintf("显著性方法: %d\n", results$summary$significant_methods))
    
    if (!is.na(results$summary$main_ivw_result)) {
      cat(sprintf("主要IVW结果 OR: %.3f\n", results$summary$main_ivw_result))
    }
    
    cat("分析完成！请查看输出目录获取详细结果。\n")
    
  }, error = function(e) {
    cat("错误:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}

# 运行主函数
if (!interactive()) {
  main()
} else {
  cat("此脚本需要在命令行运行，例如:\n")
  cat("Rscript run_example.R --exposure BMI --outcome Lung_Cancer\n")
  cat("或使用默认参数:\n")
  cat("Rscript run_example.R\n")
}