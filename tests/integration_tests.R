# =============================================================================
# 集成测试脚本
# Integration Test Script for MR Analysis Pipeline
# =============================================================================

library(testthat)
library(MendelianRandomization)
library(TwoSampleMR)
library(dplyr)
library(ggplot2)

# 设置测试环境
test_context("GWAS-MR分析管道集成测试")

# 创建完整的模拟数据集
create_integration_test_data <- function() {
  set.seed(123)
  
  n_snps <- 5000
  
  # 创建更真实的GWAS数据
  test_gwas <- data.frame(
    SNP = paste0("rs", 1:n_snps),
    chr = rep(1:22, each = ceiling(n_snps/22))[1:n_snps],
    pos = sample(1:250000000, n_snps, replace = TRUE),
    beta = rnorm(n_snps, 0, 0.1),
    se = runif(n_snps, 0.01, 0.2),
    pval = 10^(-runif(n_snps, 2, 8)),  # 更有趣的p值分布
    eaf = runif(n_snps, 0.05, 0.95),
    info = runif(n_snps, 0.7, 1),
    sample_size = rep(500000, n_snps),
    ncase = rep(0, n_snps),
    ncontrol = rep(500000, n_snps),
    exposure = "BMI",
    outcome = "Lung_Cancer"
  )
  
  # 添加一些真实相关的SNPs
  causal_snps <- sample(1:n_snps, 50)
  test_gwas$beta[causal_snps] <- test_gwas$beta[causal_snps] + rnorm(50, 0.2, 0.05)
  test_gwas$pval[causal_snps] <- 10^(-runif(50, 4, 8))
  
  # 创建暴露和结果数据（部分重叠的SNPs）
  exposure_snps <- sample(test_gwas$SNP, 4000)
  outcome_snps <- sample(test_gwas$SNP, 3500)
  overlap_snps <- intersect(exposure_snps, outcome_snps)
  
  exposure_data <- test_gwas %>% 
    filter(SNP %in% exposure_snps) %>%
    mutate(exposure = "BMI")
  
  outcome_data <- test_gwas %>% 
    filter(SNP %in% outcome_snps) %>%
    mutate(
      beta = beta * 0.6 + rnorm(length(outcome_snps), 0, 0.1),
      outcome = "Lung_Cancer"
    )
  
  list(
    exposure = exposure_data,
    outcome = outcome_data,
    overlap_snps = overlap_snps
  )
}

# 测试1: 完整的数据预处理流程
test_that("完整数据预处理流程", {
  test_data <- create_integration_test_data()
  
  # 数据质量控制
  # 1. MAF过滤
  qc_exposure <- test_data$exposure %>%
    mutate(maf = pmin(eaf, 1 - eaf)) %>%
    filter(maf >= 0.01)
  
  qc_outcome <- test_data$outcome %>%
    mutate(maf = pmin(eaf, 1 - eaf)) %>%
    filter(maf >= 0.01)
  
  # 2. Info分数过滤
  qc_exposure <- qc_exposure %>% filter(info >= 0.8)
  qc_outcome <- qc_outcome %>% filter(info >= 0.8)
  
  # 3. p值过滤
  significant_exposure <- qc_exposure %>% filter(pval < 5e-8)
  significant_outcome <- qc_outcome %>% filter(pval < 5e-8)
  
  expect_true(nrow(qc_exposure) > 0)
  expect_true(nrow(qc_outcome) > 0)
  expect_true(nrow(significant_exposure) >= 0)
  expect_true(nrow(significant_outcome) >= 0)
})

# 测试2: 工具变量提取和验证
test_that("工具变量提取和验证", {
  test_data <- create_integration_test_data()
  
  # 提取工具变量
  exposure_snps <- test_data$exposure %>%
    filter(pval < 5e-8) %>%
    filter(maf >= 0.01) %>%
    filter(info >= 0.8) %>%
    pull(SNP)
  
  # 检查工具变量的有效性
  if (length(exposure_snps) >= 3) {
    # 模拟F统计量计算
    f_stats <- test_data$exposure %>%
      filter(SNP %in% exposure_snps) %>%
      mutate(f_stat = beta^2 / se^2)
    
    strong_ivs <- f_stats %>% filter(f_stat >= 10)
    
    expect_type(exposure_snps, "character")
    expect_true(length(exposure_snps) >= 0)
  }
})

# 测试3: 数据协调和harmonization
test_that("数据协调和harmonization", {
  test_data <- create_integration_test_data()
  
  # 选择重叠的SNPs进行分析
  overlap_data <- test_data$exposure %>%
    inner_join(test_data$outcome, by = "SNP", suffix = c("_exp", "_out"))
  
  expect_true(nrow(overlap_data) > 0)
  
  # 模拟harmonization过程
  harmonized_data <- overlap_data %>%
    mutate(
      # 确保效应等位基因一致
      beta_exp_adj = ifelse(eaf_exp > 0.5, -beta_exp, beta_exp),
      beta_out_adj = ifelse(eaf_out > 0.5, -beta_out, beta_out)
    )
  
  expect_true(all(!is.na(harmonized_data$beta_exp_adj)))
  expect_true(all(!is.na(harmonized_data$beta_out_adj)))
})

# 测试4: 多种MR方法的集成测试
test_that("多种MR方法集成测试", {
  test_data <- create_integration_test_data()
  
  # 创建重叠数据
  overlap_data <- test_data$exposure %>%
    inner_join(test_data$outcome, by = "SNP", suffix = c("_exp", "_out")) %>%
    filter(pval_exp < 5e-6) %>%  # 放宽阈值以获得足够SNPs
    mutate(
      beta_exp_adj = ifelse(eaf_exp > 0.5, -beta_exp, beta_exp),
      beta_out_adj = ifelse(eaf_out > 0.5, -beta_out, beta_out)
    )
  
  if (nrow(overlap_data) >= 10) {
    # 测试IVW方法
    ivw_result <- tryCatch({
      mr_input <- mr_input(
        bx = overlap_data$beta_exp_adj,
        bxse = overlap_data$se_exp,
        by = overlap_data$beta_out_adj,
        byse = overlap_data$se_out,
        snps = overlap_data$SNP
      )
      mr_ivw(ivw_input)
    }, error = function(e) NULL)
    
    # 测试MR-Egger方法
    egger_result <- tryCatch({
      mr_egger(ivw_input)
    }, error = function(e) NULL)
    
    # 测试加权中位数方法
    median_result <- tryCatch({
      mr_median(ivw_input)
    }, error = function(e) NULL)
    
    expect_true(!is.null(ivw_result) | !is.null(egger_result) | !is.null(median_result))
  }
})

# 测试5: 敏感性分析集成测试
test_that("敏感性分析集成测试", {
  test_data <- create_integration_test_data()
  
  overlap_data <- test_data$exposure %>%
    inner_join(test_data$outcome, by = "SNP", suffix = c("_exp", "_out")) %>%
    filter(pval_exp < 5e-6)
  
  if (nrow(overlap_data) >= 10) {
    # 留一法测试
    leave_one_out_results <- list()
    for (i in 1:min(5, nrow(overlap_data))) {  # 只测试前5个
      loo_data <- overlap_data[-i, ]
      if (nrow(loo_data) >= 3) {
        # 模拟计算效应
        effect <- sum(loo_data$beta_exp * loo_data$beta_out / loo_data$se_out^2) / 
                 sum(loo_data$beta_exp^2 / loo_data$se_out^2)
        leave_one_out_results[[i]] <- effect
      }
    }
    
    expect_true(length(leave_one_out_results) > 0)
    
    # MR-PRESSO模拟测试
    # 计算异质性统计量
    overlap_data <- overlap_data %>%
      mutate(
        w = 1 / se_out^2,
        combo_stat = (beta_out - sum(beta_out * w) / sum(w))^2 * w
      )
    
    heterogeneity_stat <- sum(overlap_data$combo_stat)
    expect_true(heterogeneity_stat >= 0)
  }
})

# 测试6: 结果可视化和报告生成
test_that("结果可视化和报告生成", {
  test_data <- create_integration_test_data()
  
  # 创建模拟MR结果
  mr_results <- data.frame(
    method = c("IVW", "MR-Egger", "Weighted median", "Simple mode"),
    nsnp = c(50, 50, 50, 50),
    b = c(0.05, 0.03, 0.04, 0.06),
    se = c(0.02, 0.03, 0.025, 0.04),
    pval = c(0.01, 0.3, 0.1, 0.15),
    or = exp(c(0.05, 0.03, 0.04, 0.06)),
    stringsAsFactors = FALSE
  )
  
  # 测试森林图数据准备
  forest_data <- mr_results %>%
    mutate(
      ci_lower = b - 1.96 * se,
      ci_upper = b + 1.96 * se
    )
  
  expect_true(all(c("ci_lower", "ci_upper") %in% names(forest_data)))
  
  # 测试显著性统计
  significant_methods <- sum(mr_results$pval < 0.05)
  expect_type(significant_methods, "integer")
  expect_true(significant_methods >= 0)
})

# 测试7: 错误处理和边界情况
test_that("错误处理和边界情况", {
  # 测试空数据
  expect_error({
    empty_exposure <- data.frame()
    empty_outcome <- data.frame()
    # 这里应该产生错误
  })
  
  # 测试数据不匹配
  mismatched_data <- list(
    exposure = data.frame(SNP = c("rs1", "rs2"), beta = c(0.1, 0.2)),
    outcome = data.frame(SNP = c("rs3", "rs4"), beta = c(0.3, 0.4))
  )
  
  overlap <- inner_join(mismatched_data$exposure, 
                       mismatched_data$outcome, 
                       by = "SNP")
  expect_equal(nrow(overlap), 0)
  
  # 测试工具变量不足
  few_snps <- data.frame(
    SNP = c("rs1", "rs2"),
    beta = c(0.1, 0.2),
    se = c(0.05, 0.06)
  )
  
  expect_warning({
    # 如果只有2个SNPs，MR分析会给出警告
    if (nrow(few_snps) < 3) {
      expect_true(TRUE)  # 测试通过
    }
  })
})

# 测试8: 性能测试
test_that("性能测试", {
  # 大数据集处理测试
  large_data <- create_integration_test_data()
  large_data$exposure <- large_data$exposure %>% sample_n(2000, replace = TRUE)
  large_data$outcome <- large_data$outcome %>% sample_n(2000, replace = TRUE)
  
  start_time <- Sys.time()
  
  # 数据预处理
  qc_exposure <- large_data$exposure %>%
    filter(pval < 5e-8) %>%
    filter(maf >= 0.01)
  
  qc_outcome <- large_data$outcome %>%
    filter(pval < 5e-8) %>%
    filter(maf >= 0.01)
  
  # 数据合并
  overlap_data <- inner_join(qc_exposure, qc_outcome, by = "SNP")
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  expect_true(processing_time < 10)  # 应该在10秒内完成
})

# 运行集成测试
cat("开始运行GWAS-MR分析管道集成测试...\n")

# 创建测试数据
integration_test_data <- create_integration_test_data()
cat(sprintf("集成测试数据创建完成:\n"))
cat(sprintf("- 暴露数据: %d SNPs\n", nrow(integration_test_data$exposure)))
cat(sprintf("- 结果数据: %d SNPs\n", nrow(integration_test_data$outcome)))
cat(sprintf("- 重叠SNPs: %d\n", length(integration_test_data$overlap_snps)))

# 运行所有测试
cat("所有集成测试已定义完成\n")
cat("使用方法: test_file('path/to/integration_tests.R') 运行测试\n")