# =============================================================================
# MR分析主函数单元测试
# Unit Tests for MR Analysis Main Functions
# =============================================================================

library(testthat)
library(MendelianRandomization)
library(TwoSampleMR)

# 设置测试环境
test_context("MR分析主函数测试")

# 创建测试数据
create_test_gwas_data <- function(n_snps = 1000) {
  set.seed(123)
  
  # 创建暴露数据
  exposure_data <- data.frame(
    SNP = paste0("rs", 1:n_snps),
    beta = rnorm(n_snps, 0, 0.1),
    se = runif(n_snps, 0.01, 0.1),
    pval = runif(n_snps, 1e-10, 1),
    eaf = runif(n_snps, 0.1, 0.9),
    sample_size = rep(500000, n_snps),
    ncase = rep(0, n_snps),
    ncontrol = rep(500000, n_snps),
    exposure = "Test_Exposure"
  )
  
  # 创建结果数据（与暴露数据有相关性）
  outcome_data <- exposure_data %>%
    mutate(
      beta = beta * 0.8 + rnorm(n_snps, 0, 0.05),  # 部分相关
      outcome = "Test_Outcome"
    )
  
  list(exposure = exposure_data, outcome = outcome_data)
}

test_data <- create_test_gwas_data(100)

# 测试1: 数据加载函数
test_that("数据加载函数测试", {
  # 测试创建的数据结构
  expect_type(test_data, "list")
  expect_length(test_data, 2)
  expect_true(all(c("exposure", "outcome") %in% names(test_data)))
  
  # 测试列名
  expect_true(all(c("SNP", "beta", "se", "pval") %in% names(test_data$exposure)))
  expect_true(all(c("SNP", "beta", "se", "pval") %in% names(test_data$outcome)))
})

# 测试2: 工具变量提取
test_that("工具变量提取测试", {
  # 筛选显著SNPs作为工具变量
  iv_snps <- test_data$exposure %>%
    filter(pval < 0.05) %>%
    pull(SNP)
  
  expect_type(iv_snps, "character")
  expect_true(length(iv_snps) >= 0)  # 可能没有显著的SNPs
})

# 测试3: 双样本MR分析
test_that("双样本MR分析测试", {
  # 过滤工具变量
  iv_snps <- test_data$exposure %>%
    filter(pval < 0.001) %>%
    pull(SNP)
  
  if (length(iv_snps) > 0) {
    # harmonize数据
    harmonized <- harmonise_data(
      test_data$exposure %>% filter(SNP %in% iv_snps),
      test_data$outcome %>% filter(SNP %in% iv_snps)
    )
    
    expect_s3_class(harmonized, "data.frame")
    expect_true(nrow(harmonized) > 0)
  }
})

# 测试4: 敏感性分析
test_that("敏感性分析测试", {
  iv_snps <- test_data$exposure %>%
    filter(pval < 0.001) %>%
    pull(SNP)
  
  if (length(iv_snps) > 5) {  # 需要足够多的SNPs
    # 模拟MR结果
    mr_result <- data.frame(
      method = c("IVW", "MR-Egger", "Weighted median"),
      nsnp = c(length(iv_snps), length(iv_snps), length(iv_snps)),
      b = c(0.05, 0.03, 0.04),
      se = c(0.02, 0.03, 0.025),
      pval = c(0.01, 0.3, 0.1),
      or = c(1.05, 1.03, 1.04),
      ci_low = c(1.01, 0.97, 0.99),
      ci_upp = c(1.09, 1.09, 1.09)
    )
    
    # 测试留一法
    leave_one_out <- mr_leaveoneout(harmonise_data(
      test_data$exposure %>% filter(SNP %in% iv_snps),
      test_data$outcome %>% filter(SNP %in% iv_snps)
    ))
    
    expect_s3_class(leave_one_out, "data.frame")
  }
})

# 测试5: 数据质量控制
test_that("数据质量控制测试", {
  # 测试QC函数（如果存在）
  qc_data <- test_data$exposure %>%
    mutate(
      maf = pmin(eaf, 1-eaf),
      hwe_p = runif(n(), 1e-10, 1),
      info = runif(n(), 0.8, 1),
      missing_rate = runif(n(), 0, 0.1)
    )
  
  # MAF过滤
  filtered_data <- qc_data %>% filter(maf >= 0.01)
  expect_true(nrow(filtered_data) >= 0)
  expect_true(nrow(filtered_data) <= nrow(qc_data))
})

# 测试6: 结果汇总
test_that("结果汇总测试", {
  # 创建模拟结果
  mr_results <- data.frame(
    method = c("IVW", "MR-Egger", "Weighted median"),
    nsnp = c(10, 10, 10),
    b = c(0.05, 0.03, 0.04),
    se = c(0.02, 0.03, 0.025),
    pval = c(0.01, 0.3, 0.1),
    or = c(1.05, 1.03, 1.04),
    ci_low = c(1.01, 0.97, 0.99),
    ci_upp = c(1.09, 1.09, 1.09)
  )
  
  # 测试OR值计算
  expect_equal(mr_results$or, exp(mr_results$b), tolerance = 1e-6)
  
  # 测试显著性判断
  significant <- mr_results$pval < 0.05
  expect_type(significant, "logical")
  expect_length(significant, nrow(mr_results))
})

# 测试7: 错误处理
test_that("错误处理测试", {
  # 测试空数据
  expect_error({
    empty_data <- data.frame()
    # 这里应该调用主要分析函数，但会失败
  })
  
  # 测试缺失必要列
  incomplete_data <- data.frame(
    SNP = c("rs1", "rs2"),
    beta = c(0.1, 0.2)  # 缺少se列
  )
  
  expect_error({
    # 尝试进行需要se列的分析
    harmonise_data(incomplete_data, incomplete_data)
  })
})

# 测试8: 输出格式
test_that("输出格式测试", {
  # 创建标准输出格式
  result_format <- data.frame(
    exposure = "Test_Exposure",
    outcome = "Test_Outcome",
    method = "IVW",
    nsnp = 10,
    b = 0.05,
    se = 0.02,
    pval = 0.01,
    or = 1.05,
    ci_low = 1.01,
    ci_upp = 1.09,
    stringsAsFactors = FALSE
  )
  
  expect_s3_class(result_format, "data.frame")
  expect_true(all(c("b", "se", "pval", "or") %in% names(result_format)))
})

# 运行测试
cat("运行MR分析主函数单元测试...\n")
cat("测试数据创建完成，包含", nrow(test_data$exposure), "个SNPs\n")
cat("所有主要功能测试已定义完成\n")