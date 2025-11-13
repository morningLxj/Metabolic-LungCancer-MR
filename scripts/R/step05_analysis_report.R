############################################################################
# 生成Step05单变量MR分析报告
# 基于 step05_mr_results_summary.csv 生成综合分析报告
############################################################################

cat("生成Step05单变量MR分析报告\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 1. 加载必要的包
required_packages <- c("dplyr", "openxlsx")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 2. 创建输出目录
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 3. 读取结果数据
results_file <- "step05_mr_results_summary.csv"
if (!file.exists(results_file)) {
  stop(sprintf("错误：找不到结果文件 %s\n请先运行 step05_单变量MR分析_完整修复版.R", results_file))
}

cat("正在读取结果文件...\n")
results_df <- read.csv(results_file, stringsAsFactors = FALSE)
cat(sprintf("✓ 成功读取 %d 个分析结果\n\n", nrow(results_df)))

# 4. 数据清理和准备
results_df$significant_fdr <- as.logical(results_df$significant_fdr)
results_df$significant_nominal <- as.logical(results_df$significant_nominal)
results_df$pval <- as.numeric(results_df$pval)
results_df$fdr_pval <- as.numeric(results_df$fdr_pval)
results_df$heterogeneity_p <- as.numeric(results_df$heterogeneity_p)
results_df$pleiotropy_p <- as.numeric(results_df$pleiotropy_p)

# 计算Bonferroni阈值
bonferroni_threshold <- 0.05 / nrow(results_df)
results_df$significant_bonferroni <- results_df$pval < bonferroni_threshold

# 5. 生成分析报告
report_file <- "results/step05_analysis_report.txt"
sink(report_file)

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("单变量孟德尔随机化（Univariable MR）分析综合报告\n")
cat("步骤5：代谢/炎症性状对肺癌的因果效应\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【一、分析概述】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat("1. 研究目的\n")
cat("   检验代谢性状和炎症标志物对肺癌（整体、肺腺癌、鳞状细胞肺癌）的因果效应\n\n")

cat("2. 分析方法\n")
cat("   • 暴露: 代谢指标（BMI、血脂、血糖等）和炎症标志物（CRP、IL6等）\n")
cat("   • 结局: 肺癌（整体肺癌、肺腺癌、鳞状细胞肺癌）\n")
cat("   • 方法: 双样本孟德尔随机化（Two-sample MR）\n")
cat("   • 主要分析: IVW（逆方差加权法）\n")
cat("   • 敏感性分析: MR-Egger、加权中位数、异质性、多效性检验\n\n")

# 统计暴露和结局数量
n_exposures <- length(unique(results_df$exposure))
n_outcomes <- length(unique(results_df$outcome))
total_analyses <- nrow(results_df)

cat("3. 数据源\n")
cat(sprintf("   • 暴露因子数: %d 个\n", n_exposures))
cat(sprintf("     - 代谢性状: %d 个\n", sum(results_df$category == "Metabolic", na.rm = TRUE) / n_outcomes))
cat(sprintf("     - 炎症标志物: %d 个\n", sum(results_df$category == "Inflammatory", na.rm = TRUE) / n_outcomes))
cat(sprintf("   • 结局数: %d 个\n", n_outcomes))
cat(sprintf("   • 总分析配对数: %d 个\n\n", total_analyses))

cat("【二、分析结果统计】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat(sprintf("总分析配对数: %d\n", total_analyses))
cat(sprintf("  • 成功分析: %d (100.0%%)\n", total_analyses))

# 显著性统计
nominal_sig <- sum(results_df$significant_nominal, na.rm = TRUE)
fdr_sig <- sum(results_df$significant_fdr, na.rm = TRUE)
bonf_sig <- sum(results_df$significant_bonferroni, na.rm = TRUE)

cat("\n显著性统计:\n")
cat(sprintf("  • 名义显著 (P<0.05): %d (%.1f%%)\n",
           nominal_sig,
           100 * nominal_sig / total_analyses))
cat(sprintf("  • FDR显著 (FDR<0.05): %d (%.1f%%)\n",
           fdr_sig,
           100 * fdr_sig / total_analyses))
cat(sprintf("  • Bonferroni显著 (P<%.2e): %d (%.1f%%)\n\n",
           bonferroni_threshold,
           bonf_sig,
           100 * bonf_sig / total_analyses))

cat("【三、关键发现】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

if (fdr_sig > 0) {
  cat("发现了代谢/炎症性状对肺癌的显著因果效应（FDR校正后显著）:\n\n")
  
  sig_results <- results_df %>%
    filter(significant_fdr == TRUE) %>%
    arrange(pval)
  
  for (i in 1:min(20, nrow(sig_results))) {
    row <- sig_results[i, ]
    cat(sprintf("%d. %s → %s\n", i, row$exposure, row$outcome))
    cat(sprintf("   OR = %s\n", row$or_95ci))
    cat(sprintf("   Beta = %.4f (SE = %.4f)\n", row$beta, row$se))
    cat(sprintf("   P = %.2e, FDR = %.2e\n", row$pval, row$fdr_pval))
    cat(sprintf("   工具变量: %d SNPs (协调后: %d), 策略 = %s\n",
               row$n_snps, row$n_harmonized, row$extraction_strategy))
    if (!is.na(row$heterogeneity_p)) {
      cat(sprintf("   异质性 P = %.2e", row$heterogeneity_p))
      if (row$heterogeneity_p < 0.05) {
        cat(" ⚠ (存在异质性)")
      }
      cat("\n")
    }
    if (!is.na(row$pleiotropy_p)) {
      cat(sprintf("   多效性 P = %.2e", row$pleiotropy_p))
      if (row$pleiotropy_p < 0.05) {
        cat(" ⚠ (存在多效性)")
      }
      cat("\n")
    }
    cat("\n")
  }
  
} else if (nominal_sig > 0) {
  cat("发现了名义显著的关联，但未达到FDR阈值:\n\n")
  
  nom_sig <- results_df %>%
    filter(significant_nominal == TRUE) %>%
    arrange(pval) %>%
    head(10)
  
  for (i in 1:nrow(nom_sig)) {
    row <- nom_sig[i, ]
    cat(sprintf("%d. %s → %s\n", i, row$exposure, row$outcome))
    cat(sprintf("   OR = %s, P = %.2e (FDR = %.2e)\n\n",
               row$or_95ci, row$pval, row$fdr_pval))
  }
  
} else {
  cat("未发现显著关联（P<0.05）\n")
  cat("所有分析的P值均 > 0.05\n\n")
}

# 按结局分组统计
cat("按结局分组的显著结果:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for (outcome in unique(results_df$outcome)) {
  outcome_data <- results_df[results_df$outcome == outcome, ]
  n_sig <- sum(outcome_data$significant_fdr, na.rm = TRUE)
  n_nom <- sum(outcome_data$significant_nominal, na.rm = TRUE)
  cat(sprintf("\n%s:\n", outcome))
  cat(sprintf("  • FDR显著: %d / %d (%.1f%%)\n",
             n_sig, nrow(outcome_data),
             100 * n_sig / nrow(outcome_data)))
  cat(sprintf("  • 名义显著: %d / %d (%.1f%%)\n",
             n_nom, nrow(outcome_data),
             100 * n_nom / nrow(outcome_data)))
}
cat("\n")

# 按暴露类别分组统计
cat("按暴露类别分组的显著结果:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for (category in unique(results_df$category)) {
  if (is.na(category)) next
  category_data <- results_df[results_df$category == category, ]
  n_sig <- sum(category_data$significant_fdr, na.rm = TRUE)
  n_nom <- sum(category_data$significant_nominal, na.rm = TRUE)
  cat(sprintf("\n%s性状:\n", category))
  cat(sprintf("  • FDR显著: %d / %d (%.1f%%)\n",
             n_sig, nrow(category_data),
             100 * n_sig / nrow(category_data)))
  cat(sprintf("  • 名义显著: %d / %d (%.1f%%)\n",
             n_nom, nrow(category_data),
             100 * n_nom / nrow(category_data)))
}
cat("\n")

cat("【四、质量控制结果】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat("1. 工具变量强度\n")
cat(sprintf("   • SNP数量 (平均±SD): %.1f ± %.1f\n",
           mean(results_df$n_snps, na.rm = TRUE),
           sd(results_df$n_snps, na.rm = TRUE)))
cat(sprintf("   • 协调后SNP数量 (平均±SD): %.1f ± %.1f\n",
           mean(results_df$n_harmonized, na.rm = TRUE),
           sd(results_df$n_harmonized, na.rm = TRUE)))
cat(sprintf("   • 协调成功率: %.1f%%\n\n",
           100 * mean(results_df$n_harmonized / results_df$n_snps, na.rm = TRUE)))

cat("2. 异质性检验\n")
het_count <- sum(results_df$heterogeneity_p < 0.05, na.rm = TRUE)
het_total <- sum(!is.na(results_df$heterogeneity_p))
cat(sprintf("   • 检测到显著异质性 (P<0.05): %d / %d (%.1f%%)\n",
           het_count, het_total,
           100 * het_count / het_total))
if (het_count > 0) {
  cat("   • 提示结果可能存在异质性，建议谨慎解释\n")
}
cat("\n")

cat("3. 水平多效性检验\n")
pleio_count <- sum(results_df$pleiotropy_p < 0.05, na.rm = TRUE)
pleio_total <- sum(!is.na(results_df$pleiotropy_p))
cat(sprintf("   • 检测到显著多效性 (P<0.05): %d / %d (%.1f%%)\n",
           pleio_count, pleio_total,
           100 * pleio_count / pleio_total))
if (pleio_count > 0) {
  cat("   • 提示可能存在水平多效性，需要进一步敏感性分析\n")
}
cat("\n")

cat("4. 工具变量提取策略\n")
strategy_stats <- table(results_df$extraction_strategy)
for (strategy in names(strategy_stats)) {
  cat(sprintf("   • %s: %d (%.1f%%)\n",
             strategy, strategy_stats[strategy],
             100 * strategy_stats[strategy] / nrow(results_df)))
}
cat("\n")

# 识别需要关注的效应方向
cat("5. 效应方向分析\n")
if (nrow(results_df) > 0) {
  protective_count <- sum(results_df$or < 1 & results_df$significant_fdr, na.rm = TRUE)
  risk_count <- sum(results_df$or > 1 & results_df$significant_fdr, na.rm = TRUE)
  
  if (fdr_sig > 0) {
    cat(sprintf("   • 保护性效应 (OR<1, FDR显著): %d 个\n", protective_count))
    cat(sprintf("   • 风险性效应 (OR>1, FDR显著): %d 个\n\n", risk_count))
  }
}

cat("【五、Top显著结果详细列表】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

if (fdr_sig > 0) {
  top_results <- results_df %>%
    filter(significant_fdr == TRUE) %>%
    arrange(pval) %>%
    head(15)
  
  cat("FDR校正后显著的前15个结果:\n\n")
  for (i in 1:nrow(top_results)) {
    row <- top_results[i, ]
    cat(sprintf("%d. %s → %s\n", i, row$exposure, row$outcome))
    cat(sprintf("   类别: %s\n", row$category))
    cat(sprintf("   OR (95%% CI) = %s\n", row$or_95ci))
    cat(sprintf("   Beta = %.4f, SE = %.4f\n", row$beta, row$se))
    cat(sprintf("   P值 = %.2e, FDR = %.2e\n", row$pval, row$fdr_pval))
    cat(sprintf("   SNPs: %d (协调后: %d)\n", row$n_snps, row$n_harmonized))
    cat(sprintf("   策略: %s\n", row$extraction_strategy))
    cat("\n")
  }
} else if (nominal_sig > 0) {
  top_results <- results_df %>%
    filter(significant_nominal == TRUE) %>%
    arrange(pval) %>%
    head(10)
  
  cat("名义显著的前10个结果:\n\n")
  for (i in 1:nrow(top_results)) {
    row <- top_results[i, ]
    cat(sprintf("%d. %s → %s\n", i, row$exposure, row$outcome))
    cat(sprintf("   OR (95%% CI) = %s\n", row$or_95ci))
    cat(sprintf("   P值 = %.2e (FDR = %.2e)\n\n", row$pval, row$fdr_pval))
  }
} else {
  cat("暂无显著结果\n\n")
}

cat("【六、结论与建议】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

if (fdr_sig > 0) {
  cat("主要结论:\n")
  cat("✓ 发现了代谢/炎症性状对肺癌的显著遗传因果效应\n")
  cat(sprintf("✓ 共识别出 %d 个FDR校正后显著的关联\n", fdr_sig))
  cat("✓ 结果提示某些代谢/炎症性状可能是肺癌的因果风险因素或保护因素\n\n")
  
  cat("建议:\n")
  cat("1. 对显著结果进行敏感性分析，评估结果的稳健性\n")
  cat("2. 使用独立的GWAS数据集进行验证\n")
  cat("3. 考虑进行多变量MR分析，控制多个暴露之间的混杂\n")
  cat("4. 探索显著关联的生物学机制\n")
  cat("5. 评估工具变量强度，确保结果不受弱工具变量偏倚影响\n")
  cat("6. 对于存在异质性或多效性的结果，需要谨慎解释\n\n")
  
} else if (nominal_sig > 0) {
  cat("主要结论:\n")
  cat("✓ 发现了一些名义显著的关联，但未通过FDR多重校正\n")
  cat(sprintf("✓ 共识别出 %d 个名义显著 (P<0.05) 的关联\n", nominal_sig))
  cat("✓ 提示可能存在潜在关联，但证据强度有限\n\n")
  
  cat("建议:\n")
  cat("1. 在更大的样本中验证这些名义显著的结果\n")
  cat("2. 考虑进行多变量MR分析以提高统计功效\n")
  cat("3. 结合其他证据（如功能研究、流行病学研究）进行综合判断\n")
  cat("4. 探索是否存在亚组特异性效应\n\n")
  
} else {
  cat("主要结论:\n")
  cat("○ 未发现显著关联（P<0.05）\n")
  cat("○ 提示所检测的代谢/炎症性状可能对肺癌无显著因果效应\n\n")
  
  cat("建议:\n")
  cat("1. 考虑扩大样本量或使用更精确的工具变量\n")
  cat("2. 探索是否存在非线性的因果关联\n")
  cat("3. 评估是否存在亚组特异性效应\n")
  cat("4. 考虑其他可能的暴露因子或结局定义\n\n")
}

cat("【七、方法学说明】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

cat("1. MR方法假设\n")
cat("   • 关联性假设: 工具变量与暴露相关\n")
cat("   • 独立性假设: 工具变量与混杂因素无关\n")
cat("   • 排他性假设: 工具变量仅通过暴露影响结局（无多效性）\n\n")

cat("2. 统计分析\n")
cat("   • 主要方法: IVW (Inverse Variance Weighted)\n")
cat("   • 敏感性分析: MR-Egger, Weighted Median\n")
cat("   • 多重校正: FDR (Benjamini-Hochberg), Bonferroni\n")
cat("   • 质量控制: 异质性检验 (Cochran's Q), 多效性检验 (MR-Egger intercept)\n\n")

cat("3. 结果解释\n")
cat("   • OR > 1: 暴露增加肺癌风险\n")
cat("   • OR < 1: 暴露降低肺癌风险\n")
cat("   • FDR < 0.05: 多重校正后显著\n")
cat("   • 异质性P < 0.05: 可能存在异质性\n")
cat("   • 多效性P < 0.05: 可能存在水平多效性\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("报告生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

sink()

cat("✓ 已保存分析报告: ", report_file, "\n\n")

# 6. 生成摘要统计表
cat("生成摘要统计表...\n")
summary_stats <- data.frame(
  metric = c(
    "总分析数",
    "名义显著 (P<0.05)",
    "FDR显著 (FDR<0.05)",
    "Bonferroni显著",
    "平均SNP数",
    "平均协调后SNP数",
    "显著异质性分析数",
    "显著多效性分析数"
  ),
  value = c(
    total_analyses,
    nominal_sig,
    fdr_sig,
    bonf_sig,
    round(mean(results_df$n_snps, na.rm = TRUE), 1),
    round(mean(results_df$n_harmonized, na.rm = TRUE), 1),
    het_count,
    pleio_count
  ),
  stringsAsFactors = FALSE
)

summary_file <- "results/step05_summary_statistics.txt"
write.table(summary_stats, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("✓ 已保存摘要统计: ", summary_file, "\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Step05分析报告生成完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("生成的文件:\n")
cat("  ✓ results/step05_analysis_report.txt - 详细分析报告\n")
cat("  ✓ results/step05_summary_statistics.txt - 摘要统计表\n\n")

cat("主要发现:\n")
if (fdr_sig > 0) {
  cat(sprintf("  ✓ 发现 %d 个FDR显著关联\n", fdr_sig))
  cat(sprintf("  ✓ 发现 %d 个名义显著关联\n", nominal_sig))
} else if (nominal_sig > 0) {
  cat(sprintf("  ○ 发现 %d 个名义显著关联（未通过FDR校正）\n", nominal_sig))
} else {
  cat("  ○ 未发现显著关联\n")
}

cat("\n请查看 results/step05_analysis_report.txt 获取详细分析结果。\n\n")

