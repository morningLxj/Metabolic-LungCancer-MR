# 生成论文图表
cat("生成可视化图表...\n")

library(ggplot2)
library(forestplot)

# 生成森林图 (类似论文中的Figure 1)
generate_forest_plot <- function(mr_results) {
  # 筛选主要方法的结果
  main_results <- mr_results[mr_results$method == "Inverse variance weighted", ]
  
  # 准备森林图数据
  forest_data <- data.frame(
    trait = main_results$exposure,
    outcome = main_results$outcome,
    or = main_results$or,
    lci = main_results$or_lci95,
    uci = main_results$or_uci95,
    pvalue = main_results$pval
  )
  
  # 创建森林图
  p <- ggplot(forest_data, aes(x = or, y = trait, color = outcome)) +
    geom_point(position = position_dodge(width = 0.8)) +
    geom_errorbarh(aes(xmin = lci, xmax = uci), 
                   height = 0.2, position = position_dodge(width = 0.8)) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_x_log10() +
    labs(x = "Odds Ratio (log scale)", y = "Metabolic Trait") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# 生成异质性散点图 (类似论文中的Figure 2)
generate_heterogeneity_plot <- function(sensitivity_results) {
  heterogeneity_data <- data.frame()
  
  for (result_id in names(sensitivity_results)) {
    hetero <- sensitivity_results[[result_id]]$heterogeneity
    if (!is.null(hetero)) {
      temp_df <- data.frame(
        analysis = result_id,
        Q = hetero$Q[1],
        Q_pval = hetero$Q_pval[1]
      )
      heterogeneity_data <- rbind(heterogeneity_data, temp_df)
    }
  }
  
  p <- ggplot(heterogeneity_data, aes(x = Q, y = -log10(Q_pval))) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(x = "Cochran's Q statistic", y = "-log10(P-value)") +
    theme_minimal()
  
  return(p)
}

# 生成图表
forest_plot <- generate_forest_plot(all_mr_results)
heterogeneity_plot <- generate_heterogeneity_plot(sensitivity_results)

# 保存图表
ggsave("04_visualization/forest_plot.png", forest_plot, width = 12, height = 8)
ggsave("04_visualization/heterogeneity_plot.png", heterogeneity_plot, width = 8, height = 6)

cat("图表生成完成!\n")