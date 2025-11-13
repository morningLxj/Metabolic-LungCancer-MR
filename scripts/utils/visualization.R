# =============================================================================
# 结果可视化和图表生成工具函数
# Results Visualization and Plotting Utilities
# =============================================================================

#' 创建MR分析结果森林图
#' Create Forest Plot for MR Analysis Results
#'
#' @param mr_results data.frame MR分析结果
#' @param method_name character MR方法名称
#' @param exposure_name character 暴露变量名称
#' @param outcome_name character 结果变量名称
#' @param plot_title character 图表标题
#' @param save_path character 保存路径
#'
#' @return ggplot对象
#' @export
#'
create_mr_forest_plot <- function(mr_results, 
                                 method_name = "MR-Egger",
                                 exposure_name = "Exposure",
                                 outcome_name = "Outcome",
                                 plot_title = NULL,
                                 save_path = NULL) {
  
  # 设置默认标题
  if (is.null(plot_title)) {
    plot_title <- sprintf("%s Method: %s on %s", 
                         method_name, exposure_name, outcome_name)
  }
  
  # 创建森林图数据
  if (!"ci_lower" %in% names(mr_results)) {
    mr_results$ci_lower <- mr_results$beta - 1.96 * mr_results$se
    mr_results$ci_upper <- mr_results$beta + 1.96 * mr_results$se
  }
  
  # 排序数据
  plot_data <- mr_results %>%
    arrange(beta) %>%
    mutate(SNP_id = paste0("SNP_", row_number())) %>%
    mutate(sig = ifelse(ci_lower * ci_upper > 0, "Significant", "Non-significant"))
  
  # 创建森林图
  forest_plot <- ggplot(plot_data, aes(x = beta, y = SNP_id)) +
    geom_point(aes(color = sig), size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = sig), 
                  height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = plot_title,
      subtitle = sprintf("OR = %.3f, 95%% CI: [%.3f, %.3f], p = %.2e", 
                       exp(mr_results$beta[1]), 
                       exp(mr_results$ci_lower[1]), 
                       exp(mr_results$ci_upper[1]), 
                       mr_results$pval[1]),
      x = sprintf("Log Odds Ratio (%s)", outcome_name),
      y = "Individual SNPs",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("Significant" = "red", "Non-significant" = "blue"))
  
  # 如果指定了保存路径，则保存图片
  if (!is.null(save_path)) {
    ggsave(save_path, plot = forest_plot, width = 10, height = 8, dpi = 300)
    cat(sprintf("森林图已保存至: %s\n", save_path))
  }
  
  return(forest_plot)
}

#' 创建漏斗图检验MR分析的多效性
#' Create Funnel Plot to Test for Pleiotropy in MR Analysis
#'
#' @param mr_data data.frame MR分析数据
#' @param method_name character MR方法名称
#' @param exposure_name character 暴露变量名称
#' @param outcome_name character 结果变量名称
#' @param save_path character 保存路径
#'
#' @return ggplot对象
#' @export
#'
create_mr_funnel_plot <- function(mr_data,
                                 method_name = "IVW",
                                 exposure_name = "Exposure", 
                                 outcome_name = "Outcome",
                                 save_path = NULL) {
  
  # 计算精确度（1/SE）
  mr_data$precision <- 1 / mr_data$se
  
  # 创建漏斗图
  funnel_plot <- ggplot(mr_data, aes(x = beta, y = precision)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    # 添加对称线
    stat_function(fun = function(x) abs(x) * max(mr_data$precision, na.rm = TRUE) / 
                  max(abs(mr_data$beta), na.rm = TRUE), 
                  color = "red", linetype = "dotted") +
    labs(
      title = sprintf("Funnel Plot: %s Method", method_name),
      subtitle = sprintf("%s on %s", exposure_name, outcome_name),
      x = "Effect Size (β)",
      y = "Precision (1/SE)",
      caption = "Symmetry around zero suggests no directional pleiotropy"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, face = "italic")
    )
  
  # 添加多效性检验统计量
  if ("intercept" %in% names(mr_data) && "intercept_se" %in% names(mr_data)) {
    intercept_test <- mr_data$intercept[1]
    intercept_pval <- 2 * pnorm(-abs(intercept_test / mr_data$intercept_se[1]))
    
    funnel_plot <- funnel_plot +
      annotate("text", x = Inf, y = -Inf, 
              label = sprintf("Egger intercept = %.3f\np = %.3f", 
                            intercept_test, intercept_pval),
              hjust = 1.1, vjust = -0.1, 
              size = 3, color = "darkgreen")
  }
  
  # 保存图片
  if (!is.null(save_path)) {
    ggsave(save_path, plot = funnel_plot, width = 8, height = 6, dpi = 300)
    cat(sprintf("漏斗图已保存至: %s\n", save_path))
  }
  
  return(funnel_plot)
}

#' 创建MR分析结果对比图
#' Create MR Results Comparison Plot
#'
#' @param results_list list 多个MR方法的结果列表
#' @param exposure_name character 暴露变量名称
#' @param outcome_name character 结果变量名称
#' @param save_path character 保存路径
#'
#' @return ggplot对象
#' @export
#'
create_mr_comparison_plot <- function(results_list,
                                     exposure_name = "Exposure",
                                     outcome_name = "Outcome", 
                                     save_path = NULL) {
  
  # 整合结果数据
  comparison_data <- data.frame()
  
  for (method_name in names(results_list)) {
    if (is.data.frame(results_list[[method_name]]) && 
        nrow(results_list[[method_name]]) > 0) {
      temp_data <- results_list[[method_name]][1, ] %>%
        mutate(Method = method_name) %>%
        mutate(ci_lower = beta - 1.96 * se,
               ci_upper = beta + 1.96 * se)
      comparison_data <- rbind(comparison_data, temp_data)
    }
  }
  
  if (nrow(comparison_data) == 0) {
    stop("没有有效的MR结果数据用于比较")
  }
  
  # 创建比较图
  comparison_plot <- ggplot(comparison_data, aes(x = Method, y = beta)) +
    geom_point(size = 3, color = "darkblue") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                 width = 0.2, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = sprintf("MR Methods Comparison: %s on %s", exposure_name, outcome_name),
      subtitle = sprintf("OR Range: [%.3f, %.3f]", 
                       exp(min(comparison_data$ci_lower)), 
                       exp(max(comparison_data$ci_upper))),
      x = "MR Method",
      y = sprintf("Log Odds Ratio (%s)", outcome_name),
      caption = "Error bars represent 95% confidence intervals"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, face = "italic")
    ) +
    # 添加OR值标签
    geom_text(aes(label = sprintf("OR = %.3f", exp(beta))), 
             vjust = -0.5, size = 3)
  
  # 保存图片
  if (!is.null(save_path)) {
    ggsave(save_path, plot = comparison_plot, width = 10, height = 6, dpi = 300)
    cat(sprintf("MR比较图已保存至: %s\n", save_path))
  }
  
  return(comparison_plot)
}

#' 创建功能注释可视化
#' Create Functional Annotation Visualization
#'
#' @param annotation_data data.frame 功能注释数据
#' @param plot_type character 图表类型 ("bar", "pie", "stacked_bar")
#' @param category_col character 类别列名
#' @param save_path character 保存路径
#'
#' @return ggplot对象
#' @export
#'
create_annotation_plot <- function(annotation_data,
                                 plot_type = "bar",
                                 category_col = "category",
                                 save_path = NULL) {
  
  # 统计各类别数量
  category_counts <- annotation_data %>%
    group_by(across(all_of(category_col))) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(desc(count))
  
  # 创建不同的图表类型
  if (plot_type == "bar") {
    plot <- ggplot(category_counts, aes(x = reorder(!!sym(category_col), -count), y = count)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      labs(
        title = "Functional Annotation Distribution",
        x = "Category",
        y = "Number of SNPs"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
  } else if (plot_type == "pie") {
    plot <- ggplot(category_counts, aes(x = "", y = count, fill = !!sym(category_col))) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      labs(
        title = "Functional Annotation Distribution",
        fill = "Category",
        x = NULL,
        y = NULL
      ) +
      theme_minimal() +
      theme(legend.position = "right")
      
  } else if (plot_type == "stacked_bar") {
    # 需要染色体信息
    if ("chr" %in% names(annotation_data)) {
      plot <- annotation_data %>%
        group_by(chr, across(all_of(category_col))) %>%
        summarise(count = n(), .groups = 'drop') %>%
        ggplot(aes(x = chr, y = count, fill = !!sym(category_col))) +
        geom_bar(stat = "identity") +
        labs(
          title = "Functional Annotation by Chromosome",
          x = "Chromosome",
          y = "Number of SNPs",
          fill = "Category"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      stop("需要染色体信息来创建堆叠条形图")
    }
  }
  
  # 保存图片
  if (!is.null(save_path)) {
    ggsave(save_path, plot = plot, width = 10, height = 6, dpi = 300)
    cat(sprintf("功能注释图已保存至: %s\n", save_path))
  }
  
  return(plot)
}

#' 创建QQ图和曼哈顿图
#' Create QQ Plot and Manhattan Plot
#'
#' @param gwas_data data.frame GWAS数据
#' @param pval_col character p值列名
#' @param snp_col character SNP列名
#' @param chr_col character 染色体列名
#' @param pos_col character 位置列名
#' @param save_dir character 保存目录
#'
#' @return list 包含QQ图和曼哈顿图的列表
#' @export
#'
create_gwas_plots <- function(gwas_data,
                            pval_col = "pval",
                            snp_col = "SNP", 
                            chr_col = "chr",
                            pos_col = "pos",
                            save_dir = ".") {
  
  plots <- list()
  
  # 1. QQ图
  observed_p <- gwas_data[[pval_col]]
  expected_p <- ppoints(length(observed_p))
  qq_data <- data.frame(
    expected = -log10(expected_p),
    observed = -log10(sort(observed_p))
  )
  
  # 计算 genomic inflation
  median_chi_sq <- median(qchisq(observed_p, df = 1, lower.tail = FALSE))
  lambda <- median_chi_sq / qchisq(0.5, df = 1)
  
  qq_plot <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(
      title = "QQ Plot",
      subtitle = sprintf("λ = %.3f", lambda),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  qq_path <- file.path(save_dir, "qq_plot.png")
  ggsave(qq_path, plot = qq_plot, width = 8, height = 6, dpi = 300)
  plots$qq_plot <- qq_plot
  
  # 2. 曼哈顿图
  if (all(c(chr_col, pos_col) %in% names(gwas_data))) {
    manhattan_data <- gwas_data %>%
      mutate(log_p = -log10(!!sym(pval_col))) %>%
      arrange(!!sym(chr_col), !!sym(pos_col)) %>%
      group_by(!!sym(chr_col)) %>%
      mutate(chr_pos = cumsum(!!sym(pos_col))) %>%
      ungroup()
    
    # 设置染色体颜色
    manhattan_data$color <- ifelse(as.numeric(manhattan_data[[chr_col]]) %% 2 == 0, 
                                  "even", "odd")
    
    manhattan_plot <- ggplot(manhattan_data, aes(x = chr_pos, y = log_p)) +
      geom_point(aes(color = color), alpha = 0.6, size = 1) +
      scale_color_manual(values = c("even" = "darkblue", "odd" = "darkgreen")) +
      facet_grid(~ get(chr_col), scales = "free_x", space = "free_x") +
      labs(
        title = "Manhattan Plot",
        x = "Chromosome Position",
        y = "-log10(p)"
      ) +
      theme_minimal() +
      theme(
        strip.text.x = element_text(size = 8),
        legend.position = "none"
      ) +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red")
    
    manhattan_path <- file.path(save_dir, "manhattan_plot.png")
    ggsave(manhattan_path, plot = manhattan_plot, width = 12, height = 6, dpi = 300)
    plots$manhattan_plot <- manhattan_plot
  }
  
  cat(sprintf("GWAS图表已保存至目录: %s\n", save_dir))
  return(plots)
}

#' 创建综合结果报告图
#' Create Comprehensive Results Report Plot
#'
#' @param results_data list 包含各种分析结果的数据列表
#' @param save_path character 保存路径
#'
#' @return ggplot对象
#' @export
#'
create_comprehensive_report_plot <- function(results_data, save_path = NULL) {
  
  # 创建多面板图
  layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  
  # 准备绘图数据
  plot_list <- list()
  
  # Panel 1: MR结果摘要
  if ("mr_results" %in% names(results_data) && 
      is.data.frame(results_data$mr_results) && 
      nrow(results_data$mr_results) > 0) {
    
    mr_summary <- results_data$mr_results %>%
      mutate(method = gsub("_", " ", Method)) %>%
      arrange(beta) %>%
      head(5)  # 只显示前5个结果
    
    p1 <- ggplot(mr_summary, aes(x = reorder(method, beta), y = beta)) +
      geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7) +
      coord_flip() +
      labs(title = "MR Method Results", x = NULL, y = "Effect Size") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    
    plot_list[[1]] <- p1
  }
  
  # Panel 2: 显著性SNP分布
  if ("significant_snps" %in% names(results_data)) {
    sig_counts <- data.frame(
      Analysis = names(results_data$significant_snps),
      Count = as.numeric(results_data$significant_snps)
    )
    
    p2 <- ggplot(sig_counts, aes(x = reorder(Analysis, -Count), y = Count)) +
      geom_bar(stat = "identity", fill = "lightcoral", alpha = 0.7) +
      labs(title = "Significant SNPs by Analysis", x = NULL, y = "Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    
    plot_list[[2]] <- p2
  }
  
  # Panel 3: 工具变量质量
  if ("iv_quality" %in% names(results_data)) {
    iv_metrics <- results_data$iv_quality
    
    p3 <- ggplot(iv_metrics, aes(x = metric, y = value)) +
      geom_bar(stat = "identity", fill = "lightgreen", alpha = 0.7) +
      labs(title = "Instrument Quality Metrics", x = NULL, y = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    
    plot_list[[3]] <- p3
  }
  
  # Panel 4: 分析摘要统计
  summary_stats <- data.frame(
    Metric = c("Total SNPs", "Valid IVs", "Significant Associations", "Heritability"),
    Value = c(
      results_data$total_snps %||% 0,
      results_data$valid_ivs %||% 0, 
      results_data$sig_associations %||% 0,
      results_data$heritability %||% 0
    )
  )
  
  p4 <- ggplot(summary_stats, aes(x = Metric, y = Value)) +
    geom_bar(stat = "identity", fill = "gold", alpha = 0.7) +
    labs(title = "Analysis Summary", x = NULL, y = "Count/Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  plot_list[[4]] <- p4
  
  # 创建综合图
  if (length(plot_list) > 0) {
    # 组合所有面板
    combined_plot <- grid.arrange(
      grobs = plot_list,
      layout_matrix = layout_matrix,
      top = "GWAS-MR Analysis Comprehensive Report",
      bottom = "Generated by Metabolic-LungCancer-MR Pipeline"
    )
  } else {
    # 如果没有数据，创建空白图
    combined_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No Results Available", 
              size = 12, hjust = 0.5) +
      theme_void() +
      labs(title = "GWAS-MR Analysis Comprehensive Report")
  }
  
  # 保存图片
  if (!is.null(save_path)) {
    ggsave(save_path, plot = combined_plot, width = 12, height = 10, dpi = 300)
    cat(sprintf("综合报告图已保存至: %s\n", save_path))
  }
  
  return(combined_plot)
}

# 辅助函数：安全提取列表元素（如果不存在返回NULL）
`%||%` <- function(x, y) if (!is.null(x)) x else y