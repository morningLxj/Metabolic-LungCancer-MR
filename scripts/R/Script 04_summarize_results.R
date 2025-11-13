# ==============================================================================
# Script: 04_summarize_results.R (已修复 R 语法错误)
# Purpose: 汇总共定位结果并生成可视化
# ==============================================================================

cat("
================================================================================
Summarizing Colocalization Results
================================================================================

")

# ==============================================================================
# 加载包
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(ggrepel)
  library(ggsci)
  library(scales)
  library(logger)
})

# 设置日志
log_threshold(INFO)
log_appender(appender_tee("logs/04_summarize_results.log"))

log_info("Starting results summarization")

# ==============================================================================
# 配置参数
# ==============================================================================

RESULTS_DIR <- "results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
REPORTS_DIR <- file.path(RESULTS_DIR, "reports")

# <--- 修正: 自动创建尚不存在的目录 (与 01/02 脚本 [cite: 1, line 31; 2, line 32] 保持一致)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORTS_DIR, recursive = TRUE, showWarnings = FALSE)

# 共定位阈值
PP_H4_STRONG <- 0.8
PP_H4_MODERATE <- 0.5

# 图形主题 - 符合SCI期刊标准
theme_publication <- theme_bw(base_size = 14) +
  theme(
    # 网格和背景
    panel.grid.major = element_line(size = 0.4, color = "grey85"),
    panel.grid.minor = element_line(size = 0.2, color = "grey90"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    
    # 坐标轴和标题
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(3, "pt"),
    
    # 图例
    legend.position = "top",
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
    
    # 分面标题
    strip.background = element_rect(fill = "grey95", color = "black", size = 0.5),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    
    # 整体布局
    panel.border = element_rect(color = "black", size = 0.8, fill = NA),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(0, 0, 0.5, 0, "cm")),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0, 0, 0.5, 0, "cm")),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

theme_set(theme_publication)

# ==============================================================================
# 加载数据
# ==============================================================================

log_info("Loading colocalization results")

# 此路径与 03_run_coloc.R 的输出 (第 478 行) 完全匹配
coloc_results <- readRDS(file.path(TABLES_DIR, "coloc_results_full.rds"))

log_info(sprintf("  Total results: %d", nrow(coloc_results)))

# 只保留成功的分析
coloc_success <- coloc_results %>%
  filter(status == "success")

log_info(sprintf("  Successful analyses: %d", nrow(coloc_success)))

# ==============================================================================
# 1. 识别显著共定位的基因
# ==============================================================================

# <--- 修正: R 语法 (替换 Python 的 "..." * 80)
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Identifying Significant Colocalization")
log_info(paste(rep("=", 80), collapse = ""))

# 定义共定位类别
classify_coloc <- function(pp_h4) {
  case_when(
    pp_h4 >= PP_H4_STRONG ~ "Strong",
    pp_h4 >= PP_H4_MODERATE ~ "Moderate",
    TRUE ~ "Weak"
  )
}

coloc_classified <- coloc_success %>%
  mutate(
    # Exposure - eQTL 共定位强度
    exp_eqtl_strength = classify_coloc(exp_eqtl_PP.H4),
    # eQTL - Outcome 共定位强度
    eqtl_out_strength = classify_coloc(eqtl_out_PP.H4),
    # 整体评估：两个都要强
    overall_strength = case_when(
      exp_eqtl_PP.H4 >= PP_H4_STRONG & eqtl_out_PP.H4 >= PP_H4_STRONG ~ "Both Strong",
      exp_eqtl_PP.H4 >= PP_H4_STRONG | eqtl_out_PP.H4 >= PP_H4_STRONG ~ "One Strong",
      exp_eqtl_PP.H4 >= PP_H4_MODERATE & eqtl_out_PP.H4 >= PP_H4_MODERATE ~ "Both Moderate",
      TRUE ~ "Weak"
    )
  )

# 统计
cat("
Colocalization Strength Distribution:
")
table(coloc_classified$overall_strength) %>% print()

# 提取强共定位结果
strong_coloc <- coloc_classified %>%
  filter(overall_strength == "Both Strong") %>%
  arrange(desc(exp_eqtl_PP.H4 + eqtl_out_PP.H4))

log_info(sprintf("  Strong colocalization: %d gene-exposure-outcome combinations", 
                 nrow(strong_coloc)))

# 保存强共定位结果
fwrite(strong_coloc, 
       file.path(TABLES_DIR, "strong_colocalization.txt"),
       sep = "\t")

# 按基因汇总
gene_summary <- strong_coloc %>%
  group_by(gene) %>%
  summarise(
    n_exposures = n_distinct(exposure),
    n_outcomes = n_distinct(outcome),
    n_combinations = n(),
    exposures = paste(unique(exposure), collapse = ", "),
    outcomes = paste(unique(outcome), collapse = ", "),
    mean_exp_eqtl_PP.H4 = mean(exp_eqtl_PP.H4),
    mean_eqtl_out_PP.H4 = mean(eqtl_out_PP.H4),
    max_exp_eqtl_PP.H4 = max(exp_eqtl_PP.H4),
    max_eqtl_out_PP.H4 = max(eqtl_out_PP.H4),
    .groups = "drop"
  ) %>%
  arrange(desc(n_combinations))

log_info(sprintf("  Unique genes with strong colocalization: %d", nrow(gene_summary)))

fwrite(gene_summary,
       file.path(TABLES_DIR, "gene_summary_strong_coloc.txt"),
       sep = "\t")

cat("
Top 10 genes by number of strong colocalizations:
")
print(head(gene_summary, 10))

# ==============================================================================
# 2. 可视化：热图
# ==============================================================================

# <--- 修正: R 语法
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Generating Heatmaps")
log_info(paste(rep("=", 80), collapse = ""))

# 函数：创建共定位热图
create_coloc_heatmap <- function(data, value_col, title, filename) {
  
  # 准备矩阵数据 - 处理重复值
  mat_data <- data %>%
    select(gene, exposure, outcome, value = all_of(value_col)) %>%
    mutate(
      exp_out = paste(exposure, outcome, sep = "_"),
      value = as.numeric(value)  # 确保值为数值型
    ) %>%
    select(gene, exp_out, value) %>%
    # 处理重复的基因-暴露-结局组合，取最大值
    group_by(gene, exp_out) %>%
    summarise(value = max(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = exp_out, values_from = value) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # 手动将NA替换为0
  mat_data[is.na(mat_data)] <- 0
  
  # 只保留至少有一个强共定位的基因
  keep_genes <- rowSums(mat_data >= PP_H4_STRONG) > 0
  mat_data <- mat_data[keep_genes, , drop = FALSE]
  
  if (nrow(mat_data) == 0) {
    log_warn(sprintf("  No genes to plot for %s", title))
    return(NULL)
  }
  
  # 限制显示的基因数量
  if (nrow(mat_data) > 50) {
    # 选择PP.H4最高的50个基因
    top_genes <- names(sort(rowMeans(mat_data), decreasing = TRUE)[1:50])
    mat_data <- mat_data[top_genes, , drop = FALSE]
  }
  
  log_info(sprintf("  Creating heatmap: %s (%d genes)", title, nrow(mat_data)))
  
  # 绘制热图 - 符合期刊标准（使用系统默认字体）
  pdf(file.path(FIGURES_DIR, filename), width = 8, height = max(6, nrow(mat_data) * 0.25),
      pointsize = 12)
  
  pheatmap(
    mat_data,
    color = colorRampPalette(c("white", "#FFF7E6", "#FFB366", "#FF7043", "#D32F2F", "#B71C1C"))(100),
    breaks = seq(0, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 10,
    main = title,
    border_color = "grey60",
    cellwidth = 25,
    cellheight = 15,
    legend = TRUE,
    legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
    legend_labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")
  )
  
  dev.off()
  
  log_info(sprintf("  Saved: %s", filename))
}

# 创建热图：Exposure - eQTL
create_coloc_heatmap(
  coloc_classified,
  "exp_eqtl_PP.H4",
  "Colocalization: Exposure - Lung eQTL",
  "heatmap_exposure_eqtl.pdf"
)

# 创建热图：eQTL - Outcome
create_coloc_heatmap(
  coloc_classified,
  "eqtl_out_PP.H4",
  "Colocalization: Lung eQTL - Lung Cancer",
  "heatmap_eqtl_outcome.pdf"
)

# ==============================================================================
# 3. 可视化：散点图
# ==============================================================================

# <--- 修正: R 语法
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Generating Scatter Plots")
log_info(paste(rep("=", 80), collapse = ""))

# 散点图：Exposure-eQTL vs eQTL-Outcome PP.H4
p_scatter <- ggplot(coloc_classified, 
                    aes(x = exp_eqtl_PP.H4, y = eqtl_out_PP.H4)) +
  geom_point(aes(color = exposure, shape = outcome), 
             alpha = 0.6, size = 2) +
  geom_hline(yintercept = PP_H4_STRONG, linetype = "dashed", color = "red") +
  geom_vline(xintercept = PP_H4_STRONG, linetype = "dashed", color = "red") +
  geom_hline(yintercept = PP_H4_MODERATE, linetype = "dotted", color = "orange") +
  geom_vline(xintercept = PP_H4_MODERATE, linetype = "dotted", color = "orange") +
  scale_color_npg() +
  labs(
    x = "PP.H4 (Exposure - eQTL)",
    y = "PP.H4 (eQTL - Outcome)",
    title = "Colocalization Evidence",
    subtitle = "Red dashed: PP.H4 = 0.8 (strong); Orange dotted: PP.H4 = 0.5 (moderate)",
    color = "Exposure",
    shape = "Outcome"
  ) +
  theme_publication

# 保存高质量散点图 - 符合期刊标准
ggsave(file.path(FIGURES_DIR, "scatter_coloc_pp4.pdf"),
       p_scatter, width = 6, height = 5, dpi = 300, units = "in")

log_info("  Saved: scatter_coloc_pp4.pdf")

# 添加基因标签的版本（只标注强共定位）
p_scatter_labeled <- p_scatter +
  geom_text_repel(
    data = filter(coloc_classified,
                  exp_eqtl_PP.H4 > PP_H4_STRONG | eqtl_out_PP.H4 > PP_H4_STRONG),
    aes(label = gene),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.3,
    segment.size = 0.5,
    segment.color = "grey50"
  )

ggsave(file.path(FIGURES_DIR, "scatter_coloc_pp4_labeled.pdf"),
       p_scatter_labeled, width = 7, height = 5.5, dpi = 300, units = "in")

log_info("  Saved: scatter_coloc_pp4_labeled.pdf")

# ==============================================================================
# 4. 可视化：按暴露因素分面
# ==============================================================================

# <--- 修正: R 语法
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Generating Faceted Plots")
log_info(paste(rep("=", 80), collapse = ""))

# 按暴露因素分面的散点图
p_facet <- ggplot(coloc_classified,
                  aes(x = exp_eqtl_PP.H4, y = eqtl_out_PP.H4)) +
  geom_point(aes(color = outcome), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = PP_H4_STRONG, linetype = "dashed", 
             color = "red", size = 0.5) +
  geom_vline(xintercept = PP_H4_STRONG, linetype = "dashed", 
             color = "red", size = 0.5) +
  facet_wrap(~ exposure, ncol = 2) +
  scale_color_jco() +
  labs(
    x = "PP.H4 (Exposure - eQTL)",
    y = "PP.H4 (eQTL - Outcome)",
    title = "Colocalization by Exposure",
    color = "Outcome"
  ) +
  theme_publication

ggsave(file.path(FIGURES_DIR, "scatter_coloc_by_exposure.pdf"),
       p_facet, width = 7, height = 5.5, dpi = 300, units = "in")

log_info("  Saved: scatter_coloc_by_exposure.pdf")

# ==============================================================================
# 5. 可视化：条形图统计
# ==============================================================================

log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Generating Bar Plots")
log_info(paste(rep("=", 80), collapse = ""))

# 按暴露因素统计强共定位数量
exposure_counts <- coloc_classified %>%
  group_by(exposure, outcome) %>%
  summarise(
    total = n(),
    strong_exp_eqtl = sum(exp_eqtl_PP.H4 >= PP_H4_STRONG),
    strong_eqtl_out = sum(eqtl_out_PP.H4 >= PP_H4_STRONG),
    both_strong = sum(exp_eqtl_PP.H4 >= PP_H4_STRONG &
                      eqtl_out_PP.H4 >= PP_H4_STRONG),
    .groups = "drop"
  )

# 转换为长格式
exposure_counts_long <- exposure_counts %>%
  pivot_longer(
    cols = c(strong_exp_eqtl, strong_eqtl_out, both_strong),
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(
    category = factor(category,
                     levels = c("strong_exp_eqtl", "strong_eqtl_out", "both_strong"),
                     labels = c("Exposure-eQTL", "eQTL-Outcome", "Both"))
  )

# 高质量条形图 - 符合期刊标准
p_bar <- ggplot(exposure_counts_long,
                aes(x = exposure, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", size = 0.3) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3, fontface = "bold") +
  facet_wrap(~ outcome, ncol = 2) +
  scale_fill_manual(values = c("#E3F2FD", "#BBDEFB", "#1976D2")) +
  labs(
    x = "Exposure",
    y = "Number of Genes",
    title = "Strong Colocalization (PP.H4 ≥ 0.8)",
    fill = "Colocalization Type"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 11),
    legend.position = "top"
  )

ggsave(file.path(FIGURES_DIR, "barplot_strong_coloc_counts.pdf"),
       p_bar, width = 7, height = 5, dpi = 300, units = "in")

log_info("  Saved: barplot_strong_coloc_counts.pdf")

# ==============================================================================
# 6. 可视化：基因排名图
# ==============================================================================

# <--- 修正: R 语法
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Generating Gene Ranking Plots")
log_info(paste(rep("=", 80), collapse = ""))

# 为每个暴露-结局组合创建top基因图
create_top_genes_plot <- function(data, exp, out, top_n = 20) {
  
  plot_data <- data %>%
    filter(exposure == exp, outcome == out) %>%
    arrange(desc(exp_eqtl_PP.H4 + eqtl_out_PP.H4)) %>%
    head(top_n) %>%
    mutate(
      gene = factor(gene, levels = rev(gene)),
      combined_score = exp_eqtl_PP.H4 + eqtl_out_PP.H4
    )
  
  if (nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(y = gene)) +
    geom_segment(aes(x = 0, xend = exp_eqtl_PP.H4, yend = gene),
                 color = "steelblue", size = 3, alpha = 0.7) +
    geom_segment(aes(x = 0, xend = eqtl_out_PP.H4, yend = gene),
                 color = "coral", size = 3, alpha = 0.7,
                 position = position_nudge(y = 0.2)) +
    geom_vline(xintercept = PP_H4_STRONG, linetype = "dashed", color = "red") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = "PP.H4",
      y = "Gene",
      title = sprintf("Top %d Genes: %s -> %s", top_n, exp, out),
      subtitle = "Blue: Exposure-eQTL | Red: eQTL-Outcome"
    ) +
    theme_publication
  
  return(p)
}

# 为每个暴露-结局组合创建图
exposure_outcome_pairs <- coloc_classified %>%
  distinct(exposure, outcome)

for (i in 1:nrow(exposure_outcome_pairs)) {
  exp <- exposure_outcome_pairs$exposure[i]
  out <- exposure_outcome_pairs$outcome[i]
  
  p <- create_top_genes_plot(coloc_classified, exp, out, top_n = 20)
  
  if (!is.null(p)) {
    filename <- sprintf("top_genes_%s_%s.pdf", exp, out)
    ggsave(file.path(FIGURES_DIR, filename), p, width = 6, height = 8,
           dpi = 300, units = "in")
    log_info(sprintf("  Saved: %s", filename))
  }
}

# ==============================================================================
# 7. 创建综合表格
# ==============================================================================

# <--- 修正: R 语法
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Creating Summary Tables")
log_info(paste(rep("=", 80), collapse = ""))

# 表1：按暴露-结局汇总
summary_by_exp_out <- coloc_classified %>%
  group_by(exposure, outcome) %>%
  summarise(
    n_genes_tested = n(),
    n_strong_exp_eqtl = sum(exp_eqtl_PP.H4 >= PP_H4_STRONG),
    n_strong_eqtl_out = sum(eqtl_out_PP.H4 >= PP_H4_STRONG),
    n_both_strong = sum(exp_eqtl_PP.H4 >= PP_H4_STRONG & 
                        eqtl_out_PP.H4 >= PP_H4_STRONG),
    pct_both_strong = 100 * n_both_strong / n_genes_tested,
    mean_exp_eqtl_PP.H4 = mean(exp_eqtl_PP.H4),
    mean_eqtl_out_PP.H4 = mean(eqtl_out_PP.H4),
    median_exp_eqtl_PP.H4 = median(exp_eqtl_PP.H4),
    median_eqtl_out_PP.H4 = median(eqtl_out_PP.H4),
    .groups = "drop"
  ) %>%
  arrange(desc(n_both_strong))

fwrite(summary_by_exp_out,
       file.path(TABLES_DIR, "summary_by_exposure_outcome.txt"),
       sep = "\t")

cat("
Summary by Exposure-Outcome:
")
print(summary_by_exp_out)

# 表2：按基因汇总（所有基因）
summary_by_gene <- coloc_classified %>%
  group_by(gene) %>%
  summarise(
    n_tests = n(),
    n_exposures = n_distinct(exposure),
    n_outcomes = n_distinct(outcome),
    n_strong_exp_eqtl = sum(exp_eqtl_PP.H4 >= PP_H4_STRONG),
    n_strong_eqtl_out = sum(eqtl_out_PP.H4 >= PP_H4_STRONG),
    n_both_strong = sum(exp_eqtl_PP.H4 >= PP_H4_STRONG & 
                        eqtl_out_PP.H4 >= PP_H4_STRONG),
    max_exp_eqtl_PP.H4 = max(exp_eqtl_PP.H4),
    max_eqtl_out_PP.H4 = max(eqtl_out_PP.H4),
    mean_exp_eqtl_PP.H4 = mean(exp_eqtl_PP.H4),
    mean_eqtl_out_PP.H4 = mean(eqtl_out_PP.H4),
    .groups = "drop"
  ) %>%
  arrange(desc(n_both_strong), desc(max_exp_eqtl_PP.H4))

fwrite(summary_by_gene,
       file.path(TABLES_DIR, "summary_by_gene.txt"),
       sep = "\t")

cat("
Top 20 genes by strong colocalization:
")
print(head(summary_by_gene, 20))

# 表3：详细的强共定位结果（用于论文）- 符合期刊标准
paper_table <- strong_coloc %>%
  select(
    Gene = gene,
    Exposure = exposure,
    Outcome = outcome,
    Chr = chr,
    `N SNPs` = n_snps,
    `PP.H4 (Exp-eQTL)` = exp_eqtl_PP.H4,
    `PP.H4 (eQTL-Out)` = eqtl_out_PP.H4,
    `PP.H3 (Exp-eQTL)` = exp_eqtl_PP.H3,
    `PP.H3 (eQTL-Out)` = eqtl_out_PP.H3
  ) %>%
  mutate(
    `PP.H4 (Exp-eQTL)` = sprintf("%.3f", `PP.H4 (Exp-eQTL)`),
    `PP.H4 (eQTL-Out)` = sprintf("%.3f", `PP.H4 (eQTL-Out)`),
    `PP.H3 (Exp-eQTL)` = sprintf("%.3f", `PP.H3 (Exp-eQTL)`),
    `PP.H3 (eQTL-Out)` = sprintf("%.3f", `PP.H3 (eQTL-Out)`),
    `N SNPs` = as.integer(`N SNPs`)
  ) %>%
  arrange(Exposure, Outcome, desc(as.numeric(`PP.H4 (Exp-eQTL)`)))

# 期刊标准格式：制表符分隔，便于导入
fwrite(paper_table,
       file.path(TABLES_DIR, "table_strong_coloc_for_paper.txt"),
       sep = "\t", na = "")

# 期刊标准格式：逗号分隔，标准CSV
fwrite(paper_table,
       file.path(TABLES_DIR, "table_strong_coloc_for_paper.csv"),
       sep = ",", na = "")

# 创建LaTeX格式表格（供直接粘贴到论文中）
latex_table <- paper_table %>%
  mutate(
    Gene = paste0("\\textbf{", Gene, "}"),
    `PP.H4 (Exp-eQTL)` = paste0("\\textcolor{blue}{", `PP.H4 (Exp-eQTL)`, "}"),
    `PP.H4 (eQTL-Out)` = paste0("\\textcolor{red}{", `PP.H4 (eQTL-Out)`, "}")
  )

# 写入LaTeX格式 - 修复列数错误
latex_content <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Strong Colocalization Results (PP.H4 $\\geq$ 0.8)}\n",
  "\\label{tab:strong_coloc}\n",
  "\\begin{tabular}{lcccccc}\n",
  "\\toprule\n",
  "Gene & Exposure & Outcome & Chr & N SNPs & PP.H4 (Exp-eQTL) & PP.H4 (eQTL-Out) \\\\\n",
  "\\midrule\n",
  apply(latex_table[, -c(7, 8)], 1, function(row) {
    paste(paste0("  ", row), collapse = " & ")
  }) %>% paste(collapse = " \\\\\n"),
  " \\\\\n",
  "\\bottomrule\n",
  "\\end{tabular}\n",
  "\\end{table}\n"
)

writeLines(latex_content, file.path(TABLES_DIR, "table_strong_coloc_for_paper.tex"))

log_info("  Saved summary tables")

# ==============================================================================
# 8. 生成文本摘要
# ==============================================================================

# <--- 修正: R 语法
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Generating Text Summary")
log_info(paste(rep("=", 80), collapse = ""))

# 创建摘要报告
summary_text <- sprintf("
================================================================================
COLOCALIZATION ANALYSIS SUMMARY
================================================================================

Analysis Date: %s

OVERALL STATISTICS
------------------
Total analyses performed: %d
Successful analyses: %d (%.1f%%)
Failed analyses: %d (%.1f%%)

COLOCALIZATION RESULTS
----------------------
Strong colocalization (both PP.H4 ≥ 0.8): %d (%.1f%%)
Moderate colocalization (both PP.H4 ≥ 0.5): %d (%.1f%%)
Weak colocalization: %d (%.1f%%)

UNIQUE GENES
------------
Total genes tested: %d
Genes with strong colocalization: %d (%.1f%%)

TOP 10 GENES (by number of strong colocalizations):
%s

RESULTS BY EXPOSURE
-------------------
%s

RESULTS BY OUTCOME
------------------
%s

KEY FINDINGS
------------
1. %s exposure(s) showed strong evidence of colocalization
2. %d unique genes identified with strong colocalization
3. Most consistent findings for: %s

FILES GENERATED
---------------
Tables:
  - coloc_results_full.txt (all results)
  - strong_colocalization.txt (PP.H4 ≥ 0.8)
  - gene_summary_strong_coloc.txt (gene-level summary)
  - summary_by_exposure_outcome.txt (exposure-outcome summary)
  - summary_by_gene.txt (all genes summary)
  - table_strong_coloc_for_paper.txt (formatted for publication)

Figures:
  - heatmap_exposure_eqtl.pdf
  - heatmap_eqtl_outcome.pdf
  - scatter_coloc_pp4.pdf
  - scatter_coloc_pp4_labeled.pdf
  - scatter_coloc_by_exposure.pdf
  - barplot_strong_coloc_counts.pdf
  - top_genes_[exposure]_[outcome].pdf (multiple files)

================================================================================
",
  Sys.Date(),
  nrow(coloc_results),
  nrow(coloc_success),
  100 * nrow(coloc_success) / nrow(coloc_results),
  nrow(coloc_results) - nrow(coloc_success),
  100 * (nrow(coloc_results) - nrow(coloc_success)) / nrow(coloc_results),
  nrow(strong_coloc),
  100 * nrow(strong_coloc) / nrow(coloc_success),
  sum(coloc_classified$overall_strength == "Both Moderate"),
  100 * sum(coloc_classified$overall_strength == "Both Moderate") / nrow(coloc_success),
  sum(coloc_classified$overall_strength == "Weak"),
  100 * sum(coloc_classified$overall_strength == "Weak") / nrow(coloc_success),
  n_distinct(coloc_success$gene),
  nrow(gene_summary),
  100 * nrow(gene_summary) / n_distinct(coloc_success$gene),
  paste(capture.output(print(head(gene_summary[, c("gene", "n_combinations")], 10))), 
        collapse = "\n"),
  paste(capture.output(print(summary_by_exp_out[, c("exposure", "n_both_strong")])), 
        collapse = "\n"),
  paste(capture.output(
    coloc_classified %>% 
      group_by(outcome) %>% 
      summarise(n_strong = sum(overall_strength == "Both Strong")) %>%
      print()
  ), collapse = "\n"),
  sum(summary_by_exp_out$n_both_strong > 0),
  nrow(gene_summary),
  ifelse(nrow(gene_summary) > 0, gene_summary$gene[1], "None")
)

# 保存摘要
# <--- 修正: 使用我们创建的 REPORTS_DIR 变量
writeLines(summary_text, file.path(REPORTS_DIR, "analysis_summary.txt"))

cat(summary_text)

log_info("\n✓ Results summarization complete!")
log_info("Next step: Generate final report (05_generate_report.R)")