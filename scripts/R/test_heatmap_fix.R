# 测试热图修复
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# 加载数据
cat("Loading data...\n")
coloc_results <- readRDS("results/tables/coloc_results_full.rds")

# 过滤成功的分析
coloc_success <- coloc_results %>% filter(status == "success")

# 分类共定位
classify_coloc <- function(pp_h4) {
  case_when(
    pp_h4 >= 0.8 ~ "Strong",
    pp_h4 >= 0.5 ~ "Moderate", 
    TRUE ~ "Weak"
  )
}

coloc_classified <- coloc_success %>%
  mutate(
    exp_eqtl_strength = classify_coloc(exp_eqtl_PP.H4),
    eqtl_out_strength = classify_coloc(eqtl_out_PP.H4)
  )

cat("Data prepared. Testing heatmap function...\n")

# 测试修复后的热图函数
create_coloc_heatmap <- function(data, value_col, title, filename) {
  
  cat(sprintf("Creating heatmap: %s\n", title))
  
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
  
  cat(sprintf("Matrix dimensions: %d x %d\n", nrow(mat_data), ncol(mat_data)))
  
  # 手动将NA替换为0
  mat_data[is.na(mat_data)] <- 0
  
  # 只保留至少有一个强共定位的基因
  PP_H4_STRONG <- 0.8
  keep_genes <- rowSums(mat_data >= PP_H4_STRONG) > 0
  mat_data <- mat_data[keep_genes, , drop = FALSE]
  
  cat(sprintf("Genes with strong colocalization: %d\n", nrow(mat_data)))
  
  if (nrow(mat_data) == 0) {
    cat("No genes to plot\n")
    return(NULL)
  }
  
  # 限制显示的基因数量
  if (nrow(mat_data) > 50) {
    # 选择PP.H4最高的50个基因
    top_genes <- names(sort(rowMeans(mat_data), decreasing = TRUE)[1:50])
    mat_data <- mat_data[top_genes, , drop = FALSE]
  }
  
  cat(sprintf("Final matrix for heatmap: %d x %d\n", nrow(mat_data), ncol(mat_data)))
  
  # 创建简单的热图测试
  cat("Testing pheatmap with sample data...\n")
  
  # 创建临时测试图
  pdf("test_heatmap.pdf", width = 6, height = 4)
  
  pheatmap(
    mat_data,
    color = colorRampPalette(c("white", "#FFF7E6", "#FFB366", "#FF7043", "#D32F2F", "#B71C1C"))(100),
    breaks = seq(0, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 10,
    main = title
  )
  
  dev.off()
  cat("Test heatmap saved as test_heatmap.pdf\n")
  
  return(mat_data)
}

# 测试 Exposure - eQTL 热图
result <- create_coloc_heatmap(
  coloc_classified,
  "exp_eqtl_PP.H4",
  "Test: Exposure - Lung eQTL",
  "test_heatmap.pdf"
)

cat("Heatmap test completed successfully!\n")