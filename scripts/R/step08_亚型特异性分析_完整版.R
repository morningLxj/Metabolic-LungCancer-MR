############################################################################
# 第8步：亚型特异性分析（Subtype-Specific Analysis）
# 
# 功能：
#   1. 生成Table 6：汇总所有显著关联及其对腺癌、鳞癌、总体肺癌的效应
#   2. 生成亚型特异性森林图：并排展示每个显著性状对腺癌和鳞癌的效应
#   3. 生成热图：展示所有性状对三个结局亚型的效应矩阵
#
# 输出：
#   - Table 6: results/tables/paper_tables/Table6_Subtype_Specific_Associations.xlsx
#   - 森林图: results/figures/main_figures/FigureX_SubtypeSpecific_ForestPlot.pdf/png
#   - 热图: results/figures/main_figures/FigureX_SubtypeSpecific_Heatmap.pdf/png
############################################################################

cat("\n")
cat("=", rep("=", 78), "=\n", sep = "")
cat("第8步：亚型特异性分析（Subtype-Specific Analysis）\n")
cat("=", rep("=", 78), "=\n", sep = "")
cat("\n")

# ============================================================================
# 步骤1：加载必要的包和数据
# ============================================================================

cat("【步骤1】加载必要的包和数据\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 加载必要的包
required_packages <- c("dplyr", "openxlsx", "ggplot2", "gridExtra", 
                       "RColorBrewer", "pheatmap", "tidyr", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 创建输出目录
dirs <- c("results/tables/paper_tables", 
          "results/figures/main_figures",
          "results/figures/subtype_specific")
for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# 读取step05的结果汇总数据
results_file <- "step05_mr_results_summary.csv"
if (!file.exists(results_file)) {
  stop(sprintf("错误：找不到结果文件: %s\n请先运行step05脚本", results_file))
}

cat("✓ 读取结果数据...\n")
results_df <- read.csv(results_file, stringsAsFactors = FALSE)
cat(sprintf("  - 共 %d 条结果记录\n", nrow(results_df)))
cat(sprintf("  - 涉及 %d 个暴露因子\n", length(unique(results_df$exposure))))
cat(sprintf("  - 涉及 %d 个结局\n", length(unique(results_df$outcome))))

# ============================================================================
# 步骤2：准备数据和定义辅助函数
# ============================================================================

cat("\n【步骤2】准备数据和定义辅助函数\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 定义暴露因子标签映射（用于美化显示）
exposure_labels <- c(
  # 代谢性状
  "circulating_leptin" = "Circulating leptin",
  "vitamin_D" = "Vitamin D",
  "HbA1c" = "HbA1c",
  "ApoB" = "ApoB",
  "ApoA1" = "ApoA1",
  "IGF1" = "IGF-1",
  "ApoB_ApoA1_ratio" = "ApoB/ApoA1 ratio",
  "HDL_diameter" = "HDL diameter",
  "HDL_large" = "Large HDL",
  "remnant_cholesterol" = "Remnant cholesterol",
  "LDL_small" = "Small LDL",
  "BCAA" = "BCAA",
  "HDL_very_large" = "Very large HDL",
  "BMI" = "BMI",
  "HDL_cholesterol" = "HDL cholesterol",
  "LDL_cholesterol" = "LDL cholesterol",
  "smoking_initiation" = "Smoking initiation",
  "alcohol_drinks" = "Alcohol drinks per week",
  "fasting_glucose" = "Fasting glucose",
  "fasting_insulin" = "Fasting insulin",
  "SBP" = "Systolic BP",
  "DBP" = "Diastolic BP",
  "hypertension" = "Hypertension",
  "triglycerides" = "Triglycerides",
  "GGT" = "GGT",
  # 炎症标志物
  "CRP" = "C-reactive protein",
  "WBC" = "White blood cell count",
  "IL6" = "IL-6",
  "IL6R" = "IL-6R",
  "TNFR1" = "TNFR1"
)

# 定义结局标签映射
outcome_labels <- c(
  "Lung cancer" = "Overall lung cancer",
  "Lung adenocarcinoma" = "Lung adenocarcinoma",
  "Squamous cell lung cancer" = "Squamous cell lung cancer"
)

# 格式化P值的函数（向量化版本）
format_pvalue <- function(p) {
  # 处理缺失值
  result <- character(length(p))
  
  # 对于每个值单独处理
  for (i in seq_along(p)) {
    if (is.na(p[i])) {
      result[i] <- "NA"
    } else if (p[i] >= 0.05) {
      result[i] <- sprintf("%.3f", p[i])
    } else if (p[i] >= 0.001) {
      result[i] <- sprintf("%.3f", p[i])
    } else if (p[i] >= 0.0001) {
      result[i] <- sprintf("%.4f", p[i])
    } else {
      result[i] <- sprintf("%.2e", p[i])
    }
  }
  
  return(result)
}

# 格式化OR和CI的函数
format_or_ci <- function(or, lci, uci, digits = 3) {
  sprintf("%.*f (%.*f-%.*f)", digits, or, digits, lci, digits, uci)
}

cat("✓ 数据和函数已准备完成\n")

# ============================================================================
# 步骤3：识别显著关联的暴露因子
# ============================================================================

cat("\n【步骤3】识别显著关联的暴露因子\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 筛选显著关联（至少在一个结局中显著）
significant_exposures <- results_df %>%
  filter(significant_fdr == TRUE | significant_nominal == TRUE) %>%
  pull(exposure) %>%
  unique()

cat(sprintf("✓ 发现 %d 个显著暴露因子（至少在一个结局中显著）\n", 
            length(significant_exposures)))

if (length(significant_exposures) > 0) {
  cat("  显著暴露因子：\n")
  for (exp in significant_exposures) {
    exp_label <- ifelse(exp %in% names(exposure_labels), 
                        exposure_labels[exp], exp)
    cat(sprintf("    - %s\n", exp_label))
  }
} else {
  cat("⚠ 警告：未发现显著关联，将使用所有暴露因子进行分析\n")
  significant_exposures <- unique(results_df$exposure)
}

# ============================================================================
# 步骤4：生成Table 6 - 亚型特异性关联汇总表
# ============================================================================

cat("\n【步骤4】生成Table 6 - 亚型特异性关联汇总表\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 为每个显著暴露因子汇总三个结局的结果
table6_data <- list()

for (exp in significant_exposures) {
  exp_data <- results_df %>%
    filter(exposure == exp) %>%
    arrange(match(outcome, c("Lung cancer", "Lung adenocarcinoma", 
                             "Squamous cell lung cancer")))
  
  if (nrow(exp_data) > 0) {
    # 获取总体肺癌结果
    overall_lc <- exp_data %>% filter(outcome == "Lung cancer")
    if (nrow(overall_lc) == 0) overall_lc <- exp_data[1, ]
    
    # 获取腺癌结果
    adeno <- exp_data %>% filter(outcome == "Lung adenocarcinoma")
    if (nrow(adeno) == 0) adeno <- exp_data[1, ]
    
    # 获取鳞癌结果
    squamous <- exp_data %>% filter(outcome == "Squamous cell lung cancer")
    if (nrow(squamous) == 0) squamous <- exp_data[1, ]
    
    # 确定暴露因子类别
    category <- ifelse(nrow(overall_lc) > 0, overall_lc$category[1], 
                       ifelse(nrow(adeno) > 0, adeno$category[1], 
                              ifelse(nrow(squamous) > 0, squamous$category[1], "Unknown")))
    
    # 确定是否至少在一个结局中显著
    is_significant <- any(overall_lc$significant_fdr == TRUE, na.rm = TRUE) |
                     any(adeno$significant_fdr == TRUE, na.rm = TRUE) |
                     any(squamous$significant_fdr == TRUE, na.rm = TRUE) |
                     any(overall_lc$significant_nominal == TRUE, na.rm = TRUE) |
                     any(adeno$significant_nominal == TRUE, na.rm = TRUE) |
                     any(squamous$significant_nominal == TRUE, na.rm = TRUE)
    
    table6_data[[length(table6_data) + 1]] <- data.frame(
      # 基本信息
      Category = category,
      Exposure = exp,
      Exposure_Label = ifelse(exp %in% names(exposure_labels), 
                              exposure_labels[exp], exp),
      
      # 总体肺癌结果
      Overall_LC_OR = ifelse(nrow(overall_lc) > 0, overall_lc$or[1], NA),
      Overall_LC_OR_95CI = ifelse(nrow(overall_lc) > 0, overall_lc$or_95ci[1], NA),
      Overall_LC_P = ifelse(nrow(overall_lc) > 0, overall_lc$pval[1], NA),
      Overall_LC_FDR = ifelse(nrow(overall_lc) > 0, overall_lc$fdr_pval[1], NA),
      Overall_LC_Het_P = ifelse(nrow(overall_lc) > 0, overall_lc$heterogeneity_p[1], NA),
      Overall_LC_Significant = ifelse(nrow(overall_lc) > 0, 
                                       overall_lc$significant_fdr[1], FALSE),
      
      # 腺癌结果
      Adeno_OR = ifelse(nrow(adeno) > 0, adeno$or[1], NA),
      Adeno_OR_95CI = ifelse(nrow(adeno) > 0, adeno$or_95ci[1], NA),
      Adeno_P = ifelse(nrow(adeno) > 0, adeno$pval[1], NA),
      Adeno_FDR = ifelse(nrow(adeno) > 0, adeno$fdr_pval[1], NA),
      Adeno_Het_P = ifelse(nrow(adeno) > 0, adeno$heterogeneity_p[1], NA),
      Adeno_Significant = ifelse(nrow(adeno) > 0, adeno$significant_fdr[1], FALSE),
      
      # 鳞癌结果
      Squamous_OR = ifelse(nrow(squamous) > 0, squamous$or[1], NA),
      Squamous_OR_95CI = ifelse(nrow(squamous) > 0, squamous$or_95ci[1], NA),
      Squamous_P = ifelse(nrow(squamous) > 0, squamous$pval[1], NA),
      Squamous_FDR = ifelse(nrow(squamous) > 0, squamous$fdr_pval[1], NA),
      Squamous_Het_P = ifelse(nrow(squamous) > 0, squamous$heterogeneity_p[1], NA),
      Squamous_Significant = ifelse(nrow(squamous) > 0, 
                                     squamous$significant_fdr[1], FALSE),
      
      # 元信息
      Is_Significant = is_significant,
      
      stringsAsFactors = FALSE
    )
  }
}

if (length(table6_data) > 0) {
  table6 <- do.call(rbind, table6_data)
  
  # 按类别和显著性排序
  table6 <- table6 %>%
    arrange(Category, desc(Is_Significant), Overall_LC_P)
  
  # 创建用于Excel输出的格式化版本
  table6_formatted <- table6 %>%
    mutate(
      # 格式化总体肺癌列
      `Overall Lung Cancer` = paste0(
        ifelse(!is.na(Overall_LC_OR_95CI), Overall_LC_OR_95CI, "N/A"),
        ifelse(Overall_LC_Significant, " ***", ""),
        "\nP = ", format_pvalue(Overall_LC_P),
        ifelse(!is.na(Overall_LC_Het_P) & Overall_LC_Het_P < 0.05, 
               paste0(" [Het P = ", format_pvalue(Overall_LC_Het_P), "]"), "")
      ),
      
      # 格式化腺癌列
      `Lung Adenocarcinoma` = paste0(
        ifelse(!is.na(Adeno_OR_95CI), Adeno_OR_95CI, "N/A"),
        ifelse(Adeno_Significant, " ***", ""),
        "\nP = ", format_pvalue(Adeno_P),
        ifelse(!is.na(Adeno_Het_P) & Adeno_Het_P < 0.05, 
               paste0(" [Het P = ", format_pvalue(Adeno_Het_P), "]"), "")
      ),
      
      # 格式化鳞癌列
      `Squamous Cell Lung Cancer` = paste0(
        ifelse(!is.na(Squamous_OR_95CI), Squamous_OR_95CI, "N/A"),
        ifelse(Squamous_Significant, " ***", ""),
        "\nP = ", format_pvalue(Squamous_P),
        ifelse(!is.na(Squamous_Het_P) & Squamous_Het_P < 0.05, 
               paste0(" [Het P = ", format_pvalue(Squamous_Het_P), "]"), "")
      )
    ) %>%
    select(Category, Exposure_Label, 
           `Overall Lung Cancer`, `Lung Adenocarcinoma`, 
           `Squamous Cell Lung Cancer`)
  
  # 保存到Excel
  wb <- createWorkbook()
  
  # 详细数据表
  addWorksheet(wb, "Detailed Data")
  writeData(wb, "Detailed Data", table6, rowNames = FALSE)
  
  # 格式化表（用于论文）
  addWorksheet(wb, "Formatted Table")
  writeData(wb, "Formatted Table", table6_formatted, rowNames = FALSE)
  
  # 设置列宽
  if (ncol(table6) > 0) {
    setColWidths(wb, "Detailed Data", cols = seq_len(ncol(table6)), widths = "auto")
  }
  if (ncol(table6_formatted) > 0) {
    setColWidths(wb, "Formatted Table", cols = seq_len(ncol(table6_formatted)), 
                 widths = c(10, 25, 40, 40, 40))
  }
  
  # 保存文件
  excel_file <- "results/tables/paper_tables/Table6_Subtype_Specific_Associations.xlsx"
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  cat(sprintf("✓ Table 6 已生成: %s\n", excel_file))
  cat(sprintf("  - 包含 %d 个显著暴露因子\n", nrow(table6)))
  cat(sprintf("  - 详细数据表：%d 列，%d 行\n", ncol(table6), nrow(table6)))
  
  # 也保存CSV格式
  write.csv(table6, 
            "results/tables/paper_tables/Table6_Subtype_Specific_Associations.csv",
            row.names = FALSE)
  cat(sprintf("  - CSV格式: results/tables/paper_tables/Table6_Subtype_Specific_Associations.csv\n"))
  
} else {
  cat("⚠ 警告：没有显著关联，Table 6 为空\n")
}

# ============================================================================
# 步骤5：生成亚型特异性森林图
# ============================================================================

cat("\n【步骤5】生成亚型特异性森林图\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (length(significant_exposures) > 0 && exists("table6") && nrow(table6) > 0) {
  
  # 准备森林图数据
  forest_data <- data.frame()
  
  for (i in seq_len(nrow(table6))) {
    exp <- table6$Exposure[i]
    exp_label <- table6$Exposure_Label[i]
    category <- table6$Category[i]
    
    # 腺癌数据（从原始数据获取准确的CI）
    adeno_data_orig <- results_df %>%
      filter(exposure == exp, outcome == "Lung adenocarcinoma")
    if (nrow(adeno_data_orig) > 0 && !is.na(table6$Adeno_OR[i])) {
      forest_data <- rbind(forest_data, data.frame(
        Exposure = exp,
        Exposure_Label = exp_label,
        Category = category,
        Outcome = "Lung adenocarcinoma",
        OR = table6$Adeno_OR[i],
        OR_LCI = adeno_data_orig$or_lci[1],
        OR_UCI = adeno_data_orig$or_uci[1],
        P = table6$Adeno_P[i],
        Significant = table6$Adeno_Significant[i],
        stringsAsFactors = FALSE
      ))
    }
    
    # 鳞癌数据（从原始数据获取准确的CI）
    squamous_data <- results_df %>%
      filter(exposure == exp, outcome == "Squamous cell lung cancer")
    if (nrow(squamous_data) > 0 && !is.na(table6$Squamous_OR[i])) {
      forest_data <- rbind(forest_data, data.frame(
        Exposure = exp,
        Exposure_Label = exp_label,
        Category = category,
        Outcome = "Squamous cell lung cancer",
        OR = table6$Squamous_OR[i],
        OR_LCI = squamous_data$or_lci[1],
        OR_UCI = squamous_data$or_uci[1],
        P = table6$Squamous_P[i],
        Significant = table6$Squamous_Significant[i],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 移除缺失值
  forest_data <- forest_data %>%
    filter(!is.na(OR) & !is.na(OR_LCI) & !is.na(OR_UCI))
  
  if (nrow(forest_data) > 0) {
    
    # 设置因子顺序
    forest_data$Exposure_Label <- factor(forest_data$Exposure_Label,
                                        levels = rev(unique(forest_data$Exposure_Label)))
    forest_data$Outcome <- factor(forest_data$Outcome,
                                  levels = c("Lung adenocarcinoma", 
                                            "Squamous cell lung cancer"))
    
    # 创建森林图
    p_forest <- ggplot(forest_data, aes(x = OR, y = Exposure_Label, 
                                       color = Category, shape = Outcome)) +
      # 误差线
      geom_errorbarh(aes(xmin = OR_LCI, xmax = OR_UCI), 
                     height = 0.2, linewidth = 0.6, position = position_dodge(width = 0.6)) +
      # 点
      geom_point(size = 3, position = position_dodge(width = 0.6)) +
      # 显著标记
      geom_point(data = filter(forest_data, Significant == TRUE),
                aes(x = OR, y = Exposure_Label), 
                shape = 8, size = 2, color = "black", 
                position = position_dodge(width = 0.6)) +
      # 参考线
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      # 坐标轴
      scale_x_continuous(
        name = "Odds Ratio (95% CI)",
        trans = "log10",
        breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
        labels = c("0.5", "0.75", "1.0", "1.5", "2.0", "3.0"),
        limits = c(min(0.5, min(forest_data$OR_LCI, na.rm = TRUE) * 0.8),
                   max(3, max(forest_data$OR_UCI, na.rm = TRUE) * 1.2))
      ) +
      scale_color_manual(
        name = "Category",
        values = c("Metabolic" = "#D55E00", "Inflammatory" = "#0072B2"),
        guide = guide_legend(order = 1)
      ) +
      scale_shape_manual(
        name = "Subtype",
        values = c("Lung adenocarcinoma" = 16, "Squamous cell lung cancer" = 17),
        labels = c("Lung adenocarcinoma" = "Adenocarcinoma",
                   "Squamous cell lung cancer" = "Squamous"),
        guide = guide_legend(order = 2)
      ) +
      # 标签和主题
      labs(
        title = "Subtype-Specific Associations with Lung Cancer",
        subtitle = "Mendelian Randomization Analysis: Significant Exposures",
        y = NULL,
        caption = "Error bars represent 95% confidence intervals. Asterisks (*) denote FDR-adjusted significance (P < 0.05)."
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        plot.caption = element_text(hjust = 0, size = 9, color = "gray50"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
                panel.border = element_rect(linewidth = 0.8)
      )
    
    # 保存图片
    fig_width <- 12
    fig_height <- max(8, length(unique(forest_data$Exposure_Label)) * 0.5)
    
    ggsave("results/figures/main_figures/FigureX_SubtypeSpecific_ForestPlot.pdf",
           plot = p_forest, width = fig_width, height = fig_height,
           dpi = 600, limitsize = FALSE)
    ggsave("results/figures/main_figures/FigureX_SubtypeSpecific_ForestPlot.png",
           plot = p_forest, width = fig_width, height = fig_height,
           dpi = 600, limitsize = FALSE)
    
    cat(sprintf("✓ 亚型特异性森林图已生成\n"))
    cat(sprintf("  - PDF: results/figures/main_figures/FigureX_SubtypeSpecific_ForestPlot.pdf\n"))
    cat(sprintf("  - PNG: results/figures/main_figures/FigureX_SubtypeSpecific_ForestPlot.png\n"))
    cat(sprintf("  - 包含 %d 个暴露因子，%d 个数据点\n", 
                length(unique(forest_data$Exposure)), nrow(forest_data)))
    
  } else {
    cat("⚠ 警告：没有足够的数据生成森林图\n")
  }
  
} else {
  cat("⚠ 警告：没有显著关联，跳过森林图生成\n")
}

# ============================================================================
# 步骤6：生成亚型特异性热图
# ============================================================================

cat("\n【步骤6】生成亚型特异性热图\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 准备热图数据：所有暴露因子对三个结局的效应
heatmap_data <- results_df %>%
  select(exposure, outcome, or, pval, fdr_pval, significant_fdr, category) %>%
  filter(!is.na(or)) %>%
  mutate(
    log_or = log(or),
    significant = ifelse(significant_fdr == TRUE, "FDR<0.05",
                       ifelse(pval < 0.05, "P<0.05", "NS"))
  )

if (nrow(heatmap_data) > 0) {
  
  # 创建OR值矩阵（使用NA而不是0填充，避免影响范围计算）
  or_matrix <- heatmap_data %>%
    select(exposure, outcome, log_or) %>%
    pivot_wider(names_from = outcome, 
                values_from = log_or, 
                values_fill = NA) %>%
    arrange(match(exposure, c(rev(significant_exposures), 
                             setdiff(unique(exposure), significant_exposures))))
  
  # 创建显著性矩阵
  sig_matrix <- heatmap_data %>%
    select(exposure, outcome, significant) %>%
    mutate(sig_num = case_when(
      significant == "FDR<0.05" ~ 2,
      significant == "P<0.05" ~ 1,
      TRUE ~ 0
    )) %>%
    select(exposure, outcome, sig_num) %>%
    pivot_wider(names_from = outcome, 
                values_from = sig_num, 
                values_fill = 0)
  
  # 转换为矩阵
  or_mat <- as.matrix(or_matrix[, -1])
  rownames(or_mat) <- or_matrix$exposure
  
  sig_mat <- as.matrix(sig_matrix[, -1])
  rownames(sig_mat) <- sig_matrix$exposure
  
  # 检查并处理列名（pivot_wider可能将空格转换为点号）
  # 获取实际可用的列名
  available_cols <- colnames(or_mat)
  
  # 定义可能的列名映射（包括带空格和带点号的版本）
  col_mapping <- list(
    "Lung cancer" = "Overall\nLung Cancer",
    "Lung.cancer" = "Overall\nLung Cancer",
    "Lung adenocarcinoma" = "Lung\nAdenocarcinoma",
    "Lung.adenocarcinoma" = "Lung\nAdenocarcinoma",
    "Squamous cell lung cancer" = "Squamous Cell\nLung Cancer",
    "Squamous.cell.lung.cancer" = "Squamous Cell\nLung Cancer"
  )
  
  # 重新排序列：总体肺癌、腺癌、鳞癌
  # 尝试匹配多种可能的列名格式
  col_order <- NULL
  for (pattern in c("Lung.cancer", "Lung cancer", 
                    "Lung.adenocarcinoma", "Lung adenocarcinoma",
                    "Squamous.cell.lung.cancer", "Squamous cell lung cancer")) {
    if (any(grepl(gsub("\\.", "\\\\\\.", pattern), available_cols))) {
      matched_col <- available_cols[grepl(gsub("\\.", "\\\\\\.", pattern), available_cols)][1]
      if (!matched_col %in% col_order) {
        col_order <- c(col_order, matched_col)
      }
    }
  }
  
  # 如果找到了列，进行重排序
  if (length(col_order) > 0 && all(col_order %in% colnames(or_mat))) {
    or_mat <- or_mat[, col_order, drop = FALSE]
    sig_mat <- sig_mat[, col_order, drop = FALSE]
  }
  
  # 美化行名
  row_labels <- sapply(rownames(or_mat), function(x) {
    if (x %in% names(exposure_labels)) {
      exposure_labels[x]
    } else {
      gsub("_", " ", x)
    }
  })
  
  # 美化列名（应用映射）
  new_colnames <- sapply(colnames(or_mat), function(x) {
    if (x %in% names(col_mapping)) {
      col_mapping[[x]]
    } else {
      # 尝试匹配部分名称
      for (pattern in names(col_mapping)) {
        if (grepl(gsub("\\.", "\\\\\\.", pattern), x, ignore.case = TRUE)) {
          return(col_mapping[[pattern]])
        }
      }
      return(x)
    }
  })
  colnames(or_mat) <- new_colnames
  colnames(sig_mat) <- new_colnames
  
  # 确定颜色范围（处理NA、0和Inf情况）
  or_mat_no_na <- or_mat[!is.na(or_mat) & is.finite(or_mat)]
  if (length(or_mat_no_na) > 0 && any(or_mat_no_na != 0)) {
    or_range <- range(or_mat_no_na)
    max_abs <- max(abs(or_range))
    # 如果max_abs为0或无效，设置默认值
    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- 1.0
    }
  } else {
    # 如果没有有效值，设置默认范围
    max_abs <- 1.0
  }
  
  # 创建标注矩阵（星号）
  annotation_mat <- matrix("", nrow = nrow(or_mat), ncol = ncol(or_mat))
  rownames(annotation_mat) <- rownames(or_mat)
  colnames(annotation_mat) <- colnames(or_mat)
  
  if (nrow(sig_mat) > 0 && ncol(sig_mat) > 0) {
    for (i in seq_len(nrow(sig_mat))) {
      for (j in seq_len(ncol(sig_mat))) {
        if (sig_mat[i, j] == 2) {
          annotation_mat[i, j] <- "**"
        } else if (sig_mat[i, j] == 1) {
          annotation_mat[i, j] <- "*"
        }
      }
    }
  }
  
  # 生成热图（使用pheatmap）
  tryCatch({
    # 确保breaks参数有效
    heatmap_breaks <- tryCatch({
      seq(-max_abs, max_abs, length.out = 101)
    }, error = function(e) {
      # 如果失败，使用默认范围
      seq(-1, 1, length.out = 101)
    })
    
    # 验证breaks长度必须等于颜色数量+1
    if (length(heatmap_breaks) != 101) {
      heatmap_breaks <- seq(-max_abs, max_abs, length.out = 101)
    }
    
    pdf("results/figures/main_figures/FigureX_SubtypeSpecific_Heatmap.pdf",
        width = 8, height = max(10, nrow(or_mat) * 0.3))
    
    pheatmap(
      or_mat,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      cluster_rows = FALSE,  # 保持排序
      cluster_cols = FALSE,  # 保持列顺序
      scale = "none",
      main = "Subtype-Specific Associations with Lung Cancer\n(log Odds Ratio)",
      breaks = heatmap_breaks,
      display_numbers = annotation_mat,
      number_color = "black",
      number_format = "%.2f",
      fontsize = 9,
      fontsize_row = 8,
      fontsize_col = 10,
      fontsize_number = 10,
      border_color = "white",
      cellwidth = 40,
      cellheight = max(10, 300 / nrow(or_mat)),
      annotation_row = data.frame(
        Category = sapply(rownames(or_mat), function(x) {
          exp_data <- results_df %>% filter(exposure == x)
          if (nrow(exp_data) > 0) exp_data$category[1] else "Unknown"
        }),
        row.names = rownames(or_mat)
      ),
      annotation_colors = list(
        Category = c("Metabolic" = "#D55E00", "Inflammatory" = "#0072B2")
      )
    )
    
    dev.off()
    
    png("results/figures/main_figures/FigureX_SubtypeSpecific_Heatmap.png",
        width = 8, height = max(10, nrow(or_mat) * 0.3), 
        units = "in", res = 600)
    
    pheatmap(
      or_mat,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      scale = "none",
      main = "Subtype-Specific Associations with Lung Cancer\n(log Odds Ratio)",
      breaks = heatmap_breaks,
      display_numbers = annotation_mat,
      number_color = "black",
      number_format = "%.2f",
      fontsize = 9,
      fontsize_row = 8,
      fontsize_col = 10,
      fontsize_number = 10,
      border_color = "white",
      cellwidth = 40,
      cellheight = max(10, 300 / nrow(or_mat)),
      annotation_row = data.frame(
        Category = sapply(rownames(or_mat), function(x) {
          exp_data <- results_df %>% filter(exposure == x)
          if (nrow(exp_data) > 0) exp_data$category[1] else "Unknown"
        }),
        row.names = rownames(or_mat)
      ),
      annotation_colors = list(
        Category = c("Metabolic" = "#D55E00", "Inflammatory" = "#0072B2")
      )
    )
    
    dev.off()
    
    cat(sprintf("✓ 亚型特异性热图已生成\n"))
    cat(sprintf("  - PDF: results/figures/main_figures/FigureX_SubtypeSpecific_Heatmap.pdf\n"))
    cat(sprintf("  - PNG: results/figures/main_figures/FigureX_SubtypeSpecific_Heatmap.png\n"))
    cat(sprintf("  - 包含 %d 个暴露因子 × %d 个结局\n", 
                nrow(or_mat), ncol(or_mat)))
    cat(sprintf("  - 图例：** = FDR < 0.05, * = P < 0.05\n"))
    
  }, error = function(e) {
    cat(sprintf("⚠ 热图生成失败: %s\n", conditionMessage(e)))
  })
  
} else {
  cat("⚠ 警告：没有数据生成热图\n")
}

# ============================================================================
# 步骤7：生成总结报告
# ============================================================================

cat("\n【步骤7】生成总结报告\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 统计信息
if (exists("table6") && nrow(table6) > 0) {
  n_significant_overall <- sum(table6$Overall_LC_Significant, na.rm = TRUE)
  n_significant_adeno <- sum(table6$Adeno_Significant, na.rm = TRUE)
  n_significant_squamous <- sum(table6$Squamous_Significant, na.rm = TRUE)
  
  cat("\n【统计摘要】\n")
  cat(sprintf("  • 显著暴露因子总数：%d\n", nrow(table6)))
  cat(sprintf("  • 总体肺癌显著关联：%d\n", n_significant_overall))
  cat(sprintf("  • 腺癌显著关联：%d\n", n_significant_adeno))
  cat(sprintf("  • 鳞癌显著关联：%d\n", n_significant_squamous))
  
  # 识别亚型特异性关联（仅在某一亚型显著）
  subtype_specific <- table6 %>%
    filter(
      (Adeno_Significant == TRUE & Squamous_Significant == FALSE & 
       (is.na(Overall_LC_Significant) | Overall_LC_Significant == FALSE)) |
      (Adeno_Significant == FALSE & Squamous_Significant == TRUE & 
       (is.na(Overall_LC_Significant) | Overall_LC_Significant == FALSE))
    )
  
  if (nrow(subtype_specific) > 0) {
    cat(sprintf("\n  • 亚型特异性关联：%d\n", nrow(subtype_specific)))
    cat("    仅在一个亚型显著：\n")
    for (i in seq_len(nrow(subtype_specific))) {
      exp_label <- subtype_specific$Exposure_Label[i]
      if (subtype_specific$Adeno_Significant[i]) {
        cat(sprintf("      - %s (仅腺癌显著)\n", exp_label))
      }
      if (subtype_specific$Squamous_Significant[i]) {
        cat(sprintf("      - %s (仅鳞癌显著)\n", exp_label))
      }
    }
  }
}

cat("\n")
cat("=", rep("=", 78), "=\n", sep = "")
cat("第8步分析完成！\n")
cat("=", rep("=", 78), "=\n", sep = "")
cat("\n")

cat("【输出文件列表】\n")
cat("  1. Table 6: results/tables/paper_tables/Table6_Subtype_Specific_Associations.xlsx\n")
cat("  2. 森林图: results/figures/main_figures/FigureX_SubtypeSpecific_ForestPlot.pdf/png\n")
cat("  3. 热图: results/figures/main_figures/FigureX_SubtypeSpecific_Heatmap.pdf/png\n")
cat("\n")

