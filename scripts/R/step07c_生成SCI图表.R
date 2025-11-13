############################################################################
# 步骤7c：生成SCI 10分期刊标准图表 - v4.3增强版
# 基于第7步中介分析结果生成4个关键图表
#
# 图表列表：
# 1. 中介路径网络图（代谢-炎症-肺癌路径网络）
# 2. 亚型特异性对比图（鳞癌vs腺癌的中介效应对比）
# 3. 效应大小分布图（所有路径的间接效应分布）
# 4. 成功vs失败路径对比图（中介分析挑战的可视化）
#
# v4.3更新（2025-11-02）：
# - 🔧 新增 label_dodge 和 label_push 参数，边标签自动避开节点
# - 🔧 解决边标签被节点遮挡的问题（如-438.1%和P值）
# - ⭐ 标签现在会智能地放置在节点旁边而非节点上方
############################################################################

cat("步骤7c：生成SCI 10分期刊标准图表\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 加载必要的包
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(tidyr)
  library(scales)
  library(openxlsx)
})

# ============================================================================
# SCI可视化主题函数库（已合并）
# ============================================================================

# Okabe-Ito 色盲友好调色板
okabe_ito_colors <- list(
  orange = "#E69F00",         # 代谢性状
  sky_blue = "#56B4E9",       # 炎症标志物
  green = "#009E73",          # 正向效应/结局
  vermillion = "#D55E00",     # 负向效应
  blue = "#0072B2",           # 辅助信息
  yellow = "#F0E442",         # 高亮
  reddish_purple = "#CC79A7", # 补充类别
  black = "#000000",          # 文本/强调
  gray = "#999999"            # 次要元素
)

okabe_ito <- okabe_ito_colors

# SCI标准主题函数
theme_sci <- function(base_size = 9, 
                     base_family = "Arial",
                     base_line_size = 0.21,  # 0.6 pt in mm
                     base_rect_size = 0.21,   # 0.6 pt in mm
                     grid = TRUE) {
  
  # 检查并设置字体（针对 Windows 系统优化）
  if (.Platform$OS.type == "windows") {
    # Windows 上，如果请求 Arial 但无法确认可用，使用系统默认字体
    if (base_family == "Arial") {
      font_available <- tryCatch({
        if (requireNamespace("extrafont", quietly = TRUE)) {
          "Arial" %in% extrafont::fonts()
        } else {
          FALSE  # extrafont 不可用，使用系统默认
        }
      }, error = function(e) FALSE)
      
      if (!font_available) {
        base_family <- ""  # 使用系统默认无衬线字体
      }
    }
  } else {
    # Unix 系统上，如果字体检查失败，使用 sans
    font_available <- tryCatch({
      if (requireNamespace("extrafont", quietly = TRUE)) {
        base_family %in% extrafont::fonts()
      } else {
        TRUE  # 假设字体可用
      }
    }, error = function(e) FALSE)
    
    if (!font_available && base_family != "") {
      base_family <- "sans"
    }
  }
  
  theme_classic(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    theme(
      # === 文本元素 ===
      plot.title = element_text(
        size = base_size + 2,  # 11 pt
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5, unit = "mm"),
        color = "black"
      ),
      plot.subtitle = element_text(
        size = base_size,  # 9 pt
        hjust = 0.5,
        margin = margin(b = 8, unit = "mm"),
        color = "black"
      ),
      plot.caption = element_text(
        size = base_size - 1,  # 8 pt
        hjust = 0.5,
        margin = margin(t = 5, unit = "mm"),
        color = "gray40"
      ),
      
      # === 轴标题 ===
      axis.title = element_text(
        size = base_size,  # 9 pt
        face = "bold",
        color = "black"
      ),
      axis.title.x = element_text(margin = margin(t = 5, unit = "mm")),
      axis.title.y = element_text(margin = margin(r = 5, unit = "mm"), angle = 90),
      
      # === 轴标签 ===
      axis.text = element_text(
        size = base_size - 1,  # 8 pt
        color = "black"
      ),
      axis.text.x = element_text(margin = margin(t = 2, unit = "mm")),
      axis.text.y = element_text(margin = margin(r = 2, unit = "mm")),
      
      # === 图例 ===
      legend.title = element_text(
        size = base_size,  # 9 pt
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = base_size - 1,  # 8 pt
        color = "black"
      ),
      legend.position = "right",
      legend.justification = "center",
      legend.box.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = base_line_size
      ),
      legend.margin = margin(5, 5, 5, 5, unit = "mm"),
      legend.spacing = unit(3, "mm"),
      
      # === 网格线 ===
      panel.grid.major = if (grid) {
        element_line(
          color = "gray90",
          linewidth = 0.11,  # 0.3 pt in mm
          linetype = "solid"
        )
      } else {
        element_blank()
      },
      panel.grid.minor = element_blank(),
      
      # === 轴线 ===
      axis.line = element_line(
        color = "black",
        linewidth = base_line_size  # 0.6 pt
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = base_line_size  # 0.6 pt
      ),
      axis.ticks.length = unit(2, "mm"),
      
      # === 面板背景 ===
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      
      # === 边距 ===
      plot.margin = margin(10, 10, 10, 10, unit = "mm"),
      
      # === 分面 ===
      strip.text = element_text(
        size = base_size,  # 9 pt
        face = "bold",
        color = "black"
      ),
      strip.background = element_rect(
        fill = "gray95",
        color = "black",
        linewidth = base_line_size
      ),
      
      # === 完整主题 ===
      complete = TRUE
    )
}

# 保存SCI标准图表函数
save_sci_figure <- function(plot, 
                           filename,
                           width_mm = 174,  # 双栏标准宽度
                           height_mm = 120,
                           dpi = 600,
                           formats = c("png", "pdf")) {
  
  # 确保目录存在
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  # 保存PNG版本（在线/预览）
  if ("png" %in% formats) {
    ggsave(
      paste0(filename, ".png"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi,
      bg = "white",
      type = "cairo"  # 更好的渲染质量
    )
    cat(sprintf("✓ PNG: %s.png (%.0f × %.0f mm, %d DPI)\n", 
                filename, width_mm, height_mm, dpi))
  }
  
  # 保存PDF版本（矢量，印刷用）
  if ("pdf" %in% formats) {
    ggsave(
      paste0(filename, ".pdf"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      device = "pdf",
      useDingbats = FALSE  # 避免Dingbats字体问题
    )
    cat(sprintf("✓ PDF: %s.pdf (%.0f × %.0f mm, 矢量格式)\n", 
                filename, width_mm, height_mm))
  }
  
  # 保存TIFF版本（印刷用，可选）
  if ("tiff" %in% formats || "tif" %in% formats) {
    ggsave(
      paste0(filename, ".tiff"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi,
      device = "tiff",
      compression = "lzw"  # 无损压缩
    )
    cat(sprintf("✓ TIFF: %s.tiff (%.0f × %.0f mm, %d DPI, LZW压缩)\n", 
                filename, width_mm, height_mm, dpi))
  }
}

# P值格式化函数（符合SCI期刊标准）
format_pvalue_sci <- function(p, digits = 3, format = "scientific") {
  if (format == "scientific") {
    # 科学计数法格式
    sapply(p, function(x) {
      if (is.na(x)) return("NA")
      if (x < 0.001) {
        sprintf("P < %.0e", 0.001)
      } else if (x < 0.01) {
        sprintf("P = %.2e", x)
      } else if (x < 0.05) {
        sprintf("P = %.3f", x)
      } else {
        sprintf("P = %.3f", x)
      }
    })
  } else {
    # 小数格式
    sapply(p, function(x) {
      if (is.na(x)) return("NA")
      if (x < 0.001) {
        "P < 0.001"
      } else {
        sprintf("P = %.3f", x)
      }
    })
  }
}

# 显著性标记函数
add_significance_stars <- function(p, fdr = NULL) {
  if (!is.null(fdr)) {
    # 使用FDR
    ifelse(fdr < 0.001, "***",
           ifelse(fdr < 0.01, "**",
                  ifelse(fdr < 0.05, "*", "ns")))
  } else {
    # 使用名义P值
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")))
  }
}

# 创建输出目录
dir.create("results/figures/step07_publication", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 步骤1：加载中介分析结果
# ============================================================================

cat("【步骤1】加载中介分析结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 尝试从不同位置加载结果
result_files <- c(
  "data/step07_all_mediation_results.RData",
  "results/tables/step07_mediation_results.csv",
  "results/tables/step07_mediation_results.xlsx"
)

mediation_summary <- NULL
all_pathways <- NULL

# 尝试加载RData文件
if (file.exists("data/step07_all_mediation_results.RData")) {
  load("data/step07_all_mediation_results.RData")
  cat("✓ 已从RData文件加载结果\n")
} else if (file.exists("results/tables/step07_mediation_results.csv")) {
  mediation_summary <- read.csv("results/tables/step07_mediation_results.csv", 
                               stringsAsFactors = FALSE)
  cat("✓ 已从CSV文件加载结果\n")
} else if (file.exists("results/tables/step07_mediation_results.xlsx")) {
  mediation_summary <- read.xlsx("results/tables/step07_mediation_results.xlsx")
  cat("✓ 已从Excel文件加载结果\n")
} else {
  stop("错误：找不到中介分析结果文件。请先运行步骤7。")
}

if (is.null(mediation_summary) || nrow(mediation_summary) == 0) {
  stop("错误：中介分析结果为空。")
}

cat(sprintf("✓ 已加载 %d 条中介路径分析结果\n\n", nrow(mediation_summary)))

# 如果没有all_pathways，尝试从数据中推断
if (is.null(all_pathways)) {
  cat("  提示：无法加载all_pathways，将从mediation_summary中推断\n")
}

# ============================================================================
# 数据诊断（检查数据是否满足生成图表的要求）
# ============================================================================

# 强制刷新输出缓冲区
flush.console()
cat("\n")
cat("【数据诊断】检查数据完整性\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 首先显示数据基本信息
cat("\n>>> 数据基本信息 <<<\n")
cat(sprintf("  数据行数: %d\n", nrow(mediation_summary)))
cat(sprintf("  数据列数: %d\n", ncol(mediation_summary)))
cat(sprintf("  列名: %s\n", paste(names(mediation_summary), collapse = ", ")))
cat("\n")

# 检查必需的列
cat(">>> 必需列检查 <<<\n")
required_cols <- c("exposure", "mediator", "outcome", "indirect_effect", "indirect_effect_pval")
missing_cols <- setdiff(required_cols, names(mediation_summary))
if (length(missing_cols) > 0) {
  cat(sprintf("  ⚠ 缺少必需的列: %s\n", paste(missing_cols, collapse = ", ")))
} else {
  cat("  ✓ 所有必需列都存在\n")
}
flush.console()

# 检查数据完整性
cat("\n>>> 数据完整性检查 <<<\n")
for (col in required_cols) {
  if (col %in% names(mediation_summary)) {
    na_count <- sum(is.na(mediation_summary[[col]]))
    cat(sprintf("  %s: %d 个NA值 (%.1f%%)\n", 
                col, na_count, 100 * na_count / nrow(mediation_summary)))
  } else {
    cat(sprintf("  %s: ❌ 列不存在\n", col))
  }
}
flush.console()

# 统计显著性
cat("\n>>> 显著性统计 <<<\n")
if ("indirect_effect_pval" %in% names(mediation_summary)) {
  valid_data <- mediation_summary %>% 
    filter(!is.na(indirect_effect), !is.na(indirect_effect_pval))
  cat(sprintf("  有效数据: %d 条\n", nrow(valid_data)))
  cat(sprintf("  P < 0.05: %d 条\n", sum(valid_data$indirect_effect_pval < 0.05, na.rm = TRUE)))
  cat(sprintf("  P < 0.1:  %d 条\n", sum(valid_data$indirect_effect_pval < 0.1, na.rm = TRUE)))
  cat(sprintf("  P < 0.2:  %d 条\n", sum(valid_data$indirect_effect_pval < 0.2, na.rm = TRUE)))
  
  if ("significant_fdr" %in% names(mediation_summary)) {
    cat(sprintf("  FDR显著:  %d 条\n", sum(mediation_summary$significant_fdr, na.rm = TRUE)))
  }
} else {
  cat("  ⚠ indirect_effect_pval 列不存在\n")
}
flush.console()

# 检查亚型数据
cat("\n>>> 亚型数据检查 <<<\n")
if ("outcome" %in% names(mediation_summary)) {
  outcomes <- unique(mediation_summary$outcome)
  cat(sprintf("  所有结局类型: %s\n", paste(outcomes, collapse = ", ")))
  for (out in outcomes) {
    n <- sum(mediation_summary$outcome == out, na.rm = TRUE)
    cat(sprintf("    %s: %d 条\n", out, n))
  }
  
  adeno <- mediation_summary %>% filter(outcome == "lung_adenocarcinoma")
  squamous <- mediation_summary %>% filter(outcome == "squamous_cell_lung")
  
  cat(sprintf("\n  腺癌数据: %d 条", nrow(adeno)))
  if (nrow(adeno) > 0 && "indirect_effect_pval" %in% names(adeno)) {
    cat(sprintf(" (P < 0.05: %d 条)", sum(adeno$indirect_effect_pval < 0.05, na.rm = TRUE)))
  }
  cat("\n")
  
  cat(sprintf("  鳞癌数据: %d 条", nrow(squamous)))
  if (nrow(squamous) > 0 && "indirect_effect_pval" %in% names(squamous)) {
    cat(sprintf(" (P < 0.05: %d 条)", sum(squamous$indirect_effect_pval < 0.05, na.rm = TRUE)))
  }
  cat("\n")
} else {
  cat("  ⚠ outcome 列不存在\n")
}
flush.console()

# 检查图表生成条件
cat("\n>>> 图表生成条件检查 <<<\n")

# 图表1：网络图
cat("  [图表1] 网络图 (中介路径网络图)\n")
tryCatch({
  if ("indirect_effect" %in% names(mediation_summary) && 
      "indirect_effect_pval" %in% names(mediation_summary)) {
    network_data_p05 <- mediation_summary %>%
      filter(!is.na(indirect_effect), indirect_effect_pval < 0.05)
    
    # 如果有FDR列，也包含FDR显著的结果
    if ("significant_fdr" %in% names(mediation_summary)) {
      network_data_p05 <- mediation_summary %>%
        filter(
          !is.na(indirect_effect),
          indirect_effect_pval < 0.05 | significant_fdr
        )
    }
    cat(sprintf("    P < 0.05或FDR显著: %d 条", nrow(network_data_p05)))
    if (nrow(network_data_p05) == 0) {
      network_data_p01 <- mediation_summary %>%
        filter(!is.na(indirect_effect), indirect_effect_pval < 0.1)
      cat(sprintf(" | P < 0.1: %d 条", nrow(network_data_p01)))
      if (nrow(network_data_p01) == 0) {
        cat(" ❌ 无法生成（无显著路径）")
      } else {
        cat(" ✓ 可以使用P < 0.1的数据生成")
      }
    } else {
      cat(" ✓ 可以生成")
    }
  } else {
    cat("    ❌ 缺少必需列（indirect_effect 或 indirect_effect_pval）")
  }
}, error = function(e) {
  cat(sprintf("    ❌ 检查出错: %s", conditionMessage(e)))
})
cat("\n")

# 图表2：亚型对比
cat("  [图表2] 亚型对比 (鳞癌vs腺癌)\n")
tryCatch({
  if ("outcome" %in% names(mediation_summary)) {
    subtype_data <- mediation_summary %>%
      filter(outcome %in% c("lung_adenocarcinoma", "squamous_cell_lung"))
    cat(sprintf("    亚型数据: %d 条", nrow(subtype_data)))
    if (nrow(subtype_data) > 0 && "indirect_effect_pval" %in% names(subtype_data)) {
      sig_subtype <- subtype_data %>% filter(indirect_effect_pval < 0.05)
      cat(sprintf(" | P < 0.05: %d 条", nrow(sig_subtype)))
      if (nrow(sig_subtype) == 0) {
        cat(" ⚠ 无显著路径，可能效果不佳")
      } else {
        cat(" ✓ 可以生成")
      }
    } else if (nrow(subtype_data) == 0) {
      cat(" ❌ 无法生成（无亚型数据）")
    } else {
      cat(" ⚠ 缺少p值列")
    }
  } else {
    cat("    ❌ outcome 列不存在")
  }
}, error = function(e) {
  cat(sprintf("    ❌ 检查出错: %s", conditionMessage(e)))
})
cat("\n")

# 图表3：效应分布
cat("  [图表3] 效应分布 (间接效应分布图)\n")
tryCatch({
  if ("indirect_effect" %in% names(mediation_summary)) {
    effect_dist_data <- mediation_summary %>% filter(!is.na(indirect_effect))
    cat(sprintf("    有间接效应数据: %d 条", nrow(effect_dist_data)))
    if (nrow(effect_dist_data) == 0) {
      cat(" ❌ 无法生成（无数据）")
    } else {
      cat(" ✓ 可以生成")
    }
  } else {
    cat("    ❌ indirect_effect 列不存在")
  }
}, error = function(e) {
  cat(sprintf("    ❌ 检查出错: %s", conditionMessage(e)))
})
cat("\n")

# 图表4：成功vs失败对比
cat("  [图表4] 成功vs失败对比\n")
tryCatch({
  if (exists("mediation_results") && !is.null(mediation_results)) {
    cat(sprintf("    ✓ 完整结果对象已加载（%s）\n", class(mediation_results)[1]))
    cat("    ✓ 可以生成")
  } else if (file.exists("data/step07_all_mediation_results.RData")) {
    cat("    ⚠ 文件存在但对象未加载，将在图表4部分尝试加载")
  } else {
    cat("    ❌ 完整结果对象未加载且文件不存在，无法生成详细对比")
  }
}, error = function(e) {
  cat(sprintf("    ❌ 检查出错: %s", conditionMessage(e)))
})
cat("\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("【诊断完成】如果看到上述信息，说明诊断代码已执行\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
flush.console()

# ============================================================================
# 变量名到英文标签的映射（确保图表中全部使用英文）
# ============================================================================

# 暴露因子英文标签映射
exposure_label_mapping <- c(
  "BMI" = "BMI",
  "fasting_insulin" = "Fasting Insulin",
  "fasting_glucose" = "Fasting Glucose",
  "HbA1c" = "HbA1c",
  "HDL_cholesterol" = "HDL Cholesterol",
  "LDL_cholesterol" = "LDL Cholesterol",
  "triglycerides" = "Triglycerides",
  "ApoB" = "ApoB",
  "ApoA1" = "ApoA1",
  "HDL_large" = "Large HDL",
  "HDL_diameter" = "HDL Diameter",
  "HDL_very_large" = "Very Large HDL",
  "LDL_small" = "Small LDL",
  "remnant_cholesterol" = "Remnant Cholesterol",
  "ApoB_ApoA1_ratio" = "ApoB/ApoA1",
  "SBP" = "SBP",
  "DBP" = "DBP",
  "hypertension" = "Hypertension",
  "smoking_initiation" = "Smoking Initiation",
  "alcohol_drinks" = "Alcohol Consumption",
  "IGF1" = "IGF-1",
  "circulating_leptin" = "Leptin",
  "vitamin_D" = "Vitamin D",
  "GGT" = "GGT",
  "BCAA" = "BCAA"
)

# 中介因子英文标签映射
mediator_label_mapping <- c(
  "CRP" = "CRP",
  "IL6" = "IL-6",
  "IL6R" = "IL-6R",
  "TNFR1" = "TNF-α",
  "WBC" = "WBC",
  "IL1RA" = "IL-1RA",
  "TNF" = "TNF",
  "TNFα" = "TNF-α"
)

# 结局变量英文标签映射
outcome_label_mapping <- c(
  "lung_cancer_overall" = "Overall\nLung Cancer",
  "lung_adenocarcinoma" = "Lung\nAdenocarcinoma",
  "squamous_cell_lung" = "Squamous Cell\nLung Cancer",
  "small_cell_lung" = "Small Cell\nLung Cancer"
)

# 函数：获取英文标签（支持向量输入）
get_english_label <- function(var_name, mapping) {
  # 向量化处理：逐个查找映射
  result <- sapply(var_name, function(vn) {
    if (vn %in% names(mapping)) {
      return(mapping[[vn]])
    } else {
      return(vn)
    }
  })
  return(as.character(result))
}

# ============================================================================
# 图表1：中介路径网络图（代谢-炎症-肺癌路径网络）
# ============================================================================

cat("【图表1】生成中介路径网络图\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 检查必需的列是否存在
required_cols_network <- c("exposure", "mediator", "outcome", "indirect_effect", "indirect_effect_pval")
missing_cols_network <- setdiff(required_cols_network, names(mediation_summary))
if (length(missing_cols_network) > 0) {
  stop(sprintf("错误：缺少网络图必需的列: %s", paste(missing_cols_network, collapse = ", ")))
}

# 检查可选的列
has_exp_to_med_beta <- "exp_to_med_beta" %in% names(mediation_summary)
has_med_to_out_beta <- "med_to_out_beta" %in% names(mediation_summary)
has_exp_to_med_pval <- "exp_to_med_pval" %in% names(mediation_summary)
has_med_to_out_pval <- "med_to_out_pval" %in% names(mediation_summary)
has_significant_fdr <- "significant_fdr" %in% names(mediation_summary)
has_mediation_proportion <- "mediation_proportion_percent" %in% names(mediation_summary)

# 优化筛选：仅保留P<0.05且中介比例>10%的显著路径
# 如果数据中没有中介比例字段，则仅按P值筛选
network_data <- mediation_summary %>%
  filter(
    !is.na(indirect_effect),
    indirect_effect_pval < 0.05  # 名义显著路径
  )

# 如果有中介比例字段，添加中介比例>10%的筛选条件
if (has_mediation_proportion) {
  network_data <- network_data %>%
    filter(
      !is.na(mediation_proportion_percent),
      abs(mediation_proportion_percent) > 10  # 中介比例>10%
    )
  cat("  → 已应用中介比例>10%的筛选条件\n")
}

# 按重要性排序：先按中介比例绝对值降序，再按P值升序
if (has_mediation_proportion) {
  network_data <- network_data %>%
    arrange(desc(abs(mediation_proportion_percent)), indirect_effect_pval)
} else {
  network_data <- network_data %>%
    arrange(indirect_effect_pval)
}

# 去重并限制数量（最多保留最重要的路径）
network_data <- network_data %>%
  distinct(exposure, mediator, outcome, .keep_all = TRUE) %>%
  head(15)  # 增加显示路径数，但通常会被中介比例筛选掉大部分

cat("  → 筛选了", nrow(network_data), "条关键路径（P < 0.05")
if (has_mediation_proportion) {
  cat(", 中介比例>10%")
}
cat("）用于网络图（已去重）\n")
if (nrow(network_data) > 0) {
  cat("    路径列表：\n")
  for (i in seq_len(nrow(network_data))) {
    cat(sprintf("    %d. %s → %s → %s (P = %.4f)\n", 
                i, network_data$exposure[i], network_data$mediator[i], 
                network_data$outcome[i], network_data$indirect_effect_pval[i]))
  }
} else {
  cat("    ⚠ 警告：没有找到显著路径，将尝试使用所有p<0.1的路径\n")
}
flush.console()

# 添加计算列（安全地处理缺失的列）
if (has_exp_to_med_beta && has_med_to_out_beta) {
  network_data <- network_data %>%
    mutate(
      # 节点标签
      exposure_label = exposure,
      mediator_label = mediator,
      outcome_label = case_when(
        outcome == "lung_cancer_overall" ~ "Overall\nLung Cancer",
        outcome == "lung_adenocarcinoma" ~ "Lung\nAdenocarcinoma",
        outcome == "squamous_cell_lung" ~ "Squamous Cell\nLung Cancer",
        TRUE ~ outcome
      ),
      # 效应方向
      exp_to_med_direction = ifelse(exp_to_med_beta > 0, "Positive", "Negative"),
      med_to_out_direction = ifelse(med_to_out_beta > 0, "Positive", "Negative"),
      indirect_direction = ifelse(indirect_effect > 0, "Positive", "Negative"),
      # 效应大小（绝对值）
      exp_to_med_abs = abs(exp_to_med_beta),
      med_to_out_abs = abs(med_to_out_beta),
      indirect_abs = abs(indirect_effect),
      # 显著性
      is_significant = indirect_effect_pval < 0.05,
      is_fdr_significant = if(has_significant_fdr) significant_fdr else rep(FALSE, n()),
      # 中介比例（用于标注）
      mediation_prop_percent = if(has_mediation_proportion) mediation_proportion_percent else rep(NA_real_, n()),
      # 格式化标签：中介比例和P值
      mediation_label = ifelse(
        has_mediation_proportion & !is.na(mediation_proportion_percent),
        sprintf("%.1f%%", mediation_proportion_percent),
        ""
      ),
      pval_label = sprintf("P=%.3f", indirect_effect_pval)
    )
} else {
  cat("⚠ 警告：缺少 exp_to_med_beta 或 med_to_out_beta 列，网络图将仅显示间接效应\n")
  network_data <- network_data %>%
    mutate(
      exposure_label = exposure,
      mediator_label = mediator,
      outcome_label = case_when(
        outcome == "lung_cancer_overall" ~ "Overall\nLung Cancer",
        outcome == "lung_adenocarcinoma" ~ "Lung\nAdenocarcinoma",
        outcome == "squamous_cell_lung" ~ "Squamous Cell\nLung Cancer",
        TRUE ~ outcome
      ),
      indirect_direction = ifelse(indirect_effect > 0, "Positive", "Negative"),
      indirect_abs = abs(indirect_effect),
      is_significant = indirect_effect_pval < 0.05,
      is_fdr_significant = if(has_significant_fdr) significant_fdr else rep(FALSE, n()),
      # 中介比例（用于标注）
      mediation_prop_percent = if(has_mediation_proportion) mediation_proportion_percent else rep(NA_real_, n()),
      # 格式化标签：中介比例和P值
      mediation_label = ifelse(
        has_mediation_proportion & !is.na(mediation_proportion_percent),
        sprintf("%.1f%%", mediation_proportion_percent),
        ""
      ),
      pval_label = sprintf("P=%.3f", indirect_effect_pval),
      # 为缺失的列设置默认值
      exp_to_med_beta = NA_real_,
      med_to_out_beta = NA_real_,
      exp_to_med_pval = NA_real_,
      med_to_out_pval = NA_real_,
      exp_to_med_direction = "Unknown",
      med_to_out_direction = "Unknown",
      exp_to_med_abs = 0,
      med_to_out_abs = 0
    )
}

if (nrow(network_data) == 0) {
  cat("⚠ 警告：没有显著路径用于网络图，使用所有路径（P < 0.1）\n")
  network_data <- mediation_summary %>%
    filter(!is.na(indirect_effect), indirect_effect_pval < 0.1)
  
  # 添加计算列（安全地处理缺失的列）
  if (has_exp_to_med_beta && has_med_to_out_beta) {
    network_data <- network_data %>%
      mutate(
        exposure_label = exposure,
        mediator_label = mediator,
        outcome_label = case_when(
          outcome == "lung_cancer_overall" ~ "Overall\nLung Cancer",
          outcome == "lung_adenocarcinoma" ~ "Lung\nAdenocarcinoma",
          outcome == "squamous_cell_lung" ~ "Squamous Cell\nLung Cancer",
          TRUE ~ outcome
        ),
        exp_to_med_direction = ifelse(exp_to_med_beta > 0, "Positive", "Negative"),
        med_to_out_direction = ifelse(med_to_out_beta > 0, "Positive", "Negative"),
        indirect_direction = ifelse(indirect_effect > 0, "Positive", "Negative"),
        exp_to_med_abs = abs(exp_to_med_beta),
        med_to_out_abs = abs(med_to_out_beta),
        indirect_abs = abs(indirect_effect),
        is_significant = indirect_effect_pval < 0.05,
        is_fdr_significant = if(has_significant_fdr) significant_fdr else rep(FALSE, n())
      )
  } else {
    network_data <- network_data %>%
      mutate(
        exposure_label = exposure,
        mediator_label = mediator,
        outcome_label = case_when(
          outcome == "lung_cancer_overall" ~ "Overall\nLung Cancer",
          outcome == "lung_adenocarcinoma" ~ "Lung\nAdenocarcinoma",
          outcome == "squamous_cell_lung" ~ "Squamous Cell\nLung Cancer",
          TRUE ~ outcome
        ),
        indirect_direction = ifelse(indirect_effect > 0, "Positive", "Negative"),
        indirect_abs = abs(indirect_effect),
        is_significant = indirect_effect_pval < 0.05,
        is_fdr_significant = if(has_significant_fdr) significant_fdr else rep(FALSE, n()),
        exp_to_med_beta = NA_real_,
        med_to_out_beta = NA_real_,
        exp_to_med_pval = NA_real_,
        med_to_out_pval = NA_real_,
        exp_to_med_direction = "Unknown",
        med_to_out_direction = "Unknown",
        exp_to_med_abs = 0,
        med_to_out_abs = 0
      )
  }
}

if (nrow(network_data) > 0) {
  cat("  → 开始构建网络数据...\n")
  flush.console()
  
  # 创建边数据框（两段路径：暴露->中介 和 中介->结局）
  # 只有在有数据的情况下才创建边
  if (has_exp_to_med_beta && has_med_to_out_beta && 
      !all(is.na(network_data$exp_to_med_beta)) && 
      !all(is.na(network_data$med_to_out_beta))) {
    # 创建完整路径信息映射（用于标注）
    pathway_info <- network_data %>%
      mutate(
        pathway_id = paste(exposure, mediator, outcome, sep = "_")
      ) %>%
      select(
        pathway_id,
        mediation_label, pval_label, mediation_prop_percent,
        indirect_effect_pval
      ) %>%
      distinct(pathway_id, .keep_all = TRUE)
    
    edges_exp_med <- network_data %>%
      mutate(pathway_id = paste(exposure, mediator, outcome, sep = "_")) %>%
      select(from = exposure, to = mediator, 
             weight = exp_to_med_beta, abs_weight = exp_to_med_abs,
             direction = exp_to_med_direction, pval = exp_to_med_pval,
             pathway_id) %>%
      left_join(pathway_info %>% select(pathway_id, mediation_label, pval_label, mediation_prop_percent, 
                                       path_pval = indirect_effect_pval),
                by = "pathway_id") %>%
      filter(!is.na(weight), !is.na(abs_weight)) %>%
      # 标记边类型，便于后续识别
      mutate(edge_type = "Exposure→Mediator") %>%
      # 基于 from 和 to 去重，保留权重绝对值最大的边
      arrange(desc(abs_weight)) %>%
      distinct(from, to, .keep_all = TRUE)
    
    cat("  → 已创建暴露->中介边：", nrow(edges_exp_med), "条（已去重）\n")
    flush.console()
    
    edges_med_out <- network_data %>%
      mutate(pathway_id = paste(exposure, mediator, outcome, sep = "_")) %>%
      select(from = mediator, to = outcome, 
             weight = med_to_out_beta, abs_weight = med_to_out_abs,
             direction = med_to_out_direction, pval = med_to_out_pval,
             pathway_id) %>%
      left_join(pathway_info %>% select(pathway_id, mediation_label, pval_label, mediation_prop_percent,
                                       path_pval = indirect_effect_pval),
                by = "pathway_id") %>%
      filter(!is.na(weight), !is.na(abs_weight)) %>%
      # 标记边类型，便于后续识别
      mutate(edge_type = "Mediator→Outcome") %>%
      # 基于 from 和 to 去重，保留权重绝对值最大的边
      arrange(desc(abs_weight)) %>%
      distinct(from, to, .keep_all = TRUE)
    
      cat("  → 已创建中介->结局边：", nrow(edges_med_out), "条（已去重）\n")
    flush.console()
    
    # 检查去重效果
    if (nrow(network_data) > 0) {
      total_exp_med_combinations <- nrow(network_data)
      unique_exp_med_edges <- nrow(edges_exp_med)
      total_med_out_combinations <- nrow(network_data)
      unique_med_out_edges <- nrow(edges_med_out)
      if (unique_exp_med_edges < total_exp_med_combinations) {
        cat(sprintf("    注：暴露->中介边去重：%d条路径合并为%d条唯一边\n", 
                    total_exp_med_combinations, unique_exp_med_edges))
      }
      if (unique_med_out_edges < total_med_out_combinations) {
        cat(sprintf("    注：中介->结局边去重：%d条路径合并为%d条唯一边\n", 
                    total_med_out_combinations, unique_med_out_edges))
      }
      flush.console()
    }
  } else {
    # 如果缺少分段路径数据，仅使用间接效应创建简化的网络
    cat("  ⚠ 警告：缺少分段路径数据，将创建基于间接效应的简化网络\n")
    edges_exp_med <- data.frame(
      from = character(0),
      to = character(0),
      weight = numeric(0),
      abs_weight = numeric(0),
      direction = character(0),
      pval = numeric(0),
      pathway_id = character(0),
      stringsAsFactors = FALSE
    )
    edges_med_out <- data.frame(
      from = character(0),
      to = character(0),
      weight = numeric(0),
      abs_weight = numeric(0),
      direction = character(0),
      pval = numeric(0),
      pathway_id = character(0),
      stringsAsFactors = FALSE
    )
    cat("  → 无法创建分段边，将使用简化网络结构\n")
    flush.console()
  }
  
  edges_all <- bind_rows(edges_exp_med, edges_med_out)
  
  if (nrow(edges_all) > 0) {
    # 最终去重：确保每条 from->to 边只出现一次（保留权重最大的）
    edges_before_dedup <- nrow(edges_all)
    edges_all <- edges_all %>%
      arrange(desc(abs_weight)) %>%
      distinct(from, to, .keep_all = TRUE) %>%
      mutate(
        edge_significant = if("pval" %in% names(.)) pval < 0.05 else FALSE,
        # Edge width proportional to absolute effect size (scale_edge_width_continuous handles scaling)
        edge_width = abs_weight
      )
    
    # 诊断输出：验证去重效果
    cat("  → 边去重检查：\n")
    cat(sprintf("    去重前：%d条边\n", edges_before_dedup))
    cat(sprintf("    去重后：%d条边\n", nrow(edges_all)))
    if (edges_before_dedup > nrow(edges_all)) {
      cat(sprintf("    移除重复边：%d条\n", edges_before_dedup - nrow(edges_all)))
    } else {
      cat("    无重复边\n")
    }
    flush.console()
    
    # 创建节点数据框（确保去重）
    all_nodes <- unique(c(edges_all$from, edges_all$to))
    cat("  → 从边数据中提取了", length(all_nodes), "个唯一节点\n")
    flush.console()
  } else {
    # 如果没有边，基于network_data创建节点（确保去重）
    cat("  ⚠ 警告：没有边数据，将基于路径数据创建节点列表\n")
    all_nodes <- unique(c(
      network_data$exposure,
      network_data$mediator,
      network_data$outcome
    ))
    cat("  → 基于路径数据找到", length(all_nodes), "个唯一节点\n")
    flush.console()
  }
  
  cat("  → 总节点数：", length(all_nodes), "个\n")
  flush.console()
  
  # 分类节点
  exposures <- unique(network_data$exposure)
  mediators <- unique(network_data$mediator)
  outcomes <- unique(network_data$outcome)
  
  cat("  → 暴露数：", length(exposures), "个，中介数：", length(mediators), "个，结局数：", length(outcomes), "个\n")
  flush.console()
  
  # 创建节点数据框（确保节点名称唯一）
  nodes_df <- data.frame(
    name = all_nodes,
    stringsAsFactors = FALSE
  ) %>%
    distinct(name, .keep_all = FALSE) %>%  # 确保节点名称唯一
    mutate(
      type = case_when(
        name %in% exposures ~ "Exposure",
        name %in% mediators ~ "Mediator",
        name %in% outcomes ~ "Outcome",
        TRUE ~ "Unknown"
      ),
      # 计算节点的度（连接数）
      degree = if(nrow(edges_all) > 0) {
        sapply(name, function(n) {
          sum(edges_all$from == n) + sum(edges_all$to == n)
        })
      } else {
        # 如果没有边，基于network_data中的路径频率计算度
        sapply(name, function(n) {
          sum(network_data$exposure == n) + 
          sum(network_data$mediator == n) + 
          sum(network_data$outcome == n)
        })
      },
      # 生成英文标签（确保每个节点都有标签，即使是未映射的也显示原始名称）
      label = case_when(
        name %in% exposures ~ get_english_label(name, exposure_label_mapping),
        name %in% mediators ~ get_english_label(name, mediator_label_mapping),
        name %in% outcomes ~ get_english_label(name, outcome_label_mapping),
        TRUE ~ as.character(name)  # 确保转换为字符，避免NA
      )
    ) %>%
    # 确保所有标签都是非空的字符串（处理可能的NA或空值）
    mutate(
      label = ifelse(is.na(label) | label == "" | trimws(label) == "", 
                     as.character(name), 
                     as.character(label))
    )
  
  # 检查并处理标签重复问题（如果同一标签对应多个节点）
  # 使用节点原始名称的后缀来区分，而不是数字编号
  duplicate_labels <- nodes_df %>%
    group_by(label) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 1)
  
  if (nrow(duplicate_labels) > 0) {
    cat("  ⚠ 警告：发现", nrow(duplicate_labels), "个重复标签，将使用节点名称区分\n")
    # 为重复标签添加唯一标识符
    nodes_df <- nodes_df %>%
      group_by(label) %>%
      mutate(
        # 如果标签重复，使用节点名称的最后部分作为区分符
        name_part = gsub(".*_", "", name),  # 提取下划线后的部分
        label_unique = ifelse(
          n() > 1,
          paste0(label, "\n(", name_part, ")"),
          label
        )
      ) %>%
      ungroup() %>%
      mutate(label = label_unique) %>%
      select(-label_unique, -name_part)
  }
  
  # 验证节点和标签的唯一性和完整性
  if (length(unique(nodes_df$name)) != nrow(nodes_df)) {
    cat("  ⚠ 警告：节点名称存在重复！\n")
  }
  if (length(unique(nodes_df$label)) != nrow(nodes_df)) {
    cat("  ⚠ 警告：节点标签存在重复！\n")
  }
  
  # 检查是否有空标签或NA标签
  missing_labels <- nodes_df %>% filter(is.na(label) | label == "" | trimws(label) == "")
  if (nrow(missing_labels) > 0) {
    cat("  ⚠ 警告：发现", nrow(missing_labels), "个节点缺少标签，将使用节点名称作为标签\n")
    nodes_df <- nodes_df %>%
      mutate(label = ifelse(is.na(label) | label == "" | trimws(label) == "", 
                           as.character(name), 
                           as.character(label)))
  }
  
  cat("  → 节点数据框包含", nrow(nodes_df), "个唯一节点\n")
  cat("  → 节点标签验证：", sum(!is.na(nodes_df$label) & nodes_df$label != ""), "个有效标签\n")
  cat("  → 节点类型分布：", 
      "暴露=", sum(nodes_df$type == "Exposure"), 
      "，中介=", sum(nodes_df$type == "Mediator"),
      "，结局=", sum(nodes_df$type == "Outcome"), "\n")
  flush.console()
  
  cat("  → 创建igraph对象...\n")
  flush.console()
  
  # 创建igraph对象
  if (nrow(edges_all) > 0) {
    net <- graph_from_data_frame(d = edges_all, vertices = nodes_df, directed = TRUE)
  } else {
    # 如果没有边，创建一个只有节点的空图
    net <- graph_from_data_frame(d = data.frame(from = character(0), to = character(0)), 
                                 vertices = nodes_df, directed = TRUE)
  }
  
  cat("  → 转换为tidygraph格式...\n")
  flush.console()
  
  # 转换为tidygraph格式（ggraph需要）
  # 为边添加路径标注信息（只在中介->结局边上显示完整路径信息）
  edges_all <- edges_all %>%
    mutate(
      # 判断是否为中介->结局边（这类边显示完整路径信息）
      is_med_to_out = if("edge_type" %in% names(.)) {
        edge_type == "Mediator→Outcome"
      } else {
        # 尝试从pathway_id推断：如果to是结局变量，则可能是中介->结局边
        to %in% outcomes
      },
      # 组合标注文本（仅在中介->结局边显示完整路径信息）
      edge_label = ifelse(
        is_med_to_out & !is.na(mediation_label) & mediation_label != "" & !is.na(pval_label),
        paste0(mediation_label, "\n", pval_label),
        ifelse(
          is_med_to_out & !is.na(mediation_label) & mediation_label != "",
          mediation_label,
          ""
        )
      ),
      # 路径显著性（用于颜色和线型）
      path_significant = ifelse(!is.na(path_pval), path_pval < 0.05, 
                               ifelse(!is.na(edge_significant), edge_significant, FALSE))
    )
  
  # 创建节点名称到索引的映射
  node_mapping <- data.frame(
    name = V(net)$name,
    idx = 1:vcount(net),
    stringsAsFactors = FALSE
  )
  
  # 为edges_all添加索引列以便后续合并
  edges_with_idx <- edges_all %>%
    left_join(node_mapping, by = c("from" = "name")) %>%
    rename(from_idx = idx) %>%
    left_join(node_mapping, by = c("to" = "name")) %>%
    rename(to_idx = idx)
  
  net_tidy <- as_tbl_graph(net) %>%
    activate(edges) %>%
    # 合并边的标注信息（使用整数索引进行匹配）
    mutate(edge_id = row_number()) %>%
    left_join(
      edges_with_idx %>% 
        select(from_idx, to_idx, edge_label, path_significant, mediation_prop_percent, path_pval),
      by = c("from" = "from_idx", "to" = "to_idx")
    ) %>%
    mutate(
      # 路径显著性：优先使用路径P值，否则使用边P值
      path_sig = ifelse(!is.na(path_significant), path_significant, 
                       ifelse(!is.na(edge_significant), edge_significant, FALSE)),
      # 优化线型：路径显著用实线，否则用虚线
      linetype = ifelse(path_sig, "solid", "dashed"),
      # 优化颜色：显著路径用红色，次要路径用灰色（基于路径显著性）
      # 优先使用路径显著性，如果没有则使用边显著性
      edge_color_category = case_when(
        !is.na(path_significant) & path_significant ~ "significant_path",
        !is.na(path_significant) & !path_significant ~ "non_significant_path",
        path_sig ~ "significant_path",
        TRUE ~ "non_significant_path"
      )
    ) %>%
    activate(nodes) %>%
    mutate(
      type = type,
      degree = degree,
      label = label
    )
  
  cat("  → 已完成网络对象创建\n")
  flush.console()
  
  # 使用改进的手动分层布局（三列：左侧暴露、中间中介、右侧结局）
  # 为每个节点创建坐标，确保节点顺序稳定
  node_coords <- data.frame(
    name = all_nodes,
    x = ifelse(all_nodes %in% exposures, 0, 
               ifelse(all_nodes %in% mediators, 1, 2)),
    stringsAsFactors = FALSE
  )
  
  # 按类型和名称排序，确保布局可重复
  node_coords <- node_coords[order(node_coords$x, node_coords$name), ]
  
  # 初始化y坐标列
  node_coords$y <- 0
  
  # ========================================================================
  # 改进的节点布局算法：基于连接度的智能排序和间距分配
  # ========================================================================
  # 为每个节点计算连接度（degree）以便智能排序
  node_degrees <- data.frame(
    name = V(net_tidy)$name,
    degree = degree(net_tidy),
    stringsAsFactors = FALSE
  )
  
  # 将连接度信息加入node_coords
  node_coords <- node_coords %>%
    left_join(node_degrees, by = "name")
  
  # 在每层内智能分配y坐标
  # 策略：按连接度排序，连接度高的节点靠近中心，减少边交叉
  for (layer in 0:2) {
    layer_indices <- node_coords$x == layer
    layer_data <- node_coords[layer_indices, ]
    n_layer <- nrow(layer_data)
    
    if (n_layer > 1) {
      # 【改进1】根据节点数量自适应调整间距
      # 使用更智能的间距计算，考虑整体布局平衡
      if (n_layer <= 2) {
        spacing <- 6.0  # 2个节点：大间距
      } else if (n_layer <= 3) {
        spacing <- 5.5  # 3个节点
      } else if (n_layer <= 4) {
        spacing <- 5.0  # 4个节点
      } else if (n_layer <= 5) {
        spacing <- 4.5  # 5个节点
      } else if (n_layer <= 7) {
        spacing <- 4.0  # 6-7个节点
      } else if (n_layer <= 10) {
        spacing <- 3.5  # 8-10个节点
      } else {
        spacing <- 3.0  # 11+个节点
      }
      
      # 【改进2】按连接度排序，高连接度节点靠中心
      # 这样可以减少边的交叉，使网络更清晰
      layer_data <- layer_data[order(-layer_data$degree, layer_data$name), ]
      
      # 【改进3】计算y坐标，采用居中对称布局
      total_height <- (n_layer - 1) * spacing
      
      # 使用对称分布，奇数个节点时中间节点在y=0
      # 偶数个节点时中间两个节点对称分布在y=0两侧
      y_values <- seq(total_height / 2, -total_height / 2, length.out = n_layer)
      
      # 【改进4】为中介层节点添加微调，避免完全对齐造成的视觉混乱
      # 只对中介层（layer=1）应用，且节点数>3时启用
      if (layer == 1 && n_layer > 3) {
        # 添加小幅度随机偏移（控制在spacing的10%以内）
        set.seed(42)  # 固定随机种子确保可重复
        jitter_amount <- spacing * 0.1
        y_jitter <- runif(n_layer, -jitter_amount, jitter_amount)
        y_values <- y_values + y_jitter
      }
      
      # 应用y坐标
      layer_data$y <- y_values
      node_coords[layer_indices, ] <- layer_data
      
    } else if (n_layer == 1) {
      # 单个节点居中
      node_coords$y[layer_indices] <- 0
    }
  }
  
  # 增加层间水平间距，使三列布局更清晰
  # 将x坐标从0,1,2扩展到更大的间距
  node_coords <- node_coords %>%
    mutate(x = x * 4.0)  # 将层间距从1增加到4，给节点和标签更多水平空间
  
  # 转换为矩阵格式（ggraph需要）
  # 按节点名称排序以确保匹配
  node_coords <- node_coords[order(node_coords$name), ]
  layout_matrix <- as.matrix(node_coords[, c("x", "y")])
  rownames(layout_matrix) <- node_coords$name
  
  cat("  → 开始绘制网络图（这可能需要一些时间）...\n")
  flush.console()
  
  # 创建网络图（使用分层布局算法）
  # 使用手动分层布局（三列布局）
  p_network <- ggraph(net_tidy, layout = layout_matrix) +
    # 边（箭头）- 使用弧形避免重叠，使用路径显著性决定颜色，包含标签
    geom_edge_arc(
      aes(width = abs_weight, 
          color = edge_color_category, 
          linetype = linetype,
          label = ifelse(edge_label != "" & !is.na(edge_label), edge_label, "")),
      arrow = arrow(length = unit(3, "mm"), type = "closed", ends = "last"),
      alpha = 0.8,      # 🔧 提高箭头不透明度（原0.7→0.8），更清晰 ⭐
      strength = 0.3,   # 弧度强度：0=直线，1=最大弧度，0.3是温和的弧度
      lineend = "round",
      angle_calc = "along",
      label_colour = "black",
      label_size = 3.2, # 🔧 增大P值字体（原3.0→3.2），更易读 ⭐
      label_alpha = 1.0, # 🔧 P值标签完全不透明，确保可见 ⭐
      label_dodge = unit(6, 'mm'),  # 🔧 v4.3增强：标签横向偏移增大（4mm→6mm），更远离节点 ⭐⭐⭐
      label_push = unit(3, 'mm')    # 🔧 v4.3增强：标签纵向推移增大（2mm→3mm），避免遮挡 ⭐⭐⭐
    ) +
    scale_edge_width_continuous(
      name = "Effect Size",  # 删除换行符，简化标签
      range = c(0.5, 3),
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        label.hjust = 0.5
      )
    ) +
    scale_edge_color_manual(
      name = "Significance",  # 简化标题
      values = c(
        "significant_path" = okabe_ito$vermillion,  # 显著路径用红色
        "non_significant_path" = "gray60"  # 次要路径用灰色
      ),
      labels = c(
        "significant_path" = "P < 0.05",  # 简化标签
        "non_significant_path" = "P ≥ 0.05"  # 简化标签
      ),
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(size = 2, linetype = "solid")
      )
    ) +
    scale_edge_linetype_manual(
      name = "Edge Type",  # 简化标题
      values = c("solid" = "solid", "dashed" = "dashed"),
      labels = c("solid" = "Significant", "dashed" = "Non-significant"),  # 简化标签
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(size = 2, color = "gray40")
      )
    ) +
    # 节点
    geom_node_point(
      aes(fill = type, size = degree),
      shape = 21,
      color = "white",
      stroke = 0.8
    ) +
    scale_fill_manual(
      name = "Variable Type",
      values = c(
        "Exposure" = okabe_ito$orange,
        "Mediator" = okabe_ito$sky_blue,
        "Outcome" = okabe_ito$green
      ),
      labels = c(
        "Exposure" = "Metabolic Exposure",
        "Mediator" = "Inflammatory Mediator",
        "Outcome" = "Lung Cancer"
      ),
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(size = 5, shape = 21)
      )
    ) +
    scale_size_continuous(
      range = c(5, 14),  # Increase node size for small datasets
      guide = "none"  # 移除Connectivity图例
    ) +
    # ====================================================================
    # 智能标签系统：基于节点位置和连接度的自适应标签定位
    # ====================================================================
    geom_node_label(  # 🆕 改为 label（带背景框，标签更完整清晰）
      aes(label = label),
      # 【改进1】根据节点总数自适应调整字体大小
      size = case_when(
        length(all_nodes) <= 10 ~ 4.2,   # 少量节点：大字体
        length(all_nodes) <= 15 ~ 3.8,   # 中等节点：中等字体
        length(all_nodes) <= 20 ~ 3.4,   # 较多节点：小字体
        TRUE ~ 3.0                       # 大量节点：更小字体
      ),
      repel = TRUE,  # 启用智能避让
      # 【改进2】增加标签偏移，让标签远离节点避免遮挡
      nudge_y = 0.6,  # 垂直偏移增加到0.6（原0.5）
      nudge_x = 0.4,  # 水平偏移增加到0.4（原0.3）
      # 【改进3】连接线设置优化
      segment.size = 0.4,       # 连接线适度变细，更精致
      segment.color = "gray50", # 连接线浅灰色，不抢眼
      segment.alpha = 0.6,      # 适度透明，避免视觉干扰
      # 【改进4】显著增强避让参数
      box.padding = 1.2,        # 标签框内边距增加到1.2（原1.0）
      point.padding = 0.8,      # 节点周围避让空间增加到0.8（原0.6）
      force = 20,               # 避让力显著增强到20（原15）⭐
      force_pull = 0.005,       # 减少回拉力（原0.01），允许标签移动更远
      # 【改进5】其他优化参数
      max.overlaps = Inf,       # 确保所有标签都显示
      min.segment.length = 0,   # 允许极短连接线
      direction = "both",       # 全方向避让
      xlim = c(-Inf, Inf),      # 取消x轴限制
      ylim = c(-Inf, Inf),      # 取消y轴限制
      # 🆕【改进6】label特有参数（带背景框，确保标签完整显示）
      fill = "white",                      # 白色背景
      alpha = 0.75,                        # 🔧 降低到0.75（原0.9），让箭头和P值可见 ⭐
      label.size = 0.25,                   # 边框线宽度
      label.padding = unit(0.25, "lines"), # 🔧 减小到0.25（原0.3），减少遮挡范围
      label.r = unit(0.15, "lines"),       # 圆角半径
      seed = 42                            # 固定随机种子，确保可重复
    ) +
    # Theme and labels
    labs(
      title = "Metabolic-Inflammatory-Lung Cancer Mediation Network",
      subtitle = sprintf("Showing %d nominally significant mediation paths (P < 0.05)", 
                        nrow(network_data)),
      caption = NULL  # Remove caption, legends are clear enough
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 5, unit = "mm")),
      plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(b = 10, unit = "mm")),
      plot.caption = element_text(size = 7, hjust = 0.5, color = "gray40", margin = margin(t = 10, unit = "mm")),
      legend.position = "bottom",  # 改为底部水平排列，减少水平空间占用
      legend.box = "horizontal",   # 水平排列图例
      legend.title = element_text(size = 8, face = "bold"),  # 从9减小到8
      legend.text = element_text(size = 7),  # 从8减小到7
      legend.spacing = unit(2, "mm"),  # 从4减小到2
      legend.box.spacing = unit(3, "mm"),  # 从5减小到3
      legend.key.height = unit(3, "mm"),  # 从4减小到3
      legend.key.width = unit(5, "mm"),  # 从6减小到5
      plot.margin = margin(20, 20, 25, 20, unit = "mm")  # 🆕 增加边距（原15,15,20,15）
    ) +
    coord_fixed(clip = "off")  # 🆕 关闭边界裁剪，确保标签完整显示
  
  cat("  → 网络图对象已创建，开始保存...\n")
  flush.console()
  
  # 保存图表（根据节点数量调整尺寸，增大以确保内容完整显示）
  n_nodes <- length(all_nodes)
  if (n_nodes <= 15) {
    # 少量节点时，使用更大的尺寸确保标签完整显示
    fig_width <- 240   # 🆕 从200增加到240（+20%），给标签更多空间
    fig_height <- 180  # 🆕 从150增加到180（+20%）
  } else {
    fig_width <- 240   # 🆕 从200增加到240（+20%）
    fig_height <- 200  # 🆕 从170增加到200（+18%）
  }
  
  cat("  → 保存图表到文件...\n")
  flush.console()
  
  save_sci_figure(
    p_network,
    "results/figures/step07_publication/Figure1_Mediation_Network",
    width_mm = fig_width,
    height_mm = fig_height
  )
  
  cat(sprintf("✓ 网络图已生成：%d 个节点，%d 条边，%d 条显著路径\n", 
              length(all_nodes), nrow(edges_all), nrow(network_data)))
  cat(sprintf("  图表尺寸：%.0f × %.0f mm\n\n", fig_width, fig_height))
  flush.console()
} else {
  cat("⚠ 警告：没有网络数据可用于生成网络图\n\n")
  flush.console()
}

# ============================================================================
# 图表2：亚型特异性对比图（鳞癌vs腺癌的中介效应对比）
# ============================================================================

cat("【图表2】生成亚型特异性对比图\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 筛选腺癌和鳞癌的数据
subtype_data <- mediation_summary %>%
  filter(outcome %in% c("lung_adenocarcinoma", "squamous_cell_lung")) %>%
  mutate(
    outcome_label = ifelse(outcome == "lung_adenocarcinoma", "Adenocarcinoma", "Squamous"),
    exposure_label_eng = sapply(exposure, function(x) get_english_label(x, exposure_label_mapping)),
    mediator_label_eng = sapply(mediator, function(x) get_english_label(x, mediator_label_mapping)),
    pathway_label = paste0(exposure_label_eng, " → ", mediator_label_eng)
  )

if (nrow(subtype_data) > 0) {
  # 创建亚型比较数据
  # 对于每个暴露-中介对，比较其在两个亚型中的效应
  unique_paths <- subtype_data %>%
    select(exposure, mediator) %>%
    distinct()
  
  comparison_data <- data.frame()
  
  for (i in seq_len(nrow(unique_paths))) {
    path <- unique_paths[i, ]
    
    adeno <- subtype_data %>%
      filter(exposure == path$exposure, mediator == path$mediator, 
             outcome == "lung_adenocarcinoma")
    
    squamous <- subtype_data %>%
      filter(exposure == path$exposure, mediator == path$mediator,
             outcome == "squamous_cell_lung")
    
    if (nrow(adeno) > 0 && nrow(squamous) > 0) {
      exp_label <- get_english_label(path$exposure, exposure_label_mapping)
      med_label <- get_english_label(path$mediator, mediator_label_mapping)
      comp_row <- data.frame(
        exposure = path$exposure,
        mediator = path$mediator,
        pathway_label = paste0(exp_label, " → ", med_label),
        adeno_indirect = adeno$indirect_effect[1],
        squamous_indirect = squamous$indirect_effect[1],
        adeno_pval = adeno$indirect_effect_pval[1],
        squamous_pval = squamous$indirect_effect_pval[1],
        adeno_fdr = if("fdr_pval_indirect" %in% names(adeno)) adeno$fdr_pval_indirect[1] else NA,
        squamous_fdr = if("fdr_pval_indirect" %in% names(squamous)) squamous$fdr_pval_indirect[1] else NA,
        stringsAsFactors = FALSE
      )
      
      comparison_data <- rbind(comparison_data, comp_row)
    }
  }
  
  if (nrow(comparison_data) > 0) {
    # 筛选显著路径（至少在其中一个亚型中显著）
    sig_comparison <- comparison_data %>%
      filter(adeno_pval < 0.05 | squamous_pval < 0.05) %>%
      arrange(desc(abs(adeno_indirect) + abs(squamous_indirect))) %>%
      head(20)  # 最多显示20条路径
    
    if (nrow(sig_comparison) > 0) {
      # 准备绘图数据（长格式）
      plot_data <- sig_comparison %>%
        mutate(
          pathway_label = factor(pathway_label, levels = rev(unique(pathway_label)))
        ) %>%
        pivot_longer(
          cols = c(adeno_indirect, squamous_indirect),
          names_to = "subtype",
          values_to = "indirect_effect"
        ) %>%
        mutate(
          subtype = ifelse(subtype == "adeno_indirect", "Adenocarcinoma", "Squamous"),
          pval = ifelse(subtype == "Adenocarcinoma", adeno_pval, squamous_pval),
          significant = pval < 0.05,
          # 判断是否为亚型特异性
          is_squamous_preferred = (subtype == "Squamous" & squamous_pval < 0.05 & 
                                   (is.na(adeno_pval) | adeno_pval >= 0.05)),
          is_adeno_preferred = (subtype == "Adenocarcinoma" & adeno_pval < 0.05 & 
                               (is.na(squamous_pval) | squamous_pval >= 0.05))
        )
      
      # 创建发散条形图
      p_subtype_comparison <- ggplot(plot_data, aes(x = indirect_effect, y = pathway_label, fill = subtype)) +
        geom_bar(
          stat = "identity",
          position = "dodge",
          alpha = 0.8,
          width = 0.7
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.21) +
        scale_fill_manual(
          name = "Lung Cancer Subtype",
          values = c("Adenocarcinoma" = okabe_ito$orange, "Squamous" = okabe_ito$sky_blue),
          guide = guide_legend(title.position = "top")
        ) +
        labs(
          title = "Figure 2. Subtype-Specific Mediation Effects",
          subtitle = sprintf("Comparison of indirect effects between adenocarcinoma and squamous cell carcinoma\n(%d pathways shown)", nrow(sig_comparison)),
          x = "Indirect Effect (β)",
          y = "Mediation Pathway",
          caption = "Bars represent indirect effects through inflammatory mediators. Positive values indicate risk effects, negative values indicate protective effects."
        ) +
        theme_sci() +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom"
        )
      
      # 保存图表
      save_sci_figure(
        p_subtype_comparison,
        "results/figures/step07_publication/Figure2_Subtype_Comparison",
        width_mm = 174,
        height_mm = max(120, nrow(sig_comparison) * 6)
      )
      
      cat(sprintf("✓ 亚型对比图已生成：%d 条路径\n\n", nrow(sig_comparison)))
    } else {
      cat("⚠ 警告：没有显著路径可用于亚型对比图\n\n")
    }
  } else {
    cat("⚠ 警告：没有可比较的亚型数据\n\n")
  }
} else {
  cat("⚠ 警告：没有亚型特异性数据\n\n")
}

# ============================================================================
# 图表3：效应大小分布图（所有路径的间接效应分布）
# ============================================================================

cat("【图表3】生成效应大小分布图\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 准备数据
effect_dist_data <- mediation_summary %>%
  filter(!is.na(indirect_effect)) %>%
  mutate(
    # 暴露类别分类
    exposure_category = case_when(
      exposure %in% c("BMI", "fasting_insulin", "fasting_glucose", "HbA1c") ~ "Glucose Metabolism",
      exposure %in% c("HDL_cholesterol", "LDL_cholesterol", "triglycerides", "ApoB", "ApoA1",
                      "HDL_large", "HDL_diameter", "HDL_very_large", "LDL_small",
                      "remnant_cholesterol", "ApoB_ApoA1_ratio") ~ "Lipid Metabolism",
      exposure %in% c("SBP", "DBP", "hypertension") ~ "Blood Pressure",
      exposure %in% c("smoking_initiation", "alcohol_drinks") ~ "Lifestyle",
      exposure %in% c("IGF1", "circulating_leptin", "vitamin_D") ~ "Hormones/Growth Factors",
      exposure %in% c("GGT") ~ "Liver Function",
      exposure %in% c("BCAA") ~ "Amino Acids",
      TRUE ~ "Other"
    ),
    # 中介类别
    mediator_label = case_when(
      mediator == "CRP" ~ "CRP",
      mediator == "IL6" ~ "IL-6",
      mediator == "IL6R" ~ "IL-6R",
      mediator == "TNFR1" ~ "TNF-α",
      mediator == "WBC" ~ "WBC",
      TRUE ~ mediator
    )
  )

if (nrow(effect_dist_data) > 0) {
  # 创建直方图/密度图的组合
  p_effect_distribution <- ggplot(effect_dist_data, aes(x = indirect_effect, fill = exposure_category)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      alpha = 0.7,
      color = "white",
      linewidth = 0.21,
      position = "identity"
    ) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.21) +
    scale_fill_manual(
      name = "Exposure Category",
      values = c(
        "Glucose Metabolism" = okabe_ito$orange,
        "Lipid Metabolism" = okabe_ito$sky_blue,
        "Blood Pressure" = okabe_ito$green,
        "Lifestyle" = okabe_ito$vermillion,
        "Hormones/Growth Factors" = okabe_ito$reddish_purple,
        "Liver Function" = okabe_ito$yellow,
        "Amino Acids" = okabe_ito$blue,
        "Other" = okabe_ito$gray
      ),
      guide = guide_legend(title.position = "top", ncol = 2)
    ) +
    labs(
      title = "Figure 3. Distribution of Indirect Effects Across All Mediation Pathways",
      subtitle = sprintf("Histogram and density plot of %d mediation pathways", nrow(effect_dist_data)),
      x = "Indirect Effect (β)",
      y = "Density",
      caption = "Color indicates exposure category. Dashed line at β = 0 represents no effect."
    ) +
    theme_sci() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # 保存图表
  save_sci_figure(
    p_effect_distribution,
    "results/figures/step07_publication/Figure3_Effect_Distribution",
    width_mm = 174,
    height_mm = 120
  )
  
  cat(sprintf("✓ 效应分布图已生成：%d 条路径\n\n", nrow(effect_dist_data)))
} else {
  cat("⚠ 警告：没有数据可用于效应分布图\n\n")
}

# ============================================================================
# 图表4：成功vs失败路径对比图（中介分析挑战的可视化）
# ============================================================================

cat("【图表4】生成成功vs失败路径对比图\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 尝试加载完整的结果以获取失败路径信息
if (file.exists("data/step07_all_mediation_results.RData")) {
  load("data/step07_all_mediation_results.RData")
  
  # 统计成功和失败的路径
  # 成功：至少完成了两个步骤（暴露->中介 和 中介->结局）
  # 失败：有错误信息或无法完成分析
  
  if (!is.null(mediation_results)) {
    # 从mediation_results中提取失败信息
    failure_reasons <- data.frame(
      pathway_id = character(),
      status = character(),
      reason = character(),
      stringsAsFactors = FALSE
    )
    
    for (key in names(mediation_results)) {
      result <- mediation_results[[key]]
      
      if (!is.null(result)) {
        # 判断状态
        if (!is.null(result$exp_to_med) && !is.null(result$med_to_out)) {
          if (!is.null(result$exp_to_out)) {
            status <- "Complete Success"
          } else if (isTRUE(result$partial_success)) {
            status <- "Partial Success"
          } else {
            status <- "Partial Failure"
          }
        } else {
          status <- "Failure"
        }
        
        # 提取失败原因（重新分类为更明确的类别）
        if (status %in% c("Failure", "Partial Failure")) {
          reason <- if (!is.na(result$error_message)) {
            # 重新分类失败原因
            error_msg <- result$error_message
            
            # 判断是哪个阶段的工具不足
            if (grepl("无法获取.*暴露.*工具变量|暴露.*SNP数量", error_msg, ignore.case = TRUE) ||
                (is.null(result$exp_to_med) && grepl("工具变量|IV|SNP", error_msg, ignore.case = TRUE))) {
              "Exposure-Mediator IVs Insufficient"
            } else if (grepl("无法获取.*中介.*工具变量|中介.*SNP数量|mediator.*IV|mediator.*SNP", error_msg, ignore.case = TRUE) ||
                       (is.null(result$med_to_out) && grepl("工具变量|IV|SNP", error_msg, ignore.case = TRUE))) {
              "Mediator-Outcome IVs Insufficient"
            } else if (grepl("无法提取结局数据|Outcome.*extract|outcome.*failed", error_msg, ignore.case = TRUE)) {
              "Mediator-Outcome IVs Insufficient"  # 结局数据提取失败通常也是中介-结局阶段的问题
            } else if (grepl("数据协调|harmoniz|palindromic|strand", error_msg, ignore.case = TRUE)) {
              "Harmonization Failed"
            } else if (grepl("SNP数量不足|insufficient.*SNP", error_msg, ignore.case = TRUE)) {
              # 根据上下文判断是哪个阶段
              if (is.null(result$exp_to_med)) {
                "Exposure-Mediator IVs Insufficient"
              } else {
                "Mediator-Outcome IVs Insufficient"
              }
            } else {
              "Other Technical Issues"
            }
          } else {
            "Other Technical Issues"
          }
          
          failure_reasons <- rbind(failure_reasons, data.frame(
            pathway_id = key,
            status = status,
            reason = reason,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
      # 如果有失败信息，创建失败原因统计（使用新的分类）
      if (nrow(failure_reasons) > 0) {
        # 重新映射失败原因标签为更友好的显示名称
        failure_reasons$reason_display <- case_when(
          failure_reasons$reason == "Exposure-Mediator IVs Insufficient" ~ "Exposure-Mediator\nIVs Insufficient",
          failure_reasons$reason == "Mediator-Outcome IVs Insufficient" ~ "Mediator-Outcome\nIVs Insufficient",
          failure_reasons$reason == "Harmonization Failed" ~ "Harmonization\nFailed",
          failure_reasons$reason == "Other Technical Issues" ~ "Other Technical\nIssues",
          TRUE ~ failure_reasons$reason
        )
        
        failure_summary <- failure_reasons %>%
          group_by(reason, reason_display) %>%
          summarise(count = n(), .groups = "drop") %>%
          arrange(desc(count)) %>%
          mutate(
            reason_label = factor(reason_display, levels = unique(reason_display[order(-count)])),
            percentage = 100 * count / sum(count)
          )
      
      # 计算成功路径数
      total_paths <- if (!is.null(all_pathways)) length(all_pathways) else 
                     if (exists("total_pathways")) total_pathways else nrow(mediation_summary) * 3
      
      success_count <- nrow(mediation_summary)
      failure_count <- total_paths - success_count
      
      # 创建成功vs失败对比数据
      success_failure_data <- data.frame(
        Status = c("Success", "Failure"),
        Count = c(success_count, failure_count),
        Percentage = c(100 * success_count / total_paths, 100 * failure_count / total_paths)
      )
      
      # 创建瀑布图风格的图表
      # 计算y轴最大值，确保标签有足够空间
      max_count <- max(success_failure_data$Count)
      y_max <- max_count * 1.25  # 增加25%的空间给标签
      
          # 面板A：总体成功率
          p_success_failure <- ggplot(success_failure_data, aes(x = Status, y = Count, fill = Status)) +
            geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
            geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
                      vjust = -0.2, size = 3.5, fontface = "bold", lineheight = 0.9) +
            scale_fill_manual(
              name = "Analysis Status",
              values = c("Success" = okabe_ito$green, "Failure" = okabe_ito$vermillion),
              guide = "none"
            ) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, y_max)) +
            labs(
              title = "Panel A: Overall Analysis Status",
              subtitle = sprintf("Total pathways: %d", total_paths),
              x = "Analysis Status",
              y = "Number of Pathways",
              caption = NULL
            ) +
            theme_sci() +
            theme(
              plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 9)
            )
      
      # 如果有失败原因数据，创建失败原因分解图
      if (nrow(failure_summary) > 0) {
        # 计算y轴最大值，确保标签有足够空间
        max_count_failure <- max(failure_summary$count)
        y_max_failure <- max_count_failure * 1.25  # 增加25%的空间给标签
        
          # 面板B：失败原因细分（基于失败路径）
          p_failure_reasons <- ggplot(failure_summary, aes(x = reason_label, y = count, fill = reason_label)) +
            geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
            geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, percentage)),
                      vjust = -0.2, size = 3, fontface = "bold", lineheight = 0.9) +
            scale_fill_manual(
              name = "Failure Reason",
              values = c(
                "Exposure-Mediator\nIVs Insufficient" = okabe_ito$vermillion,
                "Mediator-Outcome\nIVs Insufficient" = okabe_ito$orange,
                "Harmonization\nFailed" = okabe_ito$yellow,
                "Other Technical\nIssues" = okabe_ito$gray
              ),
              guide = "none"
            ) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, y_max_failure)) +
            labs(
              title = sprintf("Panel B: Failure Reasons Breakdown (n=%d)", failure_count),
              subtitle = "Detailed breakdown of reasons for analysis failure",
              x = "Failure Reason",
              y = "Number of Failed Pathways",
              caption = NULL
            ) +
            theme_sci() +
            theme(
              plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 9),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
            )
        
        # 使用cowplot组合两个图
        if (requireNamespace("cowplot", quietly = TRUE)) {
          library(cowplot)
          # 修改子图标题，移除主标题，只保留副标题
          p_success_failure_sub <- p_success_failure + 
            theme(plot.margin = margin(5, 5, 5, 5, unit = "mm")) +
            labs(title = NULL)
          
          p_failure_reasons_sub <- p_failure_reasons + 
            theme(plot.margin = margin(5, 5, 5, 5, unit = "mm"))
          
          # 组合两个图
          p_combined <- plot_grid(
            p_success_failure_sub,
            p_failure_reasons_sub,
            ncol = 2,
            align = "h",
            rel_widths = c(1, 1.5)
          )
          
          # 添加统一的主标题
          title_combined <- ggdraw() + 
            draw_text("Figure 4. Mediation Analysis Success and Failure Reasons",
                     x = 0.5, hjust = 0.5, y = 0.5, size = 12, fontface = "bold")
          
          p_final <- plot_grid(
            title_combined,
            p_combined,
            ncol = 1,
            rel_heights = c(0.12, 1)  # 标题区域高度
          )
          
          save_sci_figure(
            p_final,
            "results/figures/step07_publication/Figure4_Success_Failure_Analysis",
            width_mm = 174,
            height_mm = 120
          )
        } else {
          # 如果cowplot不可用，只保存成功失败对比图
          save_sci_figure(
            p_success_failure,
            "results/figures/step07_publication/Figure4_Success_Failure_Analysis",
            width_mm = 174,
            height_mm = 100
          )
          
          save_sci_figure(
            p_failure_reasons,
            "results/figures/step07_publication/Figure4_Failure_Reasons",
            width_mm = 174,
            height_mm = 100
          )
        }
      } else {
        save_sci_figure(
          p_success_failure,
          "results/figures/step07_publication/Figure4_Success_Failure_Analysis",
          width_mm = 174,
          height_mm = 100
        )
      }
      
      cat(sprintf("✓ 成功vs失败对比图已生成\n"))
      cat(sprintf("  - 成功路径: %d (%.1f%%)\n", success_count, 100 * success_count / total_paths))
      cat(sprintf("  - 失败路径: %d (%.1f%%)\n", failure_count, 100 * failure_count / total_paths))
      if (nrow(failure_summary) > 0) {
        cat(sprintf("  - 失败原因分类: %d 类\n", nrow(failure_summary)))
      }
      cat("\n")
    } else {
      cat("⚠ 警告：无法提取失败路径信息，仅生成成功路径统计\n")
      
      # 仅基于mediation_summary创建统计
      total_analyzed <- nrow(mediation_summary)
      if (exists("total_pathways")) {
        p_success <- ggplot(data.frame(Status = "Success", Count = total_analyzed), 
                           aes(x = Status, y = Count)) +
          geom_bar(stat = "identity", fill = okabe_ito$green, alpha = 0.8) +
          labs(title = "Mediation Analysis: Successfully Analyzed Pathways",
               x = "Status", y = "Count") +
          theme_sci()
        
        save_sci_figure(
          p_success,
          "results/figures/step07_publication/Figure4_Success_Failure_Analysis",
          width_mm = 174,
          height_mm = 100
        )
      }
      cat("\n")
    }
  } else {
    cat("⚠ 警告：无法加载详细结果，跳过成功vs失败对比图\n\n")
  }
} else {
  cat("⚠ 警告：无法找到完整结果文件，基于summary数据创建简化版本\n")
  
  # 基于mediation_summary创建基本统计
  total_analyzed <- nrow(mediation_summary)
  cat(sprintf("  - 成功分析路径: %d\n\n", total_analyzed))
}

# ============================================================================
# 总结
# ============================================================================

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤7c完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【生成的图表】\n")
cat("  图表1: 中介路径网络图\n")
cat("    - results/figures/step07_publication/Figure1_Mediation_Network.png/pdf\n")
cat("  图表2: 亚型特异性对比图\n")
cat("    - results/figures/step07_publication/Figure2_Subtype_Comparison.png/pdf\n")
cat("  图表3: 效应大小分布图\n")
cat("    - results/figures/step07_publication/Figure3_Effect_Distribution.png/pdf\n")
cat("  图表4: 成功vs失败路径对比图\n")
cat("    - results/figures/step07_publication/Figure4_Success_Failure_Analysis.png/pdf\n\n")

cat("【图表特点】\n")
cat("  ✓ 符合SCI 10分期刊标准（174mm双栏宽度）\n")
cat("  ✓ 使用Okabe-Ito色盲友好调色板\n")
cat("  ✓ 高分辨率输出（600 DPI PNG + 矢量PDF）\n")
cat("  ✓ Arial字体，标准字号（9-11 pt）\n")
cat("  ✓ 完整的图例和标注\n\n")

cat("【使用建议】\n")
cat("  1. 检查生成的PNG文件（预览）\n")
cat("  2. 使用PDF文件用于论文投稿（矢量格式）\n")
cat("  3. 根据期刊要求调整尺寸或字体\n")
cat("  4. 所有图表已按照SCI标准格式化\n\n")

