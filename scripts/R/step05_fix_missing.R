############################################################################
# Step05 修复脚本：补做失败的MR分析
# 用途：针对缺失的暴露-结局对重新进行分析
############################################################################

cat("Step05 修复：补做失败的MR分析\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 1. 加载必要的包
required_packages <- c("TwoSampleMR", "dplyr", "openxlsx", "ggplot2", "gridExtra")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 创建必要的输出目录
dir.create("results/figures/supplementary_sensitivity", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/univariable_mr", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# 2. 检查是否有缺失分析数据
missing_file <- "data/step05_missing_analyses.RData"
if (!file.exists(missing_file)) {
  cat("未找到缺失分析数据文件，先运行诊断脚本...\n")
  source("step05_diagnose_missing.R")
  
  if (!file.exists(missing_file)) {
    stop("错误：诊断脚本未生成缺失分析数据文件")
  }
}

load(missing_file)  # 加载 missing_analyses 和 missing_df

if (length(missing_analyses) == 0) {
  cat("✓ 没有缺失的分析，无需修复！\n")
  quit(save = "no")
}

cat(sprintf("发现 %d 个缺失的分析，开始修复...\n\n", length(missing_analyses)))

# 3. 加载已有的结果（如果存在）
results_file <- "data/step05_all_results.RData"
all_mr_results <- list()
extraction_log <- data.frame(
  exposure = character(),
  category = character(),
  n_snps = integer(),
  strategy = character(),
  stringsAsFactors = FALSE
)

if (file.exists(results_file)) {
  load(results_file)
  cat(sprintf("✓ 已加载已有结果：%d 个分析\n", length(all_mr_results)))
} else {
  cat("⚠ 未找到已有结果文件，将创建新文件\n")
}

# 4. 加载必要的数据文件
instruments_file <- "data/step04_all_instruments.RData"
outcome_file <- "results/data/outcome_data_list.RData"

if (!file.exists(instruments_file) || !file.exists(outcome_file)) {
  stop("错误：找不到必要的数据文件\n")
}

load(instruments_file)
load(outcome_file)

cat(sprintf("✓ 已加载工具变量：%d 个暴露\n", length(all_instruments)))
cat(sprintf("✓ 已加载结局数据：%d 个结局\n\n", length(outcome_data_list)))

# 5. 定义名称映射（与Step5保持一致）
exposure_name_mapping <- list(
  circulating_leptin = "leptin",
  vitamin_D = "vitamin_d",
  HbA1c = "hba1c",
  ApoB = "apoB",
  ApoA1 = "apoA1",
  IGF1 = "igf1",
  ApoB_ApoA1_ratio = "apoB_apoA1_ratio",
  HDL_diameter = "hdl_diameter",
  HDL_large = "large_hdl",
  remnant_cholesterol = "remnant_chol",
  LDL_small = "small_ldl",
  BCAA = "bcaa",
  HDL_very_large = "very_large_hdl",
  BMI = "bmi",
  HDL_cholesterol = "hdl_chol",
  LDL_cholesterol = "ldl_chol",
  smoking_initiation = "smoking",
  alcohol_drinks = "alcohol",
  fasting_glucose = "fasting_glucose",
  fasting_insulin = "fasting_insulin",
  SBP = "sbp",
  DBP = "dbp",
  hypertension = "hypertension",
  triglycerides = "triglycerides",
  GGT = "ggt",
  CRP = "crp",
  WBC = "wbc",
  IL6 = "il6",
  IL6R = "il6r",
  TNFR1 = "tnfr1"
)

outcome_name_mapping <- list(
  lung_cancer_overall = "Lung_cancer_overall",
  lung_cancer_adenocarcinoma = "Lung_cancer_adenocarcinoma",
  lung_cancer_squamous = "Lung_cancer_squamous"
)

# 定义暴露分类
metabolic_exposures <- c("circulating_leptin", "vitamin_D", "HbA1c", "ApoB", "ApoA1", 
                        "IGF1", "ApoB_ApoA1_ratio", "HDL_diameter", "HDL_large", 
                        "remnant_cholesterol", "LDL_small", "BCAA", "HDL_very_large", 
                        "BMI", "HDL_cholesterol", "LDL_cholesterol", "fasting_glucose", 
                        "fasting_insulin", "triglycerides", "GGT", "SBP", "DBP", 
                        "hypertension", "smoking_initiation", "alcohol_drinks")
inflammatory_exposures <- c("CRP", "WBC", "IL6", "IL6R", "TNFR1")

# 6. 加载Step5中的必要函数
# 创建一个临时环境来加载Step5脚本中的函数，而不执行整个脚本
step5_env <- new.env()

# 读取Step5脚本并执行必要的部分（只加载函数定义，不执行分析）
source_lines <- readLines("step05_单变量MR分析_完整修复版.R")

# 提取函数定义部分
# 需要包含 outcome_gwas_id_mapping 和 outcome_name_to_id（在# 6之前）
# 以及 perform_mr_analysis_robust（在# 6中）
mapping_start <- grep("^# 定义结局GWAS ID映射", source_lines)
if (length(mapping_start) == 0) {
  # 如果找不到注释，尝试查找 outcome_gwas_id_mapping 定义
  mapping_start <- grep("outcome_gwas_id_mapping <- list", source_lines)
}
# 确保从 outcome_gwas_id_mapping 定义之前开始，以包含完整定义
if (length(mapping_start) > 0 && mapping_start > 1) {
  # 往前查找几行以包含前面的注释（如果有）
  mapping_start <- max(1, mapping_start - 2)
}
func_start <- grep("^# 6\\. 定义健壮的MR分析函数", source_lines)
func_end <- grep("^# 7\\. 定义暴露因子分类", source_lines) - 1

# 如果找到 mapping_start，从那里开始提取；否则从 func_start 开始
extract_start <- if (length(mapping_start) > 0 && mapping_start < func_start) {
  mapping_start
} else if (length(func_start) > 0) {
  func_start
} else {
  NULL
}

# 标志变量，跟踪是否成功加载函数
functions_loaded <- FALSE

if (!is.null(extract_start) && length(func_end) > 0 && func_end > extract_start) {
  function_section <- source_lines[extract_start:func_end]
  # 移除章节标题注释（以 "# " 或 "#\t" 开头），但保留代码中的注释
  # 保留所有包含代码的行和行内注释
  function_section <- function_section[!grepl("^\\s*#\\s+[0-9]", function_section)]
  
  # 在临时环境中执行函数定义
  tryCatch({
    eval(parse(text = paste(function_section, collapse = "\n")), envir = step5_env)
    
    # 检查函数是否存在
    if (exists("perform_mr_analysis_robust", envir = step5_env) && 
        exists("outcome_name_to_id", envir = step5_env)) {
      # 将函数复制到当前环境
      perform_mr_analysis_robust <- step5_env$perform_mr_analysis_robust
      outcome_name_to_id <- step5_env$outcome_name_to_id
      if (exists("outcome_gwas_id_mapping", envir = step5_env)) {
        outcome_gwas_id_mapping <- step5_env$outcome_gwas_id_mapping
      }
      functions_loaded <- TRUE
      cat("✓ 已加载MR分析函数和 outcome_name_to_id\n")
    } else {
      cat("⚠ 函数定义不完整，将使用简化版本\n")
    }
  }, error = function(e) {
    cat(sprintf("⚠ 加载函数时出错: %s，将使用简化版本\n", conditionMessage(e)))
  })
}

if (!functions_loaded) {
  # 如果找不到，使用简化版本
  cat("⚠ 未找到完整函数，使用简化版本\n")
  
  # 定义outcome_name_to_id函数（与主文件保持一致）
  outcome_gwas_id_mapping <- list(
    lung_cancer_overall = "ebi-a-GCST004748",    # Lung cancer
    lung_adenocarcinoma = "ieu-a-984",           # Lung adenocarcinoma  
    squamous_cell_lung = "ieu-a-989"             # Squamous cell lung cancer
  )
  
  outcome_name_to_id <- function(outcome_name) {
    # 尝试直接匹配
    if (outcome_name %in% names(outcome_gwas_id_mapping)) {
      return(outcome_gwas_id_mapping[[outcome_name]])
    }
    
    # 尝试通过名称匹配
    name_lower <- tolower(outcome_name)
    if (grepl("lung.*cancer.*overall|overall.*lung.*cancer", name_lower, ignore.case = TRUE)) {
      return("ebi-a-GCST004748")
    } else if (grepl("adenocarcinoma", name_lower, ignore.case = TRUE)) {
      return("ieu-a-984")
    } else if (grepl("squamous", name_lower, ignore.case = TRUE)) {
      return("ieu-a-989")
    }
    
    return(NULL)
  }
  
  # 定义简化版的perform_mr_analysis_robust函数
  perform_mr_analysis_robust <- function(exposure_data, outcome_data, exposure_name, outcome_name) {
    tryCatch({
      # 确保数据格式正确
      if (!"SNP" %in% names(exposure_data)) {
        if ("rsid" %in% names(exposure_data)) {
          exposure_data$SNP <- exposure_data$rsid
        } else {
          return(NULL)
        }
      }
      
      if (!"SNP" %in% names(outcome_data)) {
        if ("rsid" %in% names(outcome_data)) {
          outcome_data$SNP <- outcome_data$rsid
        } else {
          return(NULL)
        }
      }
      
      exposure_snps <- unique(exposure_data$SNP)
      outcome_snp_count <- length(unique(outcome_data$SNP))
      
      # 如果结局数据SNP太少，尝试从数据库提取
      if (outcome_snp_count < 50) {
        outcome_id <- outcome_name_to_id(outcome_name)
        if (!is.null(outcome_id) && length(exposure_snps) > 0) {
          tryCatch({
            outcome_data <- TwoSampleMR::extract_outcome_data(
              snps = exposure_snps,
              outcomes = outcome_id
            )
          }, error = function(e) {
            cat(sprintf("    提取失败: %s\n", conditionMessage(e)))
          })
        }
      }
      
      # 协调数据
      harmonized <- harmonise_data(exposure_data, outcome_data)
      if (is.null(harmonized) || nrow(harmonized) == 0) {
        return(NULL)
      }
      
      # MR分析
      n_harmonized <- nrow(harmonized)
      if (n_harmonized >= 3) {
        mr_results <- mr(harmonized)
        heterogeneity <- mr_heterogeneity(harmonized)
        pleiotropy <- mr_pleiotropy_test(harmonized)
        single_snp <- mr_singlesnp(harmonized)
        loo <- mr_leaveoneout(harmonized)
      } else if (n_harmonized == 2) {
        mr_results <- mr(harmonized, method_list = c("mr_ivw", "mr_wald_ratio"))
        heterogeneity <- NULL
        pleiotropy <- NULL
        single_snp <- mr_singlesnp(harmonized)
        loo <- NULL
      } else if (n_harmonized == 1) {
        mr_results <- mr(harmonized, method_list = "mr_wald_ratio")
        heterogeneity <- NULL
        pleiotropy <- NULL
        single_snp <- mr_singlesnp(harmonized)
        loo <- NULL
      } else {
        return(NULL)
      }
      
      return(list(
        exposure = exposure_name,
        outcome = outcome_name,
        harmonized_data = harmonized,
        mr_results = mr_results,
        heterogeneity = heterogeneity,
        pleiotropy = pleiotropy,
        single_snp = single_snp,
        loo = loo,
        n_snps = nrow(exposure_data),
        n_harmonized = n_harmonized
      ))
    }, error = function(e) {
      cat(sprintf("    ✗ 错误: %s\n", e$message))
      return(NULL)
    })
  }
}

# 7. 执行修复分析
cat("【开始修复分析】\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

success_count <- 0
for (i in seq_len(nrow(missing_df))) {
  exposure_name <- missing_df$exposure[i]
  outcome_name <- missing_df$outcome[i]
  key <- paste(exposure_name, outcome_name, sep = "_")
  
  cat(sprintf("[%d/%d] %s -> %s\n", i, nrow(missing_df), exposure_name, outcome_name))
  
  # 获取映射后的名称
  mapped_exposure_name <- if (exposure_name %in% names(exposure_name_mapping)) {
    exposure_name_mapping[[exposure_name]]
  } else {
    exposure_name
  }
  
  mapped_outcome_name <- if (outcome_name %in% names(outcome_name_mapping)) {
    outcome_name_mapping[[outcome_name]]
  } else {
    outcome_name
  }
  
  # 获取数据
  if (!mapped_exposure_name %in% names(all_instruments)) {
    cat(sprintf("  ⚠ 跳过：未找到工具变量数据\n"))
    next
  }
  
  if (!mapped_outcome_name %in% names(outcome_data_list)) {
    cat(sprintf("  ⚠ 跳过：未找到结局数据\n"))
    next
  }
  
  exposure_data <- all_instruments[[mapped_exposure_name]]
  outcome_data <- outcome_data_list[[mapped_outcome_name]]
  
  if (is.null(exposure_data) || nrow(exposure_data) == 0 ||
      is.null(outcome_data) || nrow(outcome_data) == 0) {
    cat(sprintf("  ⚠ 跳过：数据为空\n"))
    next
  }
  
  # 执行分析
  result <- perform_mr_analysis_robust(
    exposure_data = exposure_data,
    outcome_data = outcome_data,
    exposure_name = exposure_name,
    outcome_name = outcome_name
  )
  
  if (!is.null(result)) {
    all_mr_results[[key]] <- result
    success_count <- success_count + 1
    
    # 保存单个结果
    category <- ifelse(exposure_name %in% metabolic_exposures, "Metabolic", "Inflammatory")
    save(result, file = paste0("results/univariable_mr/", category, "_", key, ".RData"))
    
    cat(sprintf("  ✓ 分析成功 (%d SNPs协调)\n", result$n_harmonized))
    
    # 生成图表
    tryCatch({
      cat("    生成图表...\n")
      
      # 散点图
      scatter_plots <- mr_scatter_plot(result$mr_results, result$harmonized_data)
      if (is.list(scatter_plots) && length(scatter_plots) > 0) {
        p_scatter <- if (length(scatter_plots) == 1) {
          scatter_plots[[1]]
        } else {
          if (require("gridExtra", quietly = TRUE)) {
            gridExtra::arrangeGrob(grobs = scatter_plots, ncol = 1)
          } else {
            scatter_plots[[1]]
          }
        }
        
        ggsave(
          filename = file.path("results/figures/supplementary_sensitivity", 
                             paste0("scatter_", key, ".png")),
          plot = p_scatter,
          width = 8, height = 6, dpi = 300
        )
        cat("      ✓ 散点图已保存\n")
      }
      
      # 漏斗图
      funnel_plots <- mr_funnel_plot(result$single_snp)
      if (is.list(funnel_plots) && length(funnel_plots) > 0) {
        p_funnel <- if (length(funnel_plots) == 1) {
          funnel_plots[[1]]
        } else {
          if (require("gridExtra", quietly = TRUE)) {
            gridExtra::arrangeGrob(grobs = funnel_plots, ncol = 1)
          } else {
            funnel_plots[[1]]
          }
        }
        
        ggsave(
          filename = file.path("results/figures/supplementary_sensitivity", 
                             paste0("funnel_", key, ".png")),
          plot = p_funnel,
          width = 8, height = 6, dpi = 300
        )
        cat("      ✓ 漏斗图已保存\n")
      }
      
      # 留一法图
      if (!is.null(result$loo)) {
        loo_plots <- mr_leaveoneout_plot(result$loo)
        if (is.list(loo_plots) && length(loo_plots) > 0) {
          n_snps <- nrow(result$loo)
          plot_height <- max(6, min(30, n_snps * 0.25))
          
          p_loo <- if (length(loo_plots) == 1) {
            loo_plots[[1]]
          } else {
            if (require("gridExtra", quietly = TRUE)) {
              gridExtra::arrangeGrob(grobs = loo_plots, ncol = 1)
            } else {
              loo_plots[[1]]
            }
          }
          
          ggsave(
            filename = file.path("results/figures/supplementary_sensitivity", 
                               paste0("loo_", key, ".png")),
            plot = p_loo,
            width = 10, height = plot_height, dpi = 300, limitsize = FALSE
          )
          cat("      ✓ 留一法图已保存\n")
        }
      }
      
      cat("    ✓ 图表生成完成\n")
    }, error = function(e) {
      cat(sprintf("    ⚠ 图表生成失败: %s\n", e$message))
    })
  } else {
    cat(sprintf("  ✗ 分析失败\n"))
  }
}

# 8. 保存更新后的结果
save(all_mr_results, extraction_log, file = results_file)
cat(sprintf("\n✓ 已保存更新后的结果到: %s\n", results_file))

# 9. 重新生成汇总表
cat("\n【重新生成汇总表】\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 重新生成汇总表（从已有结果中）
results_summary <- data.frame()

for (key in names(all_mr_results)) {
  result <- all_mr_results[[key]]
  
  if (!is.null(result) && !is.null(result$mr_results)) {
    ivw <- result$mr_results[result$mr_results$method == "Inverse variance weighted", ]
    
    if (nrow(ivw) > 0) {
      or <- exp(ivw$b)
      or_lci <- exp(ivw$b - 1.96 * ivw$se)
      or_uci <- exp(ivw$b + 1.96 * ivw$se)
      
      het_p <- ifelse(!is.null(result$heterogeneity) && nrow(result$heterogeneity) > 0,
                     result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", "Q_pval"],
                     NA)
      
      pleo_p <- ifelse(!is.null(result$pleiotropy) && nrow(result$pleiotropy) > 0,
                      result$pleiotropy$pval, NA)
      
      exposure_name <- result$exposure
      category <- ifelse(exposure_name %in% metabolic_exposures, "Metabolic", "Inflammatory")
      
      extraction_info <- extraction_log[extraction_log$exposure == exposure_name, ]
      
      summary_row <- data.frame(
        category = category,
        exposure = result$exposure,
        outcome = result$outcome,
        n_snps = result$n_snps,
        n_harmonized = result$n_harmonized,
        extraction_strategy = ifelse(nrow(extraction_info) > 0,
                                    extraction_info$strategy[1], "未知"),
        beta = ivw$b,
        se = ivw$se,
        pval = ivw$pval,
        or = or,
        or_lci = or_lci,
        or_uci = or_uci,
        or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
        heterogeneity_p = het_p,
        pleiotropy_p = pleo_p,
        stringsAsFactors = FALSE
      )
      
      results_summary <- rbind(results_summary, summary_row)
    }
  }
}

# FDR校正
if (nrow(results_summary) > 0) {
  results_summary$fdr_pval <- p.adjust(results_summary$pval, method = "fdr")
  results_summary$significant_fdr <- results_summary$fdr_pval < 0.05
  results_summary$significant_nominal <- results_summary$pval < 0.05
}

# 保存汇总表
write.xlsx(results_summary, "results/tables/step05_mr_results_summary.xlsx", row.names = FALSE)
write.csv(results_summary, "results/tables/step05_mr_results_summary.csv", row.names = FALSE)
write.csv(extraction_log, "results/tables/step05_extraction_log.csv", row.names = FALSE)

cat("✓ 汇总表已更新\n")
cat("  - results/tables/step05_mr_results_summary.xlsx\n")
cat("  - results/tables/step05_mr_results_summary.csv\n\n")

# 10. 生成策略分布图（如果有extraction_log数据）
if (nrow(extraction_log) > 0 && "strategy" %in% names(extraction_log)) {
  cat("\n【生成策略分布图】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  strategy_map <- c(
    "严格" = "Strict",
    "中等" = "Moderate", 
    "宽松" = "Relaxed",
    "极宽松" = "Very Relaxed",
    "超级宽松" = "Extremely Relaxed",
    "Standard" = "Standard"
  )
  
  log_data <- extraction_log
  log_data$strategy_en <- ifelse(log_data$strategy %in% names(strategy_map),
                                strategy_map[log_data$strategy],
                                log_data$strategy)
  
  strategy_colors_en <- c(
    "Strict"   = "#0072B2",
    "Moderate" = "#E69F00",
    "Relaxed"  = "#009E73",
    "Very Relaxed" = "#CC79A7",
    "Extremely Relaxed" = "#D55E00",
    "Standard" = "#999999"
  )
  
  tryCatch({
    p_strategy <- ggplot(log_data, aes(x = strategy_en, fill = strategy_en)) +
      geom_bar(alpha = 1, color = "black", linewidth = 0.3) +
      scale_fill_manual(
        values = strategy_colors_en,
        name = "Extraction\nStrategy"
      ) +
      labs(
        title = "Distribution of Instrumental Variable Extraction Strategies",
        x = "Extraction Strategy",
        y = "Number of Exposures"
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        plot.margin = margin(10, 15, 10, 10)
      )
    
    pdf_file <- "results/figures/strategy_comparison_plot.pdf"
    png_file <- "results/figures/strategy_comparison_plot.png"
    
    ggsave(pdf_file, plot = p_strategy, width = 8, height = 5, dpi = 600)
    ggsave(png_file, plot = p_strategy, width = 8, height = 5, dpi = 600)
    
    cat("✓ 策略分布图已生成:\n")
    cat("  - PDF矢量图:", pdf_file, "\n")
    cat("  - PNG位图:", png_file, "\n")
  }, error = function(e) {
    cat("⚠ 图表保存出错:", e$message, "\n")
  })
}

cat("\n【修复完成统计】\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("修复尝试:       %d 个\n", nrow(missing_df)))
cat(sprintf("修复成功:       %d 个\n", success_count))
cat(sprintf("修复成功率:     %.1f%%\n\n", 100 * success_count / max(1, nrow(missing_df))))

cat("【保存的文件】\n")
cat("主要结果:\n")
cat("  - data/step05_all_results.RData\n")
cat("汇总表格:\n")
cat("  - results/tables/step05_mr_results_summary.xlsx\n")
cat("  - results/tables/step05_mr_results_summary.csv\n")
cat("  - results/tables/step05_extraction_log.csv\n")
cat("图表文件:\n")
cat("  - results/figures/supplementary_sensitivity/*.png (散点图、漏斗图、留一法图)\n")
if (nrow(extraction_log) > 0 && "strategy" %in% names(extraction_log)) {
  cat("  - results/figures/strategy_comparison_plot.pdf\n")
  cat("  - results/figures/strategy_comparison_plot.png\n")
}
cat("\n修复完成！\n")

