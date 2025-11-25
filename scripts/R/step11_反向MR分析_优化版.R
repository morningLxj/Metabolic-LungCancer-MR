############################################################################
# 步骤11：反向孟德尔随机化（Reverse MR）完整分析 - 优化版
# 目的：检验肺癌对代谢/炎症性状的因果效应，排除反向因果偏倚
# 优化：与第5步风格统一，生成完整图表和论文级别表格
############################################################################

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤11：反向孟德尔随机化（Reverse MR）分析 - 优化版\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ============================================================================
# 步骤1：环境准备与包加载
# ============================================================================

cat("【步骤1】环境准备与包加载\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 加载必要的R包
required_packages <- c("TwoSampleMR", "dplyr", "readr", "ggplot2", 
                       "openxlsx", "gridExtra", "tidyr", "stringr", 
                       "RColorBrewer", "pheatmap")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装包: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# 尝试加载MR-PRESSO
if (!require("MRPRESSO", character.only = TRUE, quietly = TRUE)) {
  tryCatch({
    if (!require("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org")
    }
    library(devtools)
    install_github("rondolab/MR-PRESSO")
    library(MRPRESSO)
  }, error = function(e) {
    cat("⚠ 警告：MR-PRESSO包未安装，将跳过PRESSO分析\n")
  })
}

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(openxlsx)
  library(gridExtra)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
})

cat("✓ 所有必需的R包已加载\n\n")

# 创建输出目录结构（与第5步一致）
dirs <- c(
  "data",
  "results/reverse_mr",
  "results/tables", 
  "results/tables/paper_tables",
  "results/figures/reverse_mr",
  "results/figures/reverse_mr/sensitivity",
  "results/figures/reverse_mr/main_figures"
)

for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

cat("✓ 输出目录已创建\n\n")

# ============================================================================
# 步骤2：定义完整的GWAS数据集信息
# ============================================================================

cat("【步骤2】定义完整的GWAS数据集信息\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 肺癌数据集（作为暴露）
lung_cancer_datasets <- list(
  lung_cancer_overall = "ebi-a-GCST004748",
  lung_adenocarcinoma = "ieu-a-984",
  squamous_cell_lung = "ieu-a-989"
)

# 代谢指标数据集（作为结局）
metabolic_traits <- list(
  circulating_leptin = "ebi-a-GCST90007316",
  vitamin_D = "ebi-a-GCST90000618",
  HbA1c = "ebi-a-GCST90014006",
  ApoB = "ebi-a-GCST90025952",
  ApoA1 = "ebi-a-GCST90025955",
  IGF1 = "ebi-a-GCST90025989",
  ApoB_ApoA1_ratio = "ebi-a-GCST90092810",
  HDL_diameter = "ebi-a-GCST90092828",
  HDL_large = "ebi-a-GCST90092851",
  remnant_cholesterol = "ebi-a-GCST90092943",
  LDL_small = "ebi-a-GCST90092963",
  BCAA = "ebi-a-GCST90092984",
  HDL_very_large = "ebi-a-GCST90093011",
  BMI = "ieu-b-40",
  HDL_cholesterol = "ieu-b-109",
  LDL_cholesterol = "ieu-b-110",
  smoking_initiation = "ieu-b-4877",
  alcohol_drinks = "ieu-b-73",
  fasting_glucose = "ebi-a-GCST90002232",
  fasting_insulin = "ebi-a-GCST90002238",
  SBP = "ieu-b-38",
  DBP = "ieu-b-39",
  hypertension = "ieu-b-5144",
  triglycerides = "ieu-b-111",
  GGT = "ebi-a-GCST90025966"
)

# 炎症标志物数据集（作为结局）
inflammatory_traits <- list(
  CRP = "ebi-a-GCST90029070",
  WBC = "ieu-b-30",
  IL6 = "ebi-a-GCST90012005",
  IL6R = "ebi-a-GCST90012025",
  TNFR1 = "ebi-a-GCST90012015"
)

# 合并所有结局
all_outcomes <- c(metabolic_traits, inflammatory_traits)

# 定义结局分类
outcome_category <- c(
  rep("Metabolic", length(metabolic_traits)),
  rep("Inflammatory", length(inflammatory_traits))
)
names(outcome_category) <- names(all_outcomes)

cat(sprintf("✓ 已定义 %d 个肺癌暴露数据集\n", length(lung_cancer_datasets)))
cat(sprintf("✓ 已定义 %d 个代谢性状结局\n", length(metabolic_traits)))
cat(sprintf("✓ 已定义 %d 个炎症标志物结局\n", length(inflammatory_traits)))
cat(sprintf("✓ 总计 %d 个反向MR分析配对\n\n", 
           length(lung_cancer_datasets) * length(all_outcomes)))

# ============================================================================
# 步骤3：定义增强的工具变量提取函数（多策略）
# ============================================================================

cat("【步骤3】定义增强的工具变量提取函数\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 多策略工具变量提取函数（与第5步一致）
extract_instruments_robust <- function(exposure_id, exposure_name, min_snps = 3) {
  strategies <- list(
    list(name = "严格", p = 5e-8, r2 = 0.001, kb = 10000),
    list(name = "中等", p = 5e-7, r2 = 0.001, kb = 10000),
    list(name = "宽松", p = 5e-6, r2 = 0.01, kb = 5000),
    list(name = "极宽松", p = 5e-5, r2 = 0.05, kb = 5000)
  )
  
  for (strategy in strategies) {
    tryCatch({
      instruments <- extract_instruments(
        outcomes = exposure_id,
        p1 = strategy$p,
        clump = TRUE,
        r2 = strategy$r2,
        kb = strategy$kb
      )
      
      if (!is.null(instruments) && nrow(instruments) >= min_snps) {
        instruments$extraction_strategy <- strategy$name
        return(instruments)
      }
    }, error = function(e) {
      # 静默失败，尝试下一个策略
    })
  }
  
  return(NULL)
}

cat("✓ 工具变量提取函数已定义（多策略）\n\n")

# ============================================================================
# 步骤4：辅助函数 - 计算工具变量强度（与第5步一致）
# ============================================================================

calculate_instrument_strength <- function(harmonized_data) {
  if (is.null(harmonized_data) || nrow(harmonized_data) == 0) {
    return(list(f_statistic = NA, r_squared = NA, mean_f = NA, min_f = NA))
  }
  
  # 使用eaf（effect allele frequency）作为MAF的近似
  harmonized_data$maf_approx <- ifelse(is.na(harmonized_data$eaf.exposure), 
                                      0.25,
                                      ifelse(harmonized_data$eaf.exposure > 0.5,
                                             1 - harmonized_data$eaf.exposure,
                                             harmonized_data$eaf.exposure))
  
  # 计算每个SNP的R²贡献：R² = 2 * beta² * MAF * (1-MAF)
  harmonized_data$r2_contribution <- 2 * (harmonized_data$beta.exposure^2) * 
    harmonized_data$maf_approx * (1 - harmonized_data$maf_approx)
  
  total_r2 <- sum(harmonized_data$r2_contribution, na.rm = TRUE)
  
  # F统计量 ≈ beta²/se²（每个SNP）
  harmonized_data$f_per_snp <- (harmonized_data$beta.exposure / harmonized_data$se.exposure)^2
  
  mean_f <- mean(harmonized_data$f_per_snp, na.rm = TRUE)
  min_f <- min(harmonized_data$f_per_snp, na.rm = TRUE)
  
  return(list(
    f_statistic = mean_f,
    r_squared = total_r2,
    mean_f = mean_f,
    min_f = min_f
  ))
}

# ============================================================================
# 步骤5：定义健壮的反向MR分析函数（增强版）
# ============================================================================

cat("【步骤5】定义健壮的反向MR分析函数（增强版）\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

perform_reverse_mr <- function(exposure_id, exposure_name, 
                               outcome_id, outcome_name,
                               pair_id = NULL,
                               verbose = FALSE) {
  
  # 初始化结果（与第5步结构一致）
  result <- list(
    pair_id = pair_id,
    exposure_id = exposure_id,
    exposure_name = exposure_name,
    outcome_id = outcome_id,
    outcome_name = outcome_name,
    status = "Failed",
    error_message = NA,
    extraction_strategy = NA,
    
    # SNP信息
    n_snps_extracted = 0,
    n_snps_harmonized = 0,
    
    # 工具变量强度
    iv_strength = NULL,
    
    # MR结果
    mr_results = NULL,
    harmonized_data = NULL,
    single_snp = NULL,
    loo = NULL,
    heterogeneity = NULL,
    pleiotropy = NULL,
    presso = NULL
  )
  
  tryCatch({
    # 步骤1: 提取工具变量
    exposure_instruments <- extract_instruments_robust(
      exposure_id = exposure_id,
      exposure_name = exposure_name,
      min_snps = 3
    )
    
    if (is.null(exposure_instruments)) {
      result$error_message <- "工具变量提取失败"
      return(result)
    }
    
    result$n_snps_extracted <- nrow(exposure_instruments)
    result$extraction_strategy <- exposure_instruments$extraction_strategy[1]
    
    # 步骤2: 提取结局数据
    outcome_data <- extract_outcome_data(
      snps = exposure_instruments$SNP,
      outcomes = outcome_id
    )
    
    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      result$error_message <- "结局数据提取失败"
      return(result)
    }
    
    # 步骤3: 数据协调
    harmonised_data <- harmonise_data(
      exposure_dat = exposure_instruments,
      outcome_dat = outcome_data,
      action = 2
    )
    
    # 过滤掉不保留的SNP
    harmonised_data <- harmonised_data[harmonised_data$mr_keep == TRUE, ]
    
    if (nrow(harmonised_data) < 3) {
      result$error_message <- sprintf("协调后SNP不足 (n=%d)", nrow(harmonised_data))
      return(result)
    }
    
    result$n_snps_harmonized <- nrow(harmonised_data)
    result$harmonized_data <- harmonised_data
    
    # 步骤4: 计算工具变量强度
    result$iv_strength <- calculate_instrument_strength(harmonised_data)
    
    # 步骤5: 执行MR分析（使用所有可用方法）
    n_harmonized <- nrow(harmonised_data)
    if (n_harmonized >= 3) {
      mr_methods <- c("mr_ivw", "mr_egger_regression", 
                      "mr_weighted_median", "mr_weighted_mode",
                      "mr_simple_mode")
    } else if (n_harmonized == 2) {
      mr_methods <- c("mr_ivw", "mr_wald_ratio")
    } else {
      mr_methods <- "mr_wald_ratio"
    }
    
    mr_results <- mr(harmonised_data, method_list = mr_methods)
    result$mr_results <- mr_results
    
    if (is.null(mr_results) || nrow(mr_results) == 0) {
      result$error_message <- "MR分析失败"
      return(result)
    }
    
    # 步骤6: 敏感性分析
    # 异质性检验
    if (n_harmonized >= 3) {
      het_test <- mr_heterogeneity(harmonised_data, method_list = c("mr_ivw", "mr_egger_regression"))
      result$heterogeneity <- het_test
    }
    
    # 多效性检验
    if (n_harmonized >= 4) {
      pleio_test <- mr_pleiotropy_test(harmonised_data)
      result$pleiotropy <- pleio_test
    }
    
    # Single SNP分析
    single_snp_result <- mr_singlesnp(harmonised_data)
    result$single_snp <- single_snp_result
    
    # Leave-one-out分析
    if (n_harmonized >= 4) {
      loo_result <- mr_leaveoneout(harmonised_data)
      result$loo <- loo_result
    }
    
    # MR-PRESSO分析（如果包可用且SNP数量足够）
    if (require("MRPRESSO", quietly = TRUE)) {
      tryCatch({
        if (n_harmonized >= 10) {  # PRESSO需要至少10个SNP
          presso_data <- data.frame(
            BetaOutcome = harmonised_data$beta.outcome,
            BetaExposure = harmonised_data$beta.exposure,
            SdOutcome = harmonised_data$se.outcome,
            SdExposure = harmonised_data$se.exposure
          )
          
          presso_result <- MRPRESSO::mr_presso(
            BetaOutcome = "BetaOutcome",
            BetaExposure = "BetaExposure",
            SdOutcome = "SdOutcome",
            SdExposure = "SdExposure",
            OUTLIERtest = TRUE,
            DISTORTIONtest = TRUE,
            data = presso_data,
            NbDistribution = 1000,
            SignifThreshold = 0.05
          )
          result$presso <- presso_result
        }
      }, error = function(e) {
        # 静默失败，不影响主分析
      })
    }
    
    result$status <- "Success"
    return(result)
    
  }, error = function(e) {
    result$error_message <- substr(as.character(e$message), 1, 150)
    return(result)
  })
}

cat("✓ 反向MR分析函数已定义（增强版，包含所有敏感性分析）\n\n")

# ============================================================================
# 步骤6: 执行所有反向MR分析（或加载已有结果）
# ============================================================================

cat("【步骤6】执行反向MR分析\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# 检查是否存在已保存的完整结果
results_file <- "data/step11_reverse_mr_complete_results.RData"
load_existing <- FALSE

if (file.exists(results_file)) {
  cat("⚠ 发现已保存的分析结果文件！\n")
  cat(sprintf("   文件路径: %s\n", results_file))
  
  # 检查文件修改时间
  file_info <- file.info(results_file)
  file_time <- file_info$mtime
  cat(sprintf("   保存时间: %s\n", format(file_time, "%Y-%m-%d %H:%M:%S")))
  
  cat("\n请选择操作:\n")
  cat("  1. 加载已有结果，跳过分析，直接生成图表（推荐）\n")
  cat("  2. 重新进行分析（会覆盖已有结果）\n")
  cat("\n")
  
  # 尝试交互式输入，如果不可用则默认加载
  if (interactive()) {
    user_choice <- readline("请输入选项 (1 或 2，默认为1): ")
    if (is.na(user_choice) || trimws(user_choice) == "" || trimws(user_choice) == "1") {
      load_existing <- TRUE
    }
  } else {
    # 非交互式环境，默认加载已有结果
    cat("⚠ 非交互式环境，自动加载已有结果...\n")
    load_existing <- TRUE
  }
  
  if (load_existing) {
    cat("\n正在加载已保存的分析结果...\n")
    load(results_file)
    cat("✓ 分析结果已加载\n")
    
    # 统计加载的数据
    cat("\n【加载的数据统计】\n")
    cat(sprintf("  - 总配对数量: %d\n", length(all_results)))
    success_count <- sum(sapply(all_results, function(x) x$status == "Success"))
    failed_count <- length(all_results) - success_count
    cat(sprintf("  - 成功: %d\n", success_count))
    cat(sprintf("  - 失败: %d\n", failed_count))
    cat(sprintf("  - 分析日志: %d 条记录\n", nrow(analysis_log)))
    cat("\n")
    
    # 跳过分析步骤，直接进入步骤7
    cat("✓ 跳过分析步骤，直接进入结果整理和图表生成...\n\n")
  } else {
    cat("✓ 将重新进行分析...\n\n")
    load_existing <- FALSE
  }
}

# 如果需要重新分析
if (!load_existing) {
  total_pairs <- length(lung_cancer_datasets) * length(all_outcomes)
  cat(sprintf("开始分析 %d 对反向MR配对...\n", total_pairs))
  cat("这可能需要较长时间，请耐心等待...\n\n")
  
  # 初始化结果列表和计数器
  all_results <- list()
  success_count <- 0
  failed_count <- 0
  analysis_counter <- 0

  # 创建分析日志
  analysis_log <- data.frame(
    exposure = character(),
    outcome = character(),
    category = character(),
    extraction_strategy = character(),
    n_snps = integer(),
    n_harmonized = integer(),
    status = character(),
    stringsAsFactors = FALSE
  )

  # 执行分析
  start_time <- Sys.time()

  for (exposure_name in names(lung_cancer_datasets)) {
    exposure_id <- lung_cancer_datasets[[exposure_name]]
    
    cat(sprintf("\n【%s】\n", exposure_name))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    
    for (outcome_name in names(all_outcomes)) {
      outcome_id <- all_outcomes[[outcome_name]]
      analysis_counter <- analysis_counter + 1
      
      cat(sprintf("  [%d/%d] %s -> %s ", 
                 analysis_counter, total_pairs, 
                 substr(exposure_name, 1, 20), 
                 substr(outcome_name, 1, 20)))
      
      # 执行分析
      result <- perform_reverse_mr(
        exposure_id = exposure_id,
        exposure_name = exposure_name,
        outcome_id = outcome_id,
        outcome_name = outcome_name,
        pair_id = sprintf("RM%03d", analysis_counter),
        verbose = FALSE
      )
      
      # 生成唯一键
      key <- paste(exposure_name, outcome_name, sep = "_")
      
      # 保存结果
      all_results[[key]] <- result
      
      # 更新日志
      category <- outcome_category[[outcome_name]]
      analysis_log <- rbind(analysis_log, data.frame(
        exposure = exposure_name,
        outcome = outcome_name,
        category = category,
        extraction_strategy = ifelse(is.na(result$extraction_strategy), 
                                     "N/A", result$extraction_strategy),
        n_snps = result$n_snps_extracted,
        n_harmonized = result$n_snps_harmonized,
        status = result$status,
        stringsAsFactors = FALSE
      ))
      
      # 输出状态并生成图表
      if (result$status == "Success") {
        success_count <- success_count + 1
        
        # 提取IVW结果用于显示
        ivw_result <- result$mr_results[result$mr_results$method == "Inverse variance weighted", ]
        if (nrow(ivw_result) > 0) {
          ivw_or <- exp(ivw_result$b[1])
          ivw_pval <- ivw_result$pval[1]
          cat(sprintf("✓ OR=%.3f, P=%.2e\n", ivw_or, ivw_pval))
        } else {
          cat("✓ 分析成功\n")
        }
        
        # 保存单个结果
        save(result, file = paste0("results/reverse_mr/", key, ".RData"))
        
        # 生成图表（与第5步一致）
        tryCatch({
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
              filename = file.path("results/figures/reverse_mr/sensitivity", 
                                 paste0("scatter_", key, ".png")),
              plot = p_scatter,
              width = 8, height = 6, dpi = 300
            )
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
              filename = file.path("results/figures/reverse_mr/sensitivity", 
                                 paste0("funnel_", key, ".png")),
              plot = p_funnel,
              width = 8, height = 6, dpi = 300
            )
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
                filename = file.path("results/figures/reverse_mr/sensitivity", 
                                   paste0("loo_", key, ".png")),
                plot = p_loo,
                width = 10, height = plot_height, dpi = 300, limitsize = FALSE
              )
            }
          }
        }, error = function(e) {
          # 静默失败，不影响主分析
        })
      } else {
        failed_count <- failed_count + 1
        cat(sprintf("✗ %s\n", result$error_message))
      }
      
      # API限制
      Sys.sleep(0.3)
    }
  }
  
  end_time <- Sys.time()
  time_elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("分析完成! 总用时: %.1f 分钟\n", as.numeric(time_elapsed)))
  cat(sprintf("✓ 成功: %d/%d (%.1f%%)\n", 
             success_count, total_pairs, 100 * success_count / total_pairs))
  cat(sprintf("✗ 失败: %d/%d (%.1f%%)\n\n", 
             failed_count, total_pairs, 100 * failed_count / total_pairs))
  
  # 保存完整结果
  save(all_results, analysis_log, lung_cancer_datasets, all_outcomes,
       file = "data/step11_reverse_mr_complete_results.RData")
  cat("✓ 已保存完整结果: data/step11_reverse_mr_complete_results.RData\n\n")
  
  # 保存分析日志
  write.csv(analysis_log, "results/tables/step11_analysis_log.csv", row.names = FALSE)
  cat("✓ 已保存分析日志: results/tables/step11_analysis_log.csv\n\n")
}

# ============================================================================
# 步骤7：整理和导出结果
# ============================================================================

cat("【步骤7】整理和导出结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# 提取成功的结果
successful_results <- Filter(function(x) x$status == "Success", all_results)

if (length(successful_results) > 0) {
  
  # 转换为数据框（与第5步格式一致）
  results_summary <- data.frame()
  
  for (key in names(successful_results)) {
    result <- successful_results[[key]]
    
    if (!is.null(result) && !is.null(result$mr_results)) {
      ivw <- result$mr_results[result$mr_results$method == "Inverse variance weighted", ]
      
      if (nrow(ivw) > 0) {
        or <- exp(ivw$b[1])
        or_lci <- exp(ivw$b[1] - 1.96 * ivw$se[1])
        or_uci <- exp(ivw$b[1] + 1.96 * ivw$se[1])
        
        het_p <- ifelse(!is.null(result$heterogeneity) && nrow(result$heterogeneity) > 0,
                       result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", "Q_pval"][1],
                       NA)
        
        pleo_p <- ifelse(!is.null(result$pleiotropy) && nrow(result$pleiotropy) > 0,
                        result$pleiotropy$pval[1], NA)
        
        category <- outcome_category[[result$outcome_name]]
        
        summary_row <- data.frame(
          category = category,
          exposure = result$exposure_name,
          outcome = result$outcome_name,
          n_snps = result$n_snps_extracted,
          n_harmonized = result$n_snps_harmonized,
          extraction_strategy = ifelse(is.na(result$extraction_strategy),
                                      "未知", result$extraction_strategy),
          beta = ivw$b[1],
          se = ivw$se[1],
          pval = ivw$pval[1],
          or = or,
          or_lci = or_lci,
          or_uci = or_uci,
          or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
          heterogeneity_p = het_p,
          pleiotropy_p = pleo_p,
      stringsAsFactors = FALSE
    )
        
        # 添加工具变量强度
        if (!is.null(result$iv_strength)) {
          summary_row$mean_f_statistic <- result$iv_strength$mean_f
          summary_row$r_squared <- result$iv_strength$r_squared
        } else {
          summary_row$mean_f_statistic <- NA
          summary_row$r_squared <- NA
        }
        
        results_summary <- rbind(results_summary, summary_row)
      }
    }
  }
  
  # FDR校正
  if (nrow(results_summary) > 0) {
    results_summary$fdr_pval <- p.adjust(results_summary$pval, method = "fdr")
    results_summary$significant_fdr <- results_summary$fdr_pval < 0.05
    results_summary$significant_nominal <- results_summary$pval < 0.05
  
  # 按P值排序
    results_summary <- results_summary %>% arrange(pval)
  
    cat(sprintf("✓ 成功分析: %d 条\n", nrow(results_summary)))
  cat(sprintf("  - 名义显著 (P<0.05): %d (%.1f%%)\n", 
               sum(results_summary$significant_nominal),
               100 * sum(results_summary$significant_nominal) / nrow(results_summary)))
    cat(sprintf("  - FDR显著 (FDR<0.05): %d (%.1f%%)\n\n", 
               sum(results_summary$significant_fdr),
               100 * sum(results_summary$significant_fdr) / nrow(results_summary)))
    
    # 保存汇总表
    write.xlsx(results_summary, "results/tables/step11_mr_results_summary.xlsx", row.names = FALSE)
    write.csv(results_summary, "results/tables/step11_mr_results_summary.csv", row.names = FALSE)
    write.csv(analysis_log, "results/tables/step11_extraction_log.csv", row.names = FALSE)
    
    cat("✓ 汇总表已保存:\n")
    cat("  - results/tables/step11_mr_results_summary.xlsx\n")
    cat("  - results/tables/step11_mr_results_summary.csv\n")
    cat("  - results/tables/step11_extraction_log.csv\n\n")
  }
  
} else {
  cat("⚠ 所有分析均失败，无法生成结果摘要\n\n")
  results_summary <- data.frame()
}

# ============================================================================
# 步骤8：生成增强版论文表格（与第5步一致）
# ============================================================================

if (length(successful_results) > 0 && nrow(results_summary) > 0) {
  cat("【步骤8】生成增强版论文表格...\n")
  cat(paste(rep("-", 80), collapse = ""), "\n\n")
  
  # 安全提取函数（避免"更换参数长度为零"错误）
  safe_extract <- function(x, idx = 1, default = NA) {
    tryCatch({
      if (is.null(x)) {
        return(default)
      }
      # 处理长度为0的向量
      if (length(x) == 0) {
        return(default)
      }
      # 处理索引超出范围
      if (is.null(idx) || idx < 1 || idx > length(x)) {
        return(default)
      }
      # 提取值
      val <- x[idx]
      # 处理NULL
      if (is.null(val)) {
        return(default)
      }
      # 处理长度为0
      if (length(val) == 0) {
        return(default)
      }
      # 如果值是向量，取第一个元素
      if (length(val) > 1) {
        val <- val[1]
      }
      # 处理NA（但保留NA作为有效值，只在明确需要时替换）
      # 这里不自动替换NA，因为NA可能是有效的结果
      return(val)
    }, error = function(e) {
      return(default)
    })
  }
  
  # Table S1: 完整MR结果（所有方法）
  table_s1_data <- list()
  # Table 2增强版: 主MR结果（包含所有方法列）
  table2_enhanced_data <- list()
  # Table 5: 敏感性分析汇总
  table5_data <- list()
  # Table S2: 敏感性分析详情
  table_s2_data <- list()
  # Table S3: Leave-one-out详情
  table_s3_data <- list()
  # Table S6: 工具变量强度
  table_s6_data <- list()
  
  for (key in names(successful_results)) {
    result <- successful_results[[key]]
    if (is.null(result) || is.null(result$mr_results)) next
    
    exposure <- result$exposure_name
    outcome <- result$outcome_name
    category <- outcome_category[[outcome]]
    
    extraction_info <- analysis_log[analysis_log$exposure == exposure & 
                                   analysis_log$outcome == outcome, ]
    strategy <- ifelse(nrow(extraction_info) > 0, extraction_info$strategy[1], "Unknown")
    
    # 提取所有MR方法结果
    if (nrow(result$mr_results) > 0) {
      for (i in seq_len(nrow(result$mr_results))) {
        method_row <- result$mr_results[i, ]
        method_name <- method_row$method
        
        beta <- safe_extract(method_row$b, 1, NA)
        se <- safe_extract(method_row$se, 1, NA)
        pval <- safe_extract(method_row$pval, 1, NA)
        nsnp_val <- safe_extract(method_row$nsnp, 1, result$n_snps_harmonized)
        nsnp <- ifelse(is.na(nsnp_val), result$n_snps_harmonized, nsnp_val)
        
        # 跳过无效结果
        if (is.na(beta) || is.na(se)) next
        
        or <- exp(beta)
        or_lci <- exp(beta - 1.96 * se)
        or_uci <- exp(beta + 1.96 * se)
        
        table_s1_data[[length(table_s1_data) + 1]] <- data.frame(
          category = category,
          exposure = exposure,
          outcome = outcome,
          method = method_name,
          n_snps = nsnp,
          beta = beta,
          se = se,
          pval = pval,
          or = or,
          or_lci = or_lci,
          or_uci = or_uci,
          or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
          extraction_strategy = ifelse(is.na(strategy), "Unknown", strategy),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # IVW结果作为主结果
    ivw_result <- result$mr_results[result$mr_results$method == "Inverse variance weighted", ]
    if (nrow(ivw_result) == 0) {
      ivw_result <- result$mr_results[1, ]
    }
    
    if (nrow(ivw_result) > 0) {
      ivw_beta <- safe_extract(ivw_result$b, 1, NA)
      ivw_se <- safe_extract(ivw_result$se, 1, NA)
      ivw_pval <- safe_extract(ivw_result$pval, 1, NA)
      ivw_nsnp_val <- safe_extract(ivw_result$nsnp, 1, result$n_snps_harmonized)
      ivw_nsnp <- ifelse(is.na(ivw_nsnp_val), result$n_snps_harmonized, ivw_nsnp_val)
      
      # 检查是否有效值
      if (is.na(ivw_beta) || is.na(ivw_se)) {
        next  # 跳过无效结果
      }
      
      ivw_or <- exp(ivw_beta)
      ivw_or_lci <- exp(ivw_beta - 1.96 * ivw_se)
      ivw_or_uci <- exp(ivw_beta + 1.96 * ivw_se)
      
      # 提取其他方法
      egger_result <- result$mr_results[result$mr_results$method == "MR Egger", ]
      wm_result <- result$mr_results[result$mr_results$method == "Weighted median", ]
      wmode_result <- result$mr_results[result$mr_results$method == "Weighted mode", ]
      smode_result <- result$mr_results[result$mr_results$method == "Simple mode", ]
      
      # 异质性检验
      het_p <- NA
      if (!is.null(result$heterogeneity) && nrow(result$heterogeneity) > 0) {
        ivw_het <- result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", ]
        if (nrow(ivw_het) > 0) {
          het_p <- safe_extract(ivw_het$Q_pval, 1, NA)
        }
      }
      
      # 多效性检验
      pleo_p <- NA
      pleo_intercept <- NA
      if (!is.null(result$pleiotropy) && nrow(result$pleiotropy) > 0) {
        pleo_p <- safe_extract(result$pleiotropy$pval, 1, NA)
        pleo_intercept <- safe_extract(result$pleiotropy$egger_intercept, 1, NA)
      }
      
      # MR-PRESSO结果
      presso_global_p <- NA
      presso_outlier_corrected_or <- NA
      presso_outlier_corrected_or_lci <- NA
      presso_outlier_corrected_or_uci <- NA
      presso_n_outliers <- 0
      
      if (!is.null(result$presso)) {
        tryCatch({
          presso_global_p <- result$presso$`MR-PRESSO results`$`Global Test`$Pvalue
          if (!is.null(result$presso$`MR-PRESSO results`$`Outlier Test`)) {
            presso_n_outliers <- nrow(result$presso$`MR-PRESSO results`$`Outlier Test`)
          }
          if (!is.null(result$presso$`Main MR results`)) {
            presso_beta <- result$presso$`Main MR results`[1, "Causal Estimate"]
            presso_se <- result$presso$`Main MR results`[1, "Sd"]
            presso_outlier_corrected_or <- exp(presso_beta)
            presso_outlier_corrected_or_lci <- exp(presso_beta - 1.96 * presso_se)
            presso_outlier_corrected_or_uci <- exp(presso_beta + 1.96 * presso_se)
          }
        }, error = function(e) {
          # 静默处理
        })
      }
      
      # Leave-one-out稳定性
      loo_stable <- "NA"
      loo_all_significant <- FALSE
      if (!is.null(result$loo) && nrow(result$loo) > 1) {
        loo_pvals <- result$loo$p[!is.na(result$loo$p)]
        if (length(loo_pvals) > 0) {
          all_significant <- all(loo_pvals < 0.05, na.rm = TRUE)
          all_non_significant <- all(loo_pvals >= 0.05, na.rm = TRUE)
          
          if (all_significant && ivw_pval < 0.05) {
            loo_stable <- "Stable (all significant)"
            loo_all_significant <- TRUE
          } else if (all_non_significant && ivw_pval >= 0.05) {
            loo_stable <- "Stable (all non-significant)"
          } else {
            loo_stable <- "Unstable"
          }
        }
      }
      
      # 工具变量强度
      mean_f <- ifelse(!is.null(result$iv_strength), result$iv_strength$mean_f, NA)
      r_squared <- ifelse(!is.null(result$iv_strength), result$iv_strength$r_squared, NA)
      
      # 安全提取Egger结果
      egger_beta_val <- if (nrow(egger_result) > 0) safe_extract(egger_result$b, 1, NA) else NA
      egger_se_val <- if (nrow(egger_result) > 0) safe_extract(egger_result$se, 1, NA) else NA
      egger_pval_val <- if (nrow(egger_result) > 0) safe_extract(egger_result$pval, 1, NA) else NA
      egger_or_95ci_val <- if (nrow(egger_result) > 0 && !is.na(egger_beta_val) && !is.na(egger_se_val)) {
        egger_or <- exp(egger_beta_val)
        sprintf("%.3f (%.3f-%.3f)", 
               egger_or, 
               exp(egger_beta_val - 1.96 * egger_se_val),
               exp(egger_beta_val + 1.96 * egger_se_val))
      } else {
        NA
      }
      
      # 安全提取Weighted Median结果
      wm_beta_val <- if (nrow(wm_result) > 0) safe_extract(wm_result$b, 1, NA) else NA
      wm_se_val <- if (nrow(wm_result) > 0) safe_extract(wm_result$se, 1, NA) else NA
      wm_pval_val <- if (nrow(wm_result) > 0) safe_extract(wm_result$pval, 1, NA) else NA
      wm_or_95ci_val <- if (nrow(wm_result) > 0 && !is.na(wm_beta_val) && !is.na(wm_se_val)) {
        wm_or <- exp(wm_beta_val)
        sprintf("%.3f (%.3f-%.3f)", 
               wm_or, 
               exp(wm_beta_val - 1.96 * wm_se_val),
               exp(wm_beta_val + 1.96 * wm_se_val))
      } else {
        NA
      }
      
      # Table 2增强版
      table2_enhanced_data[[length(table2_enhanced_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        method = "IVW",
        n_snps = ivw_nsnp,
        or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
        pval = ivw_pval,
        egger_or_95ci = egger_or_95ci_val,
        egger_pval = egger_pval_val,
        wm_or_95ci = wm_or_95ci_val,
        wm_pval = wm_pval_val,
        pleiotropy_p = pleo_p,
        heterogeneity_p = het_p,
        stringsAsFactors = FALSE
      )
      
      # Table 5: 敏感性分析汇总
      table5_data[[length(table5_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        n_snps = ivw_nsnp,
        ivw_or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
        ivw_pval = ivw_pval,
        egger_intercept_p = pleo_p,
        heterogeneity_q_p = het_p,
        presso_global_p = presso_global_p,
        presso_n_outliers = presso_n_outliers,
        presso_outlier_corrected_or_95ci = ifelse(!is.na(presso_outlier_corrected_or),
                                                  sprintf("%.3f (%.3f-%.3f)", 
                                                         presso_outlier_corrected_or,
                                                         presso_outlier_corrected_or_lci,
                                                         presso_outlier_corrected_or_uci),
                                                  NA),
        loo_stability = loo_stable,
        stringsAsFactors = FALSE
      )
      
      # 安全提取异质性Q值
      het_q <- NA
      if (!is.null(result$heterogeneity) && nrow(result$heterogeneity) > 0) {
        ivw_het <- result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", ]
        if (nrow(ivw_het) > 0) {
          het_q <- safe_extract(ivw_het$Q, 1, NA)
        }
      }
      
      # MR-PRESSO结果字符串
      presso_or_95ci_str <- if (!is.na(presso_outlier_corrected_or)) {
        sprintf("%.3f (%.3f-%.3f)", 
               presso_outlier_corrected_or,
               presso_outlier_corrected_or_lci,
               presso_outlier_corrected_or_uci)
      } else {
        NA
      }
      
      # Table S2: 敏感性分析详情
      table_s2_data[[length(table_s2_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        method = "IVW",
        beta = ivw_beta,
        se = ivw_se,
        pval = ivw_pval,
        or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
        egger_beta = egger_beta_val,
        egger_se = egger_se_val,
        egger_pval = egger_pval_val,
        wm_beta = wm_beta_val,
        wm_se = wm_se_val,
        wm_pval = wm_pval_val,
        pleiotropy_intercept = pleo_intercept,
        pleiotropy_p = pleo_p,
        heterogeneity_q = het_q,
        heterogeneity_q_p = het_p,
        presso_global_p = presso_global_p,
        presso_outlier_corrected_or_95ci = presso_or_95ci_str,
        stringsAsFactors = FALSE
      )
      
      # Table S6: 工具变量强度
      table_s6_data[[length(table_s6_data) + 1]] <- data.frame(
        category = category,
        exposure = exposure,
        outcome = outcome,
        n_snps = ivw_nsnp,
        mean_f_statistic = mean_f,
        r_squared = r_squared,
        extraction_strategy = ifelse(is.na(strategy), "Unknown", strategy),
        stringsAsFactors = FALSE
      )
      
      # Table S3: Leave-one-out详情
      if (!is.null(result$loo) && nrow(result$loo) > 0) {
        for (j in seq_len(nrow(result$loo))) {
          loo_row <- result$loo[j, ]
          loo_beta <- loo_row$b
          loo_or <- exp(loo_beta)
          loo_or_lci <- exp(loo_beta - 1.96 * loo_row$se)
          loo_or_uci <- exp(loo_beta + 1.96 * loo_row$se)
          loo_snp_name <- ifelse(is.null(loo_row$SNP) || is.na(loo_row$SNP), 
                                 paste("SNP", j), loo_row$SNP)
          loo_nsnp <- ifelse(is.null(loo_row$nsnp) || is.na(loo_row$nsnp), 
                            NA, loo_row$nsnp)
          
          table_s3_data[[length(table_s3_data) + 1]] <- data.frame(
            category = category,
            exposure = exposure,
            outcome = outcome,
            snp_removed = loo_snp_name,
            n_snps_remaining = loo_nsnp,
            or_95ci = sprintf("%.3f (%.3f-%.3f)", loo_or, loo_or_lci, loo_or_uci),
            pval = loo_row$p,
            significant = ifelse(loo_row$p < 0.05, "Yes", "No"),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  # 合并所有表格
  table_s1_full <- if (length(table_s1_data) > 0) {
    do.call(rbind, table_s1_data)
  } else {
    data.frame()
  }
  
  table2_enhanced_full <- if (length(table2_enhanced_data) > 0) {
    do.call(rbind, table2_enhanced_data)
  } else {
    data.frame()
  }
  
  table5_full <- if (length(table5_data) > 0) {
    do.call(rbind, table5_data)
  } else {
    data.frame()
  }
  
  table_s2_full <- if (length(table_s2_data) > 0) {
    do.call(rbind, table_s2_data)
  } else {
    data.frame()
  }
  
  table_s3_full <- if (length(table_s3_data) > 0) {
    do.call(rbind, table_s3_data)
  } else {
    data.frame()
  }
  
  table_s6_full <- if (length(table_s6_data) > 0) {
    do.call(rbind, table_s6_data)
  } else {
    data.frame()
  }
  
  # 添加显著性标记并排序
  if (nrow(table_s1_full) > 0) {
    table_s1_full$significance <- ifelse(table_s1_full$pval < 0.001, "***",
                                         ifelse(table_s1_full$pval < 0.01, "**",
                                                ifelse(table_s1_full$pval < 0.05, "*", "")))
    table_s1_full <- table_s1_full[order(table_s1_full$category, 
                                        table_s1_full$exposure, 
                                        table_s1_full$outcome,
                                        table_s1_full$method), ]
  }
  
  if (nrow(table2_enhanced_full) > 0) {
    table2_enhanced_full$ivw_significant <- ifelse(table2_enhanced_full$pval < 0.05, "Yes", "No")
    table2_enhanced_full <- table2_enhanced_full[order(table2_enhanced_full$category, 
                                                      table2_enhanced_full$exposure, 
                                                      table2_enhanced_full$outcome), ]
  }
  
  # 保存为Excel和CSV
  if (nrow(table_s1_full) > 0) {
    write.csv(table_s1_full, "results/tables/paper_tables/Table_S1_Reverse_MR_Full_Results.csv", row.names = FALSE)
  }
  if (nrow(table2_enhanced_full) > 0) {
    write.csv(table2_enhanced_full, "results/tables/paper_tables/Table_2_Reverse_MR_Enhanced.csv", row.names = FALSE)
  }
  if (nrow(table5_full) > 0) {
    write.csv(table5_full, "results/tables/paper_tables/Table_5_Reverse_MR_Sensitivity.csv", row.names = FALSE)
  }
  if (nrow(table_s2_full) > 0) {
    write.csv(table_s2_full, "results/tables/paper_tables/Table_S2_Reverse_MR_Sensitivity_Details.csv", row.names = FALSE)
  }
  if (nrow(table_s3_full) > 0) {
    write.csv(table_s3_full, "results/tables/paper_tables/Table_S3_Reverse_MR_Leave_One_Out.csv", row.names = FALSE)
  }
  if (nrow(table_s6_full) > 0) {
    write.csv(table_s6_full, "results/tables/paper_tables/Table_S6_Reverse_MR_IV_Strength.csv", row.names = FALSE)
  }
  
  # 创建多工作表Excel文件
  wb <- createWorkbook()
  if (nrow(table_s1_full) > 0) {
    addWorksheet(wb, "Table_S1_Full_MR")
    writeData(wb, "Table_S1_Full_MR", table_s1_full)
  }
  if (nrow(table2_enhanced_full) > 0) {
    addWorksheet(wb, "Table_2_Enhanced")
    writeData(wb, "Table_2_Enhanced", table2_enhanced_full)
  }
  if (nrow(table5_full) > 0) {
    addWorksheet(wb, "Table_5_Sensitivity")
    writeData(wb, "Table_5_Sensitivity", table5_full)
  }
  if (nrow(table_s2_full) > 0) {
    addWorksheet(wb, "Table_S2_Details")
    writeData(wb, "Table_S2_Details", table_s2_full)
  }
  if (nrow(table_s3_full) > 0) {
    addWorksheet(wb, "Table_S3_LOO")
    writeData(wb, "Table_S3_LOO", table_s3_full)
  }
  if (nrow(table_s6_full) > 0) {
    addWorksheet(wb, "Table_S6_IV_Strength")
    writeData(wb, "Table_S6_IV_Strength", table_s6_full)
  }
  
  saveWorkbook(wb, "results/tables/paper_tables/step11_all_paper_tables.xlsx", overwrite = TRUE)
  
  cat("✓ 论文表格已保存到 results/tables/paper_tables/\n")
  cat("  - step11_all_paper_tables.xlsx (多工作表)\n")
  cat("  - Table_S1_Reverse_MR_Full_Results.csv\n")
  cat("  - Table_2_Reverse_MR_Enhanced.csv\n")
  cat("  - Table_5_Reverse_MR_Sensitivity.csv\n")
  cat("  - Table_S2_Reverse_MR_Sensitivity_Details.csv\n")
  cat("  - Table_S3_Reverse_MR_Leave_One_Out.csv\n")
  cat("  - Table_S6_Reverse_MR_IV_Strength.csv\n\n")
}

# ============================================================================
# 步骤9：生成主文可视化（森林图和热图）
# ============================================================================

if (nrow(results_summary) > 0) {
  cat("【步骤9】生成主文可视化（森林图和热图）...\n")
  
  # 9.1 生成森林图（Forest Plot）
  forest_data <- results_summary[results_summary$significant_nominal == TRUE, ]
  
  if (nrow(forest_data) > 0) {
    # 按p值排序
    forest_data <- forest_data[order(forest_data$pval), ]
    
    # 准备森林图数据
    forest_data$label <- paste(forest_data$exposure, "→", forest_data$outcome)
    forest_data$label <- factor(forest_data$label, levels = rev(forest_data$label))
    
    # 创建森林图
    tryCatch({
      p_forest <- ggplot(forest_data, aes(x = label, y = or, ymin = or_lci, ymax = or_uci, 
                                          color = category)) +
        geom_pointrange(size = 0.5, fatten = 3) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c("Metabolic" = "#E69F00", "Inflammatory" = "#56B4E9")) +
        coord_flip() +
        labs(
          title = "Forest Plot: Significant Reverse MR Associations",
          subtitle = "Lung Cancer → Metabolic/Inflammatory Traits",
          x = "Exposure → Outcome",
          y = "Odds Ratio (95% CI)",
          color = "Category"
        ) +
        theme_classic(base_size = 10) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 9),
          legend.position = "right"
        )
      
      ggsave("results/figures/reverse_mr/main_figures/Figure1_ReverseMR_ForestPlot.pdf", 
             plot = p_forest, width = 10, height = max(8, nrow(forest_data) * 0.4), 
             dpi = 600, limitsize = FALSE)
      ggsave("results/figures/reverse_mr/main_figures/Figure1_ReverseMR_ForestPlot.png", 
             plot = p_forest, width = 10, height = max(8, nrow(forest_data) * 0.4), 
             dpi = 600, limitsize = FALSE)
      
      cat("✓ 森林图已生成:\n")
      cat("  - results/figures/reverse_mr/main_figures/Figure1_ReverseMR_ForestPlot.pdf\n")
      cat("  - results/figures/reverse_mr/main_figures/Figure1_ReverseMR_ForestPlot.png\n")
    }, error = function(e) {
      cat("⚠ 森林图生成失败:", e$message, "\n")
    })
  }
  
  # 9.2 生成热图（Heatmap）
  tryCatch({
    # 准备热图数据：按暴露分组
    heatmap_data <- results_summary %>%
      select(exposure, outcome, or, pval, category) %>%
      mutate(
        log_or = log(or),
        significant = ifelse(pval < 0.05, 1, 0)
      )
    
    if (require("tidyr", quietly = TRUE) && require("pheatmap", quietly = TRUE)) {
      # 创建OR值热图矩阵
      or_matrix <- heatmap_data %>%
        select(exposure, outcome, log_or) %>%
        pivot_wider(names_from = outcome, values_from = log_or, values_fill = NA)
      
      or_mat <- as.matrix(or_matrix[, -1])
      rownames(or_mat) <- or_matrix$exposure
      
      # 创建显著性矩阵
      sig_matrix <- heatmap_data %>%
        select(exposure, outcome, significant) %>%
        pivot_wider(names_from = outcome, values_from = significant, values_fill = 0)
      
      sig_mat <- as.matrix(sig_matrix[, -1])
      rownames(sig_mat) <- sig_matrix$exposure
      
      # 生成OR值热图
      pdf("results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_OR.pdf", width = 12, height = 6)
      pheatmap::pheatmap(
        or_mat,
        color = RColorBrewer::brewer.pal(11, "RdBu"),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        main = "Reverse MR Association Heatmap (log OR)\nLung Cancer → Metabolic/Inflammatory Traits",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 7
      )
      dev.off()
      
      png("results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_OR.png", 
          width = 12, height = 6, units = "in", res = 600)
      pheatmap::pheatmap(
        or_mat,
        color = RColorBrewer::brewer.pal(11, "RdBu"),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        main = "Reverse MR Association Heatmap (log OR)\nLung Cancer → Metabolic/Inflammatory Traits",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 7
      )
      dev.off()
      
      # 生成显著性热图
      pdf("results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_Significance.pdf", width = 12, height = 6)
      pheatmap::pheatmap(
        sig_mat,
        color = c("white", "red"),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        main = "Reverse MR Association Significance (p < 0.05)\nLung Cancer → Metabolic/Inflammatory Traits",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 7
      )
      dev.off()
      
      png("results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_Significance.png", 
          width = 12, height = 6, units = "in", res = 600)
      pheatmap::pheatmap(
        sig_mat,
        color = c("white", "red"),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "none",
        main = "Reverse MR Association Significance (p < 0.05)\nLung Cancer → Metabolic/Inflammatory Traits",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 7
      )
      dev.off()
      
      cat("✓ 热图已生成:\n")
      cat("  - results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_OR.pdf\n")
      cat("  - results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_OR.png\n")
      cat("  - results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_Significance.pdf\n")
      cat("  - results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_Significance.png\n")
    }
  }, error = function(e) {
    cat("⚠ 热图生成失败:", e$message, "\n")
  })
  cat("\n")
}

# ============================================================================
# 步骤10：生成分析报告
# ============================================================================

cat("【步骤10】生成分析报告...\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

capture.output({
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("反向MR分析完整报告\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  cat("【分析覆盖度】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
  cat(sprintf("预期总分析数:   %d 个 (%d暴露 × %d结局)\n", 
              total_pairs, length(lung_cancer_datasets), length(all_outcomes)))
  cat(sprintf("实际尝试分析:   %d 个\n", analysis_counter))
  cat(sprintf("成功完成分析:   %d 个\n", success_count))
  cat(sprintf("成功率:         %.1f%%\n\n", 100 * success_count / max(1, analysis_counter)))
  
  cat("【分析结果统计】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  if (nrow(results_summary) > 0) {
    cat(sprintf("完整结果数:     %d 个\n", nrow(results_summary)))
    cat(sprintf("FDR显著 (<0.05): %d 个\n", sum(results_summary$significant_fdr, na.rm = TRUE)))
    cat(sprintf("名义显著 (<0.05): %d 个\n\n", sum(results_summary$significant_nominal, na.rm = TRUE)))
    
    # FDR显著关联
    sig_fdr <- results_summary[results_summary$significant_fdr == TRUE, ]
    if (nrow(sig_fdr) > 0) {
      cat("【FDR显著关联详情】\n")
      cat(paste(rep("-", 80), collapse = ""), "\n")
      sig_fdr <- sig_fdr[order(sig_fdr$fdr_pval), ]
      for (i in seq_len(min(10, nrow(sig_fdr)))) {
        cat(sprintf("%d. %s → %s\n", i, sig_fdr$exposure[i], sig_fdr$outcome[i]))
        cat(sprintf("   OR = %.3f (95%%CI: %.3f-%.3f)\n", sig_fdr$or[i], sig_fdr$or_lci[i], sig_fdr$or_uci[i]))
        cat(sprintf("   P = %.2e, FDR = %.2e\n", sig_fdr$pval[i], sig_fdr$fdr_pval[i]))
        cat(sprintf("   SNPs: %d (协调: %d), 策略: %s\n\n", 
                   sig_fdr$n_snps[i], sig_fdr$n_harmonized[i], sig_fdr$extraction_strategy[i]))
      }
    }
  }
  
    cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("报告生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
}, file = "results/step11_analysis_report.txt")

cat("✓ 分析报告已保存: results/step11_analysis_report.txt\n\n")

# ============================================================================
# 最终总结
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("第11步分析完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【最终统计】\n")
cat(sprintf("预期分析:       %d 个\n", total_pairs))
cat(sprintf("实际尝试:       %d 个\n", analysis_counter))
cat(sprintf("成功完成:       %d 个\n", success_count))
cat(sprintf("成功率:         %.1f%%\n", 100 * success_count / max(1, analysis_counter)))
if (nrow(results_summary) > 0) {
  cat(sprintf("FDR显著关联:    %d 个\n", sum(results_summary$significant_fdr, na.rm = TRUE)))
  cat(sprintf("名义显著关联:   %d 个\n\n", sum(results_summary$significant_nominal, na.rm = TRUE)))
}

cat("【保存的文件】\n")
cat("主要结果:\n")
cat("  - data/step11_reverse_mr_complete_results.RData\n")
cat("汇总表格:\n")
cat("  - results/tables/step11_mr_results_summary.xlsx\n")
cat("  - results/tables/step11_mr_results_summary.csv\n")
cat("  - results/tables/step11_extraction_log.csv\n")
cat("详细报告:\n")
cat("  - results/step11_analysis_report.txt\n")
cat("图表文件:\n")
cat("  - results/figures/reverse_mr/main_figures/Figure1_ReverseMR_ForestPlot.pdf/png\n")
cat("  - results/figures/reverse_mr/main_figures/Figure2_ReverseMR_Heatmap_*.pdf/png\n")
cat("敏感性分析图表:\n")
cat("  - results/figures/reverse_mr/sensitivity/*.png\n")
cat("论文表格:\n")
cat("  - results/tables/paper_tables/step11_all_paper_tables.xlsx\n")
cat("  - results/tables/paper_tables/*.csv\n\n")

cat("分析完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n提示: 可以使用以下命令查看结果:\n")
cat("  - load('data/step11_reverse_mr_complete_results.RData')\n")
cat("  - View(results_summary)\n")
cat("  - readLines('results/step11_analysis_report.txt') %>% cat(sep='\\n')\n")

# ============================================================================
# 步骤11：正向vs反向MR对比分析（SCI期刊10分标准图表）
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("【步骤11】正向vs反向MR对比分析（SCI期刊10分标准）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 加载SCI标准图表函数库（从step07c）
if (!require("scales", quietly = TRUE)) {
  install.packages("scales", repos = "https://cloud.r-project.org")
}
library(scales)

# Okabe-Ito 色盲友好调色板
okabe_ito_colors <- list(
  orange = "#E69F00",         # 正向MR
  sky_blue = "#56B4E9",       # 反向MR
  green = "#009E73",          # 双向显著
  vermillion = "#D55E00",     # 无显著
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
                     base_line_size = 0.21,
                     base_rect_size = 0.21,
                     grid = TRUE) {
  
  if (.Platform$OS.type == "windows") {
    if (base_family == "Arial") {
      font_available <- tryCatch({
        if (requireNamespace("extrafont", quietly = TRUE)) {
          "Arial" %in% extrafont::fonts()
        } else {
          FALSE
        }
      }, error = function(e) FALSE)
      
      if (!font_available) {
        base_family <- ""
      }
    }
  } else {
    font_available <- tryCatch({
      if (requireNamespace("extrafont", quietly = TRUE)) {
        base_family %in% extrafont::fonts()
      } else {
        TRUE
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
      plot.title = element_text(
        size = base_size + 2,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5, unit = "mm"),
        color = "black"
      ),
      plot.subtitle = element_text(
        size = base_size,
        hjust = 0.5,
        margin = margin(b = 8, unit = "mm"),
        color = "black"
      ),
      plot.caption = element_text(
        size = base_size - 1,
        hjust = 0.5,
        margin = margin(t = 5, unit = "mm"),
        color = "gray40"
      ),
      axis.title = element_text(
        size = base_size,
        face = "bold",
        color = "black"
      ),
      axis.title.x = element_text(margin = margin(t = 5, unit = "mm")),
      axis.title.y = element_text(margin = margin(r = 5, unit = "mm"), angle = 90),
      axis.text = element_text(
        size = base_size - 1,
        color = "black"
      ),
      axis.text.x = element_text(margin = margin(t = 2, unit = "mm")),
      axis.text.y = element_text(margin = margin(r = 2, unit = "mm")),
      legend.title = element_text(
        size = base_size,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = base_size - 1,
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
      panel.grid.major = if (grid) {
        element_line(
          color = "gray90",
          linewidth = 0.11,
          linetype = "solid"
        )
      } else {
        element_blank()
      },
      panel.grid.minor = element_blank(),
      axis.line = element_line(
        color = "black",
        linewidth = base_line_size
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = base_line_size
      ),
      axis.ticks.length = unit(2, "mm"),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      plot.margin = margin(10, 10, 10, 10, unit = "mm"),
      strip.text = element_text(
        size = base_size,
        face = "bold",
        color = "black"
      ),
      strip.background = element_rect(
        fill = "gray95",
        color = "black",
        linewidth = base_line_size
      ),
      complete = TRUE
    )
}

# 保存SCI标准图表函数
save_sci_figure <- function(plot, 
                           filename,
                           width_mm = 174,
                           height_mm = 120,
                           dpi = 600,
                           formats = c("png", "pdf")) {
  
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  if ("png" %in% formats) {
    ggsave(
      paste0(filename, ".png"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi,
      bg = "white",
      type = "cairo"
    )
    cat(sprintf("✓ PNG: %s.png (%.0f × %.0f mm, %d DPI)\n", 
                filename, width_mm, height_mm, dpi))
  }
  
  if ("pdf" %in% formats) {
    ggsave(
      paste0(filename, ".pdf"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      device = "pdf",
      useDingbats = FALSE
    )
    cat(sprintf("✓ PDF: %s.pdf (%.0f × %.0f mm, 矢量格式)\n", 
                filename, width_mm, height_mm))
  }
}

cat("【步骤11.1】加载正向MR（第5步）结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 加载正向MR结果
forward_mr_results <- NULL
tryCatch({
  if (file.exists("results/tables/step05_mr_results_summary.csv")) {
    forward_mr_results <- read.csv("results/tables/step05_mr_results_summary.csv",
                                   stringsAsFactors = FALSE)
    cat(sprintf("✓ 已加载正向MR结果: %d 个分析\n", nrow(forward_mr_results)))
  } else {
    cat("⚠ 警告：未找到正向MR结果文件（results/tables/step05_mr_results_summary.csv）\n")
    cat("  将跳过对比分析\n")
  }
}, error = function(e) {
  cat(sprintf("⚠ 加载正向MR结果失败: %s\n", e$message))
})

# 检查反向MR结果
if (nrow(results_summary) == 0) {
  cat("⚠ 警告：反向MR结果为空，无法进行对比分析\n\n")
} else if (is.null(forward_mr_results) || nrow(forward_mr_results) == 0) {
  cat("⚠ 警告：正向MR结果为空，无法进行对比分析\n\n")
} else {
  
  cat("\n【步骤11.2】准备对比数据\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # 标准化暴露和结局名称以进行匹配
  # 正向MR: exposure -> outcome (代谢/炎症 -> 肺癌)
  # 反向MR: exposure -> outcome (肺癌 -> 代谢/炎症)
  
  # 创建正向MR配对（标准化名称）
  forward_mr_clean <- forward_mr_results %>%
    mutate(
      # 提取IVW结果
      forward_or = exp(beta),
      forward_or_lci = exp(beta - 1.96 * se),
      forward_or_uci = exp(beta + 1.96 * se),
      forward_pval = pval,
      forward_sig = pval < 0.05,
      forward_fdr_sig = if("significant_fdr" %in% names(.)) significant_fdr else FALSE,
      # 创建配对键：正向是 暴露_结局
      pair_key_forward = paste(exposure, outcome, sep = "_")
    ) %>%
    select(exposure, outcome, pair_key_forward, forward_or, forward_or_lci, forward_or_uci,
           forward_pval, forward_sig, forward_fdr_sig, n_snps, n_harmonized)
  
  # 创建反向MR配对（标准化名称）
  reverse_mr_clean <- results_summary %>%
    mutate(
      # 反向MR的暴露是肺癌，结局是代谢/炎症
      # 所以配对键应该是：正向的结局_正向的暴露
      # 即：肺癌类型_代谢性状
      pair_key_reverse = paste(exposure, outcome, sep = "_"),
      reverse_or = or,
      reverse_or_lci = or_lci,
      reverse_or_uci = or_uci,
      reverse_pval = pval,
      reverse_sig = significant_nominal,
      reverse_fdr_sig = if("significant_fdr" %in% names(.)) significant_fdr else FALSE
    ) %>%
    select(exposure, outcome, pair_key_reverse, reverse_or, reverse_or_lci, reverse_or_uci,
           reverse_pval, reverse_sig, reverse_fdr_sig, n_snps, n_harmonized)
  
  # 匹配策略：找到共同的暴露-结局配对
  # 正向: 代谢性状 -> 肺癌类型
  # 反向: 肺癌类型 -> 代谢性状
  # 所以需要：正向的(exposure, outcome) 与 反向的(exposure, outcome)互换匹配
  
  # 创建匹配键
  forward_match <- forward_mr_clean %>%
    mutate(
      match_key = paste(outcome, exposure, sep = "_")  # 肺癌_代谢性状
    )
  
  reverse_match <- reverse_mr_clean %>%
    mutate(
      match_key = paste(exposure, outcome, sep = "_")  # 肺癌_代谢性状
    )
  
  # 合并数据
  comparison_data <- forward_match %>%
    full_join(reverse_match, by = "match_key", suffix = c("_forward", "_reverse")) %>%
    filter(!is.na(exposure_forward) | !is.na(exposure_reverse)) %>%
    mutate(
      # 标准化变量名
      metabolic_trait = ifelse(!is.na(exposure_forward), exposure_forward, outcome_reverse),
      lung_cancer = ifelse(!is.na(outcome_forward), outcome_forward, exposure_reverse),
      # 完整的配对键
      full_pair_key = match_key,
      # 判断显著性类型
      both_significant = (!is.na(forward_sig) && forward_sig) && (!is.na(reverse_sig) && reverse_sig),
      forward_only_sig = (!is.na(forward_sig) && forward_sig) && (is.na(reverse_sig) || !reverse_sig),
      reverse_only_sig = (is.na(forward_sig) || !forward_sig) && (!is.na(reverse_sig) && reverse_sig),
      neither_sig = (is.na(forward_sig) || !forward_sig) && (is.na(reverse_sig) || !reverse_sig),
      # 显著性类别
      significance_category = case_when(
        both_significant ~ "Bidirectional",
        forward_only_sig ~ "Forward only",
        reverse_only_sig ~ "Reverse only",
        neither_sig ~ "Neither",
        TRUE ~ "Unknown"
      )
    )
  
  cat(sprintf("✓ 已创建对比数据集: %d 个配对\n", nrow(comparison_data)))
  cat(sprintf("  - 双向显著: %d 个\n", sum(comparison_data$both_significant, na.rm = TRUE)))
  cat(sprintf("  - 仅正向显著: %d 个\n", sum(comparison_data$forward_only_sig, na.rm = TRUE)))
  cat(sprintf("  - 仅反向显著: %d 个\n", sum(comparison_data$reverse_only_sig, na.rm = TRUE)))
  cat(sprintf("  - 均不显著: %d 个\n", sum(comparison_data$neither_sig, na.rm = TRUE)))
  
  cat("\n【步骤11.3】生成对比图表（SCI期刊10分标准）\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # 创建输出目录
  dir.create("results/figures/reverse_mr/comparison", showWarnings = FALSE, recursive = TRUE)
  
  # 图表1：效应值对比散点图
  cat("  生成图表1：效应值对比散点图...\n")
  tryCatch({
    comparison_plot_data <- comparison_data %>%
      filter(!is.na(forward_or), !is.na(reverse_or)) %>%
      mutate(
        forward_log_or = log(forward_or),
        reverse_log_or = log(reverse_or)
      )
    
    if (nrow(comparison_plot_data) > 0) {
      p_comparison_scatter <- ggplot(comparison_plot_data, 
                                    aes(x = forward_log_or, y = reverse_log_or,
                                        color = significance_category)) +
        geom_point(alpha = 0.7, size = 2.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.21) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.21) +
        geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "gray70", linewidth = 0.21) +
        scale_color_manual(
          name = "Significance",
          values = c(
            "Bidirectional" = okabe_ito$green,
            "Forward only" = okabe_ito$orange,
            "Reverse only" = okabe_ito$sky_blue,
            "Neither" = okabe_ito$gray,
            "Unknown" = okabe_ito$gray
          ),
          guide = guide_legend(title.position = "top")
        ) +
        labs(
          title = "Figure 3. Forward vs Reverse MR Effect Comparison",
          subtitle = sprintf("Comparison of OR values between forward (Metabolic→Lung Cancer) and reverse (Lung Cancer→Metabolic) MR analyses\n(n = %d pairs)", nrow(comparison_plot_data)),
          x = "Forward MR log(OR)",
          y = "Reverse MR log(OR)",
          caption = "Points on the diagonal line indicate similar effect magnitudes. Bidirectional associations suggest potential pleiotropy or sample overlap."
        ) +
        theme_sci() +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right"
        )
      
      save_sci_figure(
        p_comparison_scatter,
        "results/figures/reverse_mr/comparison/Figure3_Forward_vs_Reverse_Scatter",
        width_mm = 174,
        height_mm = 140
      )
      
      cat("    ✓ 散点图已保存\n")
    } else {
      cat("    ⚠ 无数据可用于散点图\n")
    }
  }, error = function(e) {
    cat(sprintf("    ⚠ 散点图生成失败: %s\n", e$message))
  })
  
  # 图表2：显著性对比条形图
  cat("  生成图表2：显著性对比条形图...\n")
  tryCatch({
    sig_summary <- comparison_data %>%
    summarise(
        Bidirectional = sum(both_significant, na.rm = TRUE),
        `Forward only` = sum(forward_only_sig, na.rm = TRUE),
        `Reverse only` = sum(reverse_only_sig, na.rm = TRUE),
        Neither = sum(neither_sig, na.rm = TRUE)
      ) %>%
      pivot_longer(everything(), names_to = "Category", values_to = "Count") %>%
      mutate(
        Category = factor(Category, levels = c("Bidirectional", "Forward only", "Reverse only", "Neither")),
        Percentage = 100 * Count / sum(Count)
      )
    
    p_sig_comparison <- ggplot(sig_summary, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
      geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
                vjust = -0.2, size = 3.5, fontface = "bold", lineheight = 0.9) +
      scale_fill_manual(
        name = "Significance Category",
        values = c(
          "Bidirectional" = okabe_ito$green,
          "Forward only" = okabe_ito$orange,
          "Reverse only" = okabe_ito$sky_blue,
          "Neither" = okabe_ito$gray
        ),
        guide = "none"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "Figure 4. Significance Pattern Comparison",
        subtitle = "Distribution of significance patterns in forward and reverse MR analyses",
        x = "Significance Category",
        y = "Number of Associations",
        caption = NULL
      ) +
      theme_sci() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    save_sci_figure(
      p_sig_comparison,
      "results/figures/reverse_mr/comparison/Figure4_Significance_Comparison",
      width_mm = 174,
      height_mm = 100
    )
    
    cat("    ✓ 显著性对比图已保存\n")
  }, error = function(e) {
    cat(sprintf("    ⚠ 显著性对比图生成失败: %s\n", e$message))
  })
  
  # 图表3：双向关联检测森林图（仅显示双向显著的）
  cat("  生成图表3：双向关联检测森林图...\n")
  tryCatch({
    bidirectional_data <- comparison_data %>%
      filter(both_significant) %>%
      arrange(desc(forward_pval)) %>%
      head(20) %>%
      mutate(
        label = paste(metabolic_trait, "↔", lung_cancer),
        label = factor(label, levels = rev(label)),
        forward_log_or = log(forward_or),
        reverse_log_or = log(reverse_or)
      )
    
    if (nrow(bidirectional_data) > 0) {
      # 准备长格式数据用于绘图
      plot_data_long <- bidirectional_data %>%
        select(label, forward_log_or, reverse_log_or, forward_pval, reverse_pval) %>%
        pivot_longer(
          cols = c(forward_log_or, reverse_log_or),
          names_to = "direction",
          values_to = "log_or"
        ) %>%
        mutate(
          direction_label = ifelse(direction == "forward_log_or", "Forward MR", "Reverse MR"),
          pval = ifelse(direction == "forward_log_or", forward_pval, reverse_pval),
          significant = pval < 0.05
        )
      
      p_bidirectional <- ggplot(plot_data_long, 
                                aes(x = log_or, y = label, 
                                    color = direction_label, shape = direction_label)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.21) +
        scale_color_manual(
          name = "MR Direction",
          values = c(
            "Forward MR" = okabe_ito$orange,
            "Reverse MR" = okabe_ito$sky_blue
          )
        ) +
        scale_shape_manual(
          name = "MR Direction",
          values = c("Forward MR" = 16, "Reverse MR" = 17)
        ) +
        labs(
          title = "Figure 5. Bidirectional MR Associations",
          subtitle = sprintf("Associations significant in both directions (P < 0.05)\n(n = %d associations)", nrow(bidirectional_data)),
          x = "log(OR)",
          y = "Association",
          caption = "Points represent log(OR) values. Associations showing significance in both directions may indicate pleiotropy or sample overlap."
        ) +
        theme_sci() +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right"
        )
      
      save_sci_figure(
        p_bidirectional,
        "results/figures/reverse_mr/comparison/Figure5_Bidirectional_Associations",
        width_mm = 174,
        height_mm = max(120, nrow(bidirectional_data) * 8)
      )
      
      cat(sprintf("    ✓ 双向关联图已保存 (%d 个双向显著关联)\n", nrow(bidirectional_data)))
} else {
      cat("    ⚠ 无双向显著关联，跳过双向关联图\n")
    }
  }, error = function(e) {
    cat(sprintf("    ⚠ 双向关联图生成失败: %s\n", e$message))
  })
  
  # 图表4：效应大小分布对比（小提琴图）
  cat("  生成图表4：效应大小分布对比...\n")
  tryCatch({
    effect_dist_data <- comparison_data %>%
      filter(!is.na(forward_or), !is.na(reverse_or)) %>%
      select(forward_log_or = forward_or, reverse_log_or = reverse_or,
             forward_pval, reverse_pval) %>%
      mutate(
        forward_log_or = log(forward_log_or),
        reverse_log_or = log(reverse_log_or)
      ) %>%
      pivot_longer(
        cols = c(forward_log_or, reverse_log_or),
        names_to = "direction",
        values_to = "log_or"
      ) %>%
      mutate(
        direction_label = ifelse(direction == "forward_log_or", "Forward MR\n(Metabolic→Lung Cancer)", 
                                 "Reverse MR\n(Lung Cancer→Metabolic)")
      )
    
    if (nrow(effect_dist_data) > 0) {
      p_effect_dist <- ggplot(effect_dist_data, aes(x = direction_label, y = log_or, fill = direction_label)) +
        geom_violin(alpha = 0.6, trim = FALSE, linewidth = 0.21) +
        geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.21) +
        scale_fill_manual(
          name = "MR Direction",
          values = c(
            "Forward MR\n(Metabolic→Lung Cancer)" = okabe_ito$orange,
            "Reverse MR\n(Lung Cancer→Metabolic)" = okabe_ito$sky_blue
          ),
          guide = "none"
        ) +
        labs(
          title = "Figure 6. Effect Size Distribution Comparison",
          subtitle = sprintf("Distribution of log(OR) values in forward vs reverse MR analyses\n(n = %d pairs)", 
                           nrow(comparison_data %>% filter(!is.na(forward_or), !is.na(reverse_or)))),
          x = "MR Direction",
          y = "log(OR)",
          caption = "Violin plots show the distribution of effect sizes. Boxplots show median and quartiles."
        ) +
        theme_sci() +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      save_sci_figure(
        p_effect_dist,
        "results/figures/reverse_mr/comparison/Figure6_Effect_Distribution",
        width_mm = 174,
        height_mm = 120
      )
      
      cat("    ✓ 效应分布对比图已保存\n")
    } else {
      cat("    ⚠ 无数据可用于效应分布图\n")
    }
  }, error = function(e) {
    cat(sprintf("    ⚠ 效应分布图生成失败: %s\n", e$message))
  })
  
  # 保存对比数据
  write.csv(comparison_data, 
           "results/tables/step11_forward_reverse_comparison.csv", 
           row.names = FALSE)
  cat("\n✓ 对比数据已保存: results/tables/step11_forward_reverse_comparison.csv\n")
  
  cat("\n【步骤11.4】生成对比分析报告\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  capture.output({
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("正向vs反向MR对比分析报告\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    
    cat("【数据概览】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("总配对数量:        %d 个\n", nrow(comparison_data)))
    cat(sprintf("有正向MR数据:      %d 个\n", sum(!is.na(comparison_data$forward_or))))
    cat(sprintf("有反向MR数据:      %d 个\n", sum(!is.na(comparison_data$reverse_or))))
    cat(sprintf("双向均有数据:      %d 个\n\n", 
               sum(!is.na(comparison_data$forward_or) & !is.na(comparison_data$reverse_or))))
    
    cat("【显著性模式分析】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("双向显著:          %d 个 (%.1f%%)\n", 
               sum(comparison_data$both_significant, na.rm = TRUE),
               100 * sum(comparison_data$both_significant, na.rm = TRUE) / nrow(comparison_data)))
    cat(sprintf("仅正向显著:        %d 个 (%.1f%%)\n", 
               sum(comparison_data$forward_only_sig, na.rm = TRUE),
               100 * sum(comparison_data$forward_only_sig, na.rm = TRUE) / nrow(comparison_data)))
    cat(sprintf("仅反向显著:        %d 个 (%.1f%%)\n", 
               sum(comparison_data$reverse_only_sig, na.rm = TRUE),
               100 * sum(comparison_data$reverse_only_sig, na.rm = TRUE) / nrow(comparison_data)))
    cat(sprintf("均不显著:          %d 个 (%.1f%%)\n\n", 
               sum(comparison_data$neither_sig, na.rm = TRUE),
               100 * sum(comparison_data$neither_sig, na.rm = TRUE) / nrow(comparison_data)))
    
    if (sum(comparison_data$both_significant, na.rm = TRUE) > 0) {
      cat("【双向显著关联详情】\n")
      cat(paste(rep("-", 80), collapse = ""), "\n")
      bidirectional <- comparison_data %>%
        filter(both_significant) %>%
        arrange(forward_pval) %>%
        head(10)
      
      for (i in seq_len(nrow(bidirectional))) {
        cat(sprintf("%d. %s ↔ %s\n", i, 
                   bidirectional$metabolic_trait[i], 
                   bidirectional$lung_cancer[i]))
        cat(sprintf("   正向MR: OR=%.3f (95%%CI: %.3f-%.3f), P=%.2e\n",
                   bidirectional$forward_or[i],
                   bidirectional$forward_or_lci[i],
                   bidirectional$forward_or_uci[i],
                   bidirectional$forward_pval[i]))
        cat(sprintf("   反向MR: OR=%.3f (95%%CI: %.3f-%.3f), P=%.2e\n\n",
                   bidirectional$reverse_or[i],
                   bidirectional$reverse_or_lci[i],
                   bidirectional$reverse_or_uci[i],
                   bidirectional$reverse_pval[i]))
      }
    }
    
    cat("\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("报告生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
  }, file = "results/step11_comparison_report.txt")
  
  cat("✓ 对比分析报告已保存: results/step11_comparison_report.txt\n\n")
  
  cat("【对比图表文件】\n")
  cat("  - results/figures/reverse_mr/comparison/Figure3_Forward_vs_Reverse_Scatter.png/pdf\n")
  cat("  - results/figures/reverse_mr/comparison/Figure4_Significance_Comparison.png/pdf\n")
  cat("  - results/figures/reverse_mr/comparison/Figure5_Bidirectional_Associations.png/pdf\n")
  cat("  - results/figures/reverse_mr/comparison/Figure6_Effect_Distribution.png/pdf\n\n")
  
  cat("【图表特点】\n")
  cat("  ✓ 符合SCI期刊10分标准（174mm双栏宽度）\n")
  cat("  ✓ 使用Okabe-Ito色盲友好调色板\n")
  cat("  ✓ 高分辨率输出（600 DPI PNG + 矢量PDF）\n")
  cat("  ✓ Arial字体，标准字号（9-11 pt）\n")
  cat("  ✓ 完整的图例和标注\n\n")
}
