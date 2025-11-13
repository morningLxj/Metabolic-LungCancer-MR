############################################################################
# 第6步：多变量孟德尔随机化分析与可视化（完整修复版）
# 功能：评估暴露因子对肺癌的独立因果效应并生成出版级图表
# 修复内容：
#   1. 添加category字段生成逻辑
#   2. 改进错误处理和数据验证
#   3. 完善标签映射系统
#   4. 优化MVMR分析流程
#   5. 修复热图生成问题
############################################################################

cat("第6步：开始多变量孟德尔随机化分析与可视化（完整修复版）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 1. 加载必要的包
cat("【步骤1/13】加载必要的R包\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

required_packages <- c("TwoSampleMR", "MVMR", "dplyr", "tidyr", "readr", 
                       "ggplot2", "stringr", "gridExtra", "grid", "scales",
                       "reshape2", "pheatmap", "openxlsx", "RColorBrewer")

# 可选包（用于备用分析方法）
optional_packages <- list(
  "MR.RAPS" = c("MR.RAPS", FALSE),  # MR-RAPS方法（通过TwoSampleMR已包含）
  "bmr" = c("bmr", TRUE)  # 贝叶斯多变量MR（需要从GitHub安装）
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装缺失的包: %s\n", pkg))
    if (pkg %in% c("TwoSampleMR", "MVMR")) {
      if (!require("remotes", quietly = TRUE)) install.packages("remotes")
      if (pkg == "TwoSampleMR") {
        remotes::install_github("MRCIEU/TwoSampleMR")
      } else {
        remotes::install_github("WSpiller/MVMR")
      }
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# 尝试加载可选包
optional_loaded <- list()
for (pkg_name in names(optional_packages)) {
  pkg_info <- optional_packages[[pkg_name]]
  pkg <- pkg_info[1]
  from_github <- as.logical(pkg_info[2])
  
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    optional_loaded[[pkg_name]] <- TRUE
    cat(sprintf("✓ 可选包已加载: %s\n", pkg))
  } else {
    optional_loaded[[pkg_name]] <- FALSE
    if (from_github) {
      cat(sprintf("ℹ 可选包未安装: %s (可从GitHub安装)\n", pkg))
    } else {
      cat(sprintf("ℹ 可选包未安装: %s\n", pkg))
    }
  }
}

cat("✓ 所有必需的包已加载\n")
if (any(unlist(optional_loaded))) {
  cat("✓ 部分可选包已加载，备用分析方法可用\n")
}
cat("\n")

# 2. 创建输出目录
cat("【步骤2/13】创建输出目录\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

dirs <- c("results/tables", "results/tables/paper_tables", 
          "results/figures", "results/figures/main_figures",
          "results/figures/mvmr", "results/figures/smoking_adjusted", 
          "data/processed/mvmr", "data/instruments", "data/outcome")
for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}
cat("✓ 输出目录已创建\n\n")

# 3. 加载前置数据
cat("【步骤3/12】加载前置数据\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 加载第2步数据（GWAS数据集定义）
if (file.exists("data/step02_gwas_datasets_latest.RData")) {
  load("data/step02_gwas_datasets_latest.RData")
  cat("✓ 已加载第2步数据（GWAS数据集）\n")
} else {
  cat("⚠ 未找到第2步数据，将使用内置定义\n")
  # 定义默认数据集
  outcomes <- list(
    lung_cancer_overall = "ebi-a-GCST004748",
    lung_adenocarcinoma = "ieu-a-984",
    squamous_cell_lung = "ieu-a-989"
  )
  
  metabolic_traits <- list(
    BMI = "ieu-b-40", HDL_cholesterol = "ieu-b-109", LDL_cholesterol = "ieu-b-110",
    triglycerides = "ieu-b-111", fasting_glucose = "ebi-a-GCST90002232",
    fasting_insulin = "ebi-a-GCST90002238", HbA1c = "ebi-a-GCST90014006",
    SBP = "ieu-b-38", DBP = "ieu-b-39", smoking_initiation = "ieu-b-4877",
    alcohol_drinks = "ieu-b-73", vitamin_D = "ebi-a-GCST90000618",
    circulating_leptin = "ebi-a-GCST90007316", ApoB = "ebi-a-GCST90025952",
    ApoA1 = "ebi-a-GCST90025955", IGF1 = "ebi-a-GCST90025989"
  )
  
  inflammatory_traits <- list(
    CRP = "ebi-a-GCST90029070", WBC = "ieu-b-30", IL6 = "ebi-a-GCST90012005",
    IL6R = "ebi-a-GCST90012025", TNFR1 = "ebi-a-GCST90012015"
  )
}

# 加载第4步准备数据
if (file.exists("data/step04_step6_prep_latest.RData")) {
  load("data/step04_step6_prep_latest.RData")
  cat("✓ 已加载第4步准备数据\n")
} else {
  cat("ℹ 未找到第4步准备数据，将使用第2步数据\n")
  step6_prep <- list(
    metabolic_traits = metabolic_traits,
    inflammatory_traits = inflammatory_traits,
    outcomes = outcomes
  )
}
cat("\n")

# 4. 检查并加载单变量MR结果
cat("【步骤4/12】检查单变量MR结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

univariate_results <- NULL
significant_exposures <- NULL

if (file.exists("results/tables/step05_mr_results_summary.csv")) {
  univariate_results <- read_csv("results/tables/step05_mr_results_summary.csv", 
                                 show_col_types = FALSE)
  cat(sprintf("✓ 已加载%d个单变量MR结果\n", nrow(univariate_results)))
  
  # 筛选显著暴露因子
  significant_exposures <- univariate_results %>%
    filter(pval < 0.05) %>%
    pull(exposure) %>%
    unique()
  
  cat(sprintf("✓ 发现%d个名义显著暴露因子（P<0.05）\n", length(significant_exposures)))
  
  # 添加category字段（如果不存在）
  if (!"category" %in% colnames(univariate_results)) {
    cat("  添加category字段...\n")
    univariate_results <- univariate_results %>%
      mutate(
        category = case_when(
          exposure %in% names(step6_prep$metabolic_traits) ~ "Metabolic",
          exposure %in% names(step6_prep$inflammatory_traits) ~ "Inflammatory",
          TRUE ~ "Other"
        )
      )
  }
  
  # 添加exposure_order字段
  if (!"exposure_order" %in% colnames(univariate_results)) {
    univariate_results <- univariate_results %>%
      mutate(
        exposure_order = case_when(
          !is.na(fdr_pval) & fdr_pval < 0.05 ~ 1,
          pval < 0.05 ~ 2,
          TRUE ~ 3
        )
      )
  }
  
} else {
  warning("⚠ 未找到单变量MR结果文件！请先运行第5步分析。")
  cat("  继续使用所有可用的暴露因子进行MVMR分析...\n")
}
cat("\n")

# 5. 定义MVMR分析组合
cat("【步骤5/12】定义MVMR分析组合\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 使用第4步准备的数据结构
metabolic_traits <- step6_prep$metabolic_traits
inflammatory_traits <- step6_prep$inflammatory_traits
outcomes <- step6_prep$outcomes

cat(sprintf("可用暴露因子：\n"))
cat(sprintf("  代谢性暴露: %d个\n", length(metabolic_traits)))
cat(sprintf("  炎症性暴露: %d个\n", length(inflammatory_traits)))
cat(sprintf("  结局变量: %d个\n", length(outcomes)))

# 精简的MVMR组合（每组2-3个暴露因子）
mvmr_combinations_reduced <- list(
  # 代谢组合
  metabolic_bmi_hdl = c("BMI", "HDL_cholesterol"),
  metabolic_bmi_ldl = c("BMI", "LDL_cholesterol"),
  metabolic_lipids = c("HDL_cholesterol", "LDL_cholesterol"),
  metabolic_glucose = c("fasting_glucose", "HbA1c"),
  
  # 炎症组合
  inflammatory_crp_il6r = c("CRP", "IL6R"),
  inflammatory_crp_wbc = c("CRP", "WBC"),
  
  # 生活方式组合
  lifestyle_smoking_alcohol = c("smoking_initiation", "alcohol_drinks"),
  lifestyle_smoking_bmi = c("smoking_initiation", "BMI"),
  
  # 混合组合
  mixed_bmi_crp = c("BMI", "CRP"),
  mixed_hdl_crp = c("HDL_cholesterol", "CRP")
)

# 吸烟调整组合（包含吸烟作为协变量）
smoking_adjusted_combinations <- list(
  smoking_adj_metabolic = c("smoking_initiation", "BMI", "HDL_cholesterol"),
  smoking_adj_inflammatory = c("smoking_initiation", "CRP", "IL6R"),
  smoking_adj_lipids = c("smoking_initiation", "HDL_cholesterol", "LDL_cholesterol")
)

# 筛选实际可用的组合
available_combinations <- list()
smoking_adjusted_available <- list()

cat("\n可用的基础MVMR组合：\n")
for (combo_name in names(mvmr_combinations_reduced)) {
  exposures <- mvmr_combinations_reduced[[combo_name]]
  available_exposures <- c()
  
  # 检查哪些暴露因子实际可用
  for (exp in exposures) {
    if (exp %in% names(metabolic_traits) || exp %in% names(inflammatory_traits)) {
      available_exposures <- c(available_exposures, exp)
    }
  }
  
  # 如果要求显著暴露，进一步筛选
  if (!is.null(significant_exposures)) {
    sig_available <- intersect(available_exposures, significant_exposures)
    if (length(sig_available) >= 2) {
      available_exposures <- sig_available
    }
  }
  
  # 至少需要2个暴露因子
  if (length(available_exposures) >= 2) {
    available_combinations[[combo_name]] <- available_exposures
    cat(sprintf("  ✓ %s: %s\n", combo_name, paste(available_exposures, collapse=" + ")))
  }
}

cat("\n吸烟调整的MVMR组合：\n")
for (combo_name in names(smoking_adjusted_combinations)) {
  exposures <- smoking_adjusted_combinations[[combo_name]]
  available_exposures <- c()
  
  for (exp in exposures) {
    if (exp %in% names(metabolic_traits) || exp %in% names(inflammatory_traits)) {
      available_exposures <- c(available_exposures, exp)
    }
  }
  
  if (length(available_exposures) >= 2) {
    smoking_adjusted_available[[combo_name]] <- available_exposures
    cat(sprintf("  ✓ %s: %s\n", combo_name, paste(available_exposures, collapse=" + ")))
  }
}

# 如果没有足够的组合，创建默认组合
if (length(available_combinations) == 0) {
  cat("\n⚠ 警告：没有找到足够的暴露因子组合，创建默认组合\n")
  default_exposures <- c(
    head(names(metabolic_traits), 2),
    head(names(inflammatory_traits), 1)
  )
  available_combinations[["default_combo"]] <- default_exposures
  cat(sprintf("  ✓ 默认组合: %s\n", paste(default_exposures, collapse=" + ")))
}

cat(sprintf("\n总计：%d个基础组合 + %d个吸烟调整组合\n\n", 
            length(available_combinations), 
            length(smoking_adjusted_available)))

# 6. 提取工具变量数据
cat("【步骤6/12】提取工具变量数据\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 清理旧缓存（可选）
cache_files <- list.files("data/instruments/", pattern = "_instruments.rds", full.names = TRUE)
if (length(cache_files) > 10) {
  cat("  清理旧的工具变量缓存...\n")
  file.remove(cache_files)
}

# 工具变量提取函数（改进版，使用多策略提取）
extract_exposure_instruments_safe <- function(exposure_id, exposure_name, use_cache = TRUE) {
  cache_file <- paste0("data/instruments/", exposure_name, "_instruments.rds")
  
  # 尝试从缓存加载
  if (use_cache && file.exists(cache_file)) {
    tryCatch({
      instruments <- readRDS(cache_file)
      if (!is.null(instruments) && is.data.frame(instruments) && nrow(instruments) >= 3) {
        # 验证缓存数据的完整性
        required_cols <- c("SNP", "beta.exposure", "se.exposure")
        if (all(required_cols %in% colnames(instruments))) {
        return(instruments)
        } else {
          cat(sprintf("    ⚠ %s: 缓存数据不完整，重新提取\n", exposure_name))
        }
      }
    }, error = function(e) {
      cat(sprintf("    ⚠ %s: 缓存读取失败: %s\n", exposure_name, conditionMessage(e)))
    })
  }
  
  # 使用多策略提取（参考step07的robust方式）
  strategies <- list(
    list(p1 = 5e-8, r2 = 0.001, kb = 10000, name = "严格"),
    list(p1 = 1e-6, r2 = 0.001, kb = 10000, name = "中等"),
    list(p1 = 5e-6, r2 = 0.01, kb = 5000, name = "宽松"),
    list(p1 = 1e-5, r2 = 0.05, kb = 5000, name = "极宽松")
  )
  
  for (strategy in strategies) {
    tryCatch({
      instruments <- TwoSampleMR::extract_instruments(
        outcomes = exposure_id,
        p1 = strategy$p1,
        clump = TRUE,
        r2 = strategy$r2,
        kb = strategy$kb
      )
      
      # 检查返回值
      if (!is.null(instruments) && is.data.frame(instruments) && nrow(instruments) >= 3) {
        # 验证必要列（至少需要基本列）
        required_cols <- c("SNP", "beta.exposure", "se.exposure")
      if (all(required_cols %in% colnames(instruments))) {
          # 如果缺少allele信息，尝试从其他列获取或标记为NA
          if (!"effect_allele" %in% colnames(instruments)) {
            if ("ea" %in% colnames(instruments)) {
              instruments$effect_allele <- instruments$ea
            } else {
              instruments$effect_allele <- NA_character_
            }
          }
          if (!"other_allele" %in% colnames(instruments)) {
            if ("nea" %in% colnames(instruments)) {
              instruments$other_allele <- instruments$nea
            } else {
              instruments$other_allele <- NA_character_
            }
          }
          
        # 保存到缓存
          tryCatch({
        saveRDS(instruments, cache_file)
          }, error = function(e) {
            cat(sprintf("    ⚠ %s: 缓存保存失败，但不影响使用\n", exposure_name))
          })
          
        return(instruments)
        } else {
          missing_cols <- setdiff(required_cols, colnames(instruments))
          cat(sprintf("    ⚠ %s: 策略'%s'缺少必要列: %s\n", 
                     exposure_name, strategy$name, paste(missing_cols, collapse=", ")))
      }
      } else if (!is.null(instruments)) {
        cat(sprintf("    ⚠ %s: 策略'%s' SNP不足 (%d个)\n", 
                   exposure_name, strategy$name, ifelse(is.data.frame(instruments), nrow(instruments), 0)))
    }
    
      # 批次之间添加延迟，避免API限流
      Sys.sleep(0.5)
    
  }, error = function(e) {
      error_msg <- conditionMessage(e)
      # 只在第一个策略失败时输出错误，避免过多输出
      if (strategy$name == "严格") {
        cat(sprintf("    ⚠ %s: 策略'%s'失败 - %s\n", 
                   exposure_name, strategy$name, 
                   ifelse(nchar(error_msg) > 100, substr(error_msg, 1, 100), error_msg)))
      }
    return(NULL)
  })
  }
  
  # 所有策略都失败
  cat(sprintf("    ✗ %s: 所有提取策略均失败\n", exposure_name))
  return(NULL)
}

# 提取所有需要的暴露因子的工具变量
exposure_instruments <- list()
all_combinations <- c(available_combinations, smoking_adjusted_available)

# 收集所有需要的暴露因子
all_exposures_needed <- unique(unlist(all_combinations))

cat(sprintf("需要提取工具变量的暴露因子：%d个\n\n", length(all_exposures_needed)))

for (exp_name in all_exposures_needed) {
  exp_id <- NULL
  
  # 查找暴露因子ID
  if (exp_name %in% names(metabolic_traits)) {
    exp_id <- metabolic_traits[[exp_name]]
  } else if (exp_name %in% names(inflammatory_traits)) {
    exp_id <- inflammatory_traits[[exp_name]]
  }
  
  if (!is.null(exp_id)) {
    cat(sprintf("  提取 %s (ID: %s)...", exp_name, exp_id))
    flush.console()  # 立即输出
    
    instruments <- extract_exposure_instruments_safe(exp_id, exp_name)
    
    if (!is.null(instruments) && is.data.frame(instruments) && nrow(instruments) >= 3) {
      exposure_instruments[[exp_name]] <- instruments
      cat(sprintf(" ✓ %d个SNP\n", nrow(instruments)))
    } else {
      # 提供更详细的失败信息
      if (is.null(instruments)) {
        cat(sprintf(" ✗ 失败（返回NULL）\n"))
      } else if (!is.data.frame(instruments)) {
        cat(sprintf(" ✗ 失败（返回类型错误: %s）\n", class(instruments)))
      } else if (nrow(instruments) < 3) {
        cat(sprintf(" ✗ 失败（SNP数量不足: %d个）\n", nrow(instruments)))
      } else {
        cat(sprintf(" ✗ 失败（未知原因）\n"))
      }
    }
  } else {
    cat(sprintf("  ⚠ %s: 未找到GWAS ID\n", exp_name))
  }
  
  # 添加短暂延迟，避免API限流
  Sys.sleep(0.3)
}

cat(sprintf("\n✓ 成功提取%d个暴露因子的工具变量（共需%d个）\n\n", 
            length(exposure_instruments), 
            length(all_exposures_needed)))

# 如果提取失败，尝试从之前的步骤加载已有的工具变量
if (length(exposure_instruments) == 0) {
  cat("\n⚠ 警告：无法提取工具变量，尝试从之前的步骤加载...\n")
  
  # 尝试从第4步加载
  if (file.exists("data/step04_all_instruments.RData")) {
    tryCatch({
      load("data/step04_all_instruments.RData")
      if (exists("instruments_results")) {
        if (!is.null(instruments_results$metabolic)) {
          for (name in names(instruments_results$metabolic)) {
            if (name %in% all_exposures_needed && !name %in% names(exposure_instruments)) {
              exposure_instruments[[name]] <- instruments_results$metabolic[[name]]
              cat(sprintf("  ✓ 从第4步加载: %s\n", name))
            }
          }
        }
        if (!is.null(instruments_results$inflammatory)) {
          for (name in names(instruments_results$inflammatory)) {
            if (name %in% all_exposures_needed && !name %in% names(exposure_instruments)) {
              exposure_instruments[[name]] <- instruments_results$inflammatory[[name]]
              cat(sprintf("  ✓ 从第4步加载: %s\n", name))
            }
          }
        }
      }
    }, error = function(e) {
      cat(sprintf("  ✗ 从第4步加载失败: %s\n", conditionMessage(e)))
    })
  }
  
  # 尝试从缓存目录加载
  for (exp_name in all_exposures_needed) {
    if (!exp_name %in% names(exposure_instruments)) {
      cache_file <- paste0("data/instruments/", exp_name, "_instruments.rds")
      if (file.exists(cache_file)) {
        tryCatch({
          cached_instruments <- readRDS(cache_file)
          if (!is.null(cached_instruments) && is.data.frame(cached_instruments) && nrow(cached_instruments) >= 3) {
            exposure_instruments[[exp_name]] <- cached_instruments
            cat(sprintf("  ✓ 从缓存加载: %s\n", exp_name))
          }
        }, error = function(e) {
          # 静默失败
        })
      }
    }
  }
  
  cat(sprintf("\n总计：从各种来源加载了%d个暴露因子的工具变量\n\n", length(exposure_instruments)))
  
  if (length(exposure_instruments) == 0) {
    stop(paste("\n错误：无法提取或加载任何工具变量，无法继续MVMR分析",
               "建议：",
               "  1. 检查网络连接",
               "  2. 确认GWAS ID是否正确",
               "  3. 先运行第4步提取工具变量",
               "  4. 检查data/instruments/目录是否有缓存文件",
               sep = "\n"))
  }
}

# 7. 提取结局数据
cat("【步骤7/12】提取结局数据\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

outcome_data_list <- list()

for (outcome_name in names(outcomes)) {
  outcome_id <- outcomes[[outcome_name]]
  
  cat(sprintf("  提取 %s...", outcome_name))
  
  tryCatch({
    # 检查缓存
    cache_file <- paste0("data/outcome/", outcome_name, "_outcome.rds")
    if (file.exists(cache_file)) {
      outcome_data <- readRDS(cache_file)
      outcome_data_list[[outcome_name]] <- outcome_data
      cat(sprintf(" ✓ 从缓存加载 (%d个SNP)\n", nrow(outcome_data)))
      next
    }
    
    # 收集所有需要的SNP
    all_snps <- unique(unlist(lapply(exposure_instruments, function(x) x$SNP)))
    
    if (length(all_snps) > 0) {
      outcome_data <- TwoSampleMR::extract_outcome_data(
        snps = all_snps,
        outcomes = outcome_id
      )
      
      if (!is.null(outcome_data) && nrow(outcome_data) > 0) {
        saveRDS(outcome_data, cache_file)
        outcome_data_list[[outcome_name]] <- outcome_data
        cat(sprintf(" ✓ %d个SNP\n", nrow(outcome_data)))
      } else {
        cat(sprintf(" ✗ 无数据\n"))
      }
    } else {
      cat(sprintf(" ✗ 无可用SNP\n"))
    }
    
  }, error = function(e) {
    cat(sprintf(" ✗ 失败: %s\n", conditionMessage(e)))
  })
}

if (length(outcome_data_list) == 0) {
  stop("错误：无法提取任何结局数据，无法继续MVMR分析")
}

cat(sprintf("\n✓ 成功提取%d个结局的数据\n\n", length(outcome_data_list)))

# 8. 定义MVMR分析函数
cat("【步骤8/12】定义MVMR分析函数\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 8.1 弱工具变量检测函数
check_weak_instruments <- function(harmonized_data_list, exposure_list, 
                                   weak_threshold = 10) {
  # 检测哪些暴露因子有弱工具变量问题
  weak_exposures <- list()
  
  for (exp_name in exposure_list) {
    if (exp_name %in% names(harmonized_data_list)) {
      data <- harmonized_data_list[[exp_name]]
      
      # 计算F统计量
      if ("beta.exposure" %in% names(data) && "se.exposure" %in% names(data)) {
        data$f_stat <- (data$beta.exposure^2) / (data$se.exposure^2)
        mean_f <- mean(data$f_stat, na.rm = TRUE)
        median_f <- median(data$f_stat, na.rm = TRUE)
        min_f <- min(data$f_stat, na.rm = TRUE)
        
        # 判断是否为弱工具变量（特别是IL-6等小样本性状）
        is_small_sample <- grepl("IL[_-]?6|IL6|interleukin.*6", exp_name, ignore.case = TRUE)
        threshold_adjusted <- ifelse(is_small_sample, weak_threshold * 0.8, weak_threshold)
        
        is_weak <- (mean_f < threshold_adjusted) || (median_f < threshold_adjusted)
        
        weak_exposures[[exp_name]] <- list(
          is_weak = is_weak,
          mean_f = mean_f,
          median_f = median_f,
          min_f = min_f,
          n_snps = nrow(data),
          is_small_sample = is_small_sample
        )
      }
    }
  }
  
  return(weak_exposures)
}

# 8.2 增强的敏感性分析函数
perform_enhanced_sensitivity <- function(harmonized_data_list, outcome_name, 
                                        exposure_list, common_snps) {
  sensitivity_results <- list()
  
  # 1. 异质性检验（Cochran's Q）
  tryCatch({
    for (exp_name in names(harmonized_data_list)) {
      data <- harmonized_data_list[[exp_name]] %>%
        filter(SNP %in% common_snps)
      
      if (nrow(data) >= 3) {
        # 单变量MR异质性检验
        het <- mr_heterogeneity(data)
        if (!is.null(het)) {
          sensitivity_results[[exp_name]]$heterogeneity <- het
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("  异质性检验失败: %s\n", conditionMessage(e)))
  })
  
  # 2. 水平多效性检验（MR-Egger截距）
  tryCatch({
    for (exp_name in names(harmonized_data_list)) {
      data <- harmonized_data_list[[exp_name]] %>%
        filter(SNP %in% common_snps)
      
      if (nrow(data) >= 3) {
        pleio <- mr_pleiotropy_test(data)
        if (!is.null(pleio)) {
          sensitivity_results[[exp_name]]$pleiotropy <- pleio
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("  多效性检验失败: %s\n", conditionMessage(e)))
  })
  
  # 3. 留一法分析（Leave-One-Out）
  tryCatch({
    for (exp_name in names(harmonized_data_list)) {
      data <- harmonized_data_list[[exp_name]] %>%
        filter(SNP %in% common_snps)
      
      if (nrow(data) >= 4) {
        loo <- mr_leaveoneout(data)
        if (!is.null(loo)) {
          sensitivity_results[[exp_name]]$leave_one_out <- loo
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("  留一法分析失败: %s\n", conditionMessage(e)))
  })
  
  # 4. 单SNP分析
  tryCatch({
    for (exp_name in names(harmonized_data_list)) {
      data <- harmonized_data_list[[exp_name]] %>%
        filter(SNP %in% common_snps)
      
      if (nrow(data) >= 2) {
        single_snp <- mr_singlesnp(data)
        if (!is.null(single_snp)) {
          sensitivity_results[[exp_name]]$single_snp <- single_snp
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("  单SNP分析失败: %s\n", conditionMessage(e)))
  })
  
  return(sensitivity_results)
}

# 辅助函数：重新提取暴露数据（使用更宽松的参数）
re_extract_exposure_data <- function(exp_name, exp_id, use_loose_strategy = TRUE) {
  tryCatch({
    if (use_loose_strategy) {
      # 使用非常宽松的策略：更大的p值、更小的r2、更大的kb
      instruments <- TwoSampleMR::extract_instruments(
        outcomes = exp_id,
        p1 = 5e-6,  # 更宽松的p值
        clump = TRUE,
        r2 = 0.001,  # 更小的r2（允许更紧密的连锁）
        kb = 10000   # 更大的kb（允许更大的物理距离）
      )
    } else {
      # 使用标准宽松策略
      instruments <- TwoSampleMR::extract_instruments(
        outcomes = exp_id,
        p1 = 5e-8,
        clump = TRUE,
        r2 = 0.001,
        kb = 10000
      )
    }
    
    if (!is.null(instruments) && is.data.frame(instruments) && nrow(instruments) >= 3) {
      required_cols <- c("SNP", "beta.exposure", "se.exposure")
      if (all(required_cols %in% colnames(instruments))) {
        # 补充allele信息
        if (!"effect_allele" %in% colnames(instruments)) {
          if ("ea" %in% colnames(instruments)) {
            instruments$effect_allele <- instruments$ea
          } else {
            instruments$effect_allele <- NA_character_
          }
        }
        if (!"other_allele" %in% colnames(instruments)) {
          if ("nea" %in% colnames(instruments)) {
            instruments$other_allele <- instruments$nea
          } else {
            instruments$other_allele <- NA_character_
          }
        }
        return(instruments)
      }
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

run_mvmr_analysis_robust <- function(exposure_list, outcome_name, 
                                     exposure_instruments, outcome_data,
                                     analysis_type = "standard",
                                     retry_with_re_extraction = TRUE) {
  
  analysis_id <- paste(paste(exposure_list, collapse="_"), outcome_name, sep="_to_")
  
  tryCatch({
    # 步骤1：收集暴露数据
    exp_data_list <- list()
    valid_exposures <- c()
    exposure_ids <- list()  # 保存ID用于重新提取
    
    for (exp_name in exposure_list) {
      if (exp_name %in% names(exposure_instruments)) {
        exp_data <- exposure_instruments[[exp_name]]
        
        required_cols <- c("SNP", "beta.exposure", "se.exposure")
        optional_cols <- c("effect_allele", "other_allele")
        
        # 检查必要列
        if (all(required_cols %in% colnames(exp_data))) {
          # 补充可选列
          if (!"effect_allele" %in% colnames(exp_data)) {
            if ("ea" %in% colnames(exp_data)) {
              exp_data$effect_allele <- exp_data$ea
            } else {
              exp_data$effect_allele <- NA_character_
            }
          }
          if (!"other_allele" %in% colnames(exp_data)) {
            if ("nea" %in% colnames(exp_data)) {
              exp_data$other_allele <- exp_data$nea
            } else {
              exp_data$other_allele <- NA_character_
            }
          }
          
          exp_data_list[[exp_name]] <- exp_data
          valid_exposures <- c(valid_exposures, exp_name)
          
          # 保存ID用于可能的重新提取
          if (exp_name %in% names(metabolic_traits)) {
            exposure_ids[[exp_name]] <- metabolic_traits[[exp_name]]
          } else if (exp_name %in% names(inflammatory_traits)) {
            exposure_ids[[exp_name]] <- inflammatory_traits[[exp_name]]
          }
        }
      }
    }
    
    if (length(valid_exposures) < 2) {
      return(list(error = sprintf("有效暴露因子不足: %d/2", length(valid_exposures))))
    }
    
    # 步骤2：数据协调
    harmonized_data <- list()
    
    for (exp_name in names(exp_data_list)) {
      harmonized <- TwoSampleMR::harmonise_data(
        exposure_dat = exp_data_list[[exp_name]],
        outcome_dat = outcome_data,
        action = 2  # 保留回文SNP
      )
      
      if (!is.null(harmonized) && nrow(harmonized) > 0) {
        # 移除不兼容的SNP
        harmonized <- harmonized %>%
          filter(!is.na(beta.exposure), !is.na(se.exposure),
                 !is.na(beta.outcome), !is.na(se.outcome),
                 se.exposure > 0, se.outcome > 0)
        
        if (nrow(harmonized) > 0) {
          harmonized_data[[exp_name]] <- harmonized
        }
      }
    }
    
    if (length(harmonized_data) < 2) {
      return(list(error = "协调后的暴露因子不足2个"))
    }
    
    # 步骤3：找到共同的SNP
    all_snps <- lapply(harmonized_data, function(x) x$SNP)
    common_snps <- Reduce(intersect, all_snps)
    
    # 放宽要求：对于2个暴露因子，至少需要2个SNP（自由度 >= 1）
    # 对于3个或更多暴露因子，至少需要暴露因子数+1个SNP
    n_exp_actual <- length(harmonized_data)
    min_snps_required <- max(3, n_exp_actual + 1)  # 至少3个SNP以保证有自由度
    
    if (length(common_snps) < min_snps_required) {
      # 如果SNP不足且允许重试，尝试重新提取更完整的暴露数据
      if (retry_with_re_extraction && length(exposure_ids) == length(valid_exposures)) {
        # 尝试重新提取所有暴露因子的数据
        cat(sprintf("      [重试] 共同SNP不足(%d/%d)，尝试重新提取暴露数据...\n", 
                    length(common_snps), min_snps_required))
        
        re_extracted_count <- 0
        for (exp_name in valid_exposures) {
          if (exp_name %in% names(exposure_ids)) {
            exp_id <- exposure_ids[[exp_name]]
            new_data <- re_extract_exposure_data(exp_name, exp_id, use_loose_strategy = TRUE)
            if (!is.null(new_data) && nrow(new_data) > nrow(exp_data_list[[exp_name]])) {
              exp_data_list[[exp_name]] <- new_data
              re_extracted_count <- re_extracted_count + 1
            }
          }
          Sys.sleep(0.3)  # 避免API限流
        }
        
        if (re_extracted_count > 0) {
          # 重新进行协调
          harmonized_data <- list()
          for (exp_name in names(exp_data_list)) {
            harmonized <- TwoSampleMR::harmonise_data(
              exposure_dat = exp_data_list[[exp_name]],
              outcome_dat = outcome_data,
              action = 2
            )
            
            if (!is.null(harmonized) && nrow(harmonized) > 0) {
              harmonized <- harmonized %>%
                filter(!is.na(beta.exposure), !is.na(se.exposure),
                       !is.na(beta.outcome), !is.na(se.outcome),
                       se.exposure > 0, se.outcome > 0)
              
              if (nrow(harmonized) > 0) {
                harmonized_data[[exp_name]] <- harmonized
              }
            }
          }
          
          if (length(harmonized_data) >= 2) {
            # 重新计算共同SNP
            all_snps <- lapply(harmonized_data, function(x) x$SNP)
            common_snps <- Reduce(intersect, all_snps)
            n_exp_actual <- length(harmonized_data)
            min_snps_required <- max(3, n_exp_actual + 1)
            
            if (length(common_snps) >= min_snps_required) {
              cat(sprintf("      ✓ 重新提取后获得%d个共同SNP\n", length(common_snps)))
              # 继续后续分析
            } else {
              return(list(error = sprintf("共同SNP不足: %d/%d (重试后仍不足)", 
                                         length(common_snps), min_snps_required)))
            }
          } else {
            return(list(error = sprintf("共同SNP不足: %d/%d (重试后协调失败)", 
                                       length(common_snps), min_snps_required)))
          }
        } else {
          return(list(error = sprintf("共同SNP不足: %d/%d (重新提取未获得更多SNP)", 
                                     length(common_snps), min_snps_required)))
        }
      } else {
        # 不允许重试或没有ID信息，直接返回错误
        if (length(common_snps) == 0) {
          return(list(error = sprintf("共同SNP不足: %d/%d (完全无重叠)", 
                                     length(common_snps), min_snps_required)))
        } else {
          return(list(error = sprintf("共同SNP不足: %d/%d", 
                                     length(common_snps), min_snps_required)))
        }
      }
    }
    
    # 步骤4：构建矩阵
    valid_exposures_final <- names(harmonized_data)
    n_exp <- length(valid_exposures_final)
    n_snp <- length(common_snps)
    
    beta_matrix <- matrix(NA, nrow = n_snp, ncol = n_exp)
    se_matrix <- matrix(NA, nrow = n_snp, ncol = n_exp)
    colnames(beta_matrix) <- valid_exposures_final
    colnames(se_matrix) <- valid_exposures_final
    rownames(beta_matrix) <- common_snps
    rownames(se_matrix) <- common_snps
    
    for (i in seq_along(valid_exposures_final)) {
      exp_name <- valid_exposures_final[i]
      exp_data <- harmonized_data[[exp_name]] %>%
        filter(SNP %in% common_snps) %>%
        arrange(match(SNP, common_snps))
      
      beta_matrix[, i] <- exp_data$beta.exposure
      se_matrix[, i] <- exp_data$se.exposure
    }
    
    # 提取结局数据
    outcome_snps_data <- outcome_data %>%
      filter(SNP %in% common_snps) %>%
      arrange(match(SNP, common_snps))
    
    # 数据质量检查
    if (any(is.na(beta_matrix)) || any(is.na(se_matrix)) || 
        any(is.na(outcome_snps_data$beta.outcome)) || 
        any(is.na(outcome_snps_data$se.outcome))) {
      return(list(error = "数据矩阵包含缺失值"))
    }
    
    # 步骤4.5：检测弱工具变量
    weak_instrument_check <- check_weak_instruments(harmonized_data, valid_exposures_final)
    has_weak_instruments <- any(sapply(weak_instrument_check, function(x) x$is_weak))
    weak_exposure_names <- names(weak_instrument_check)[sapply(weak_instrument_check, function(x) x$is_weak)]
    
    # 步骤5：运行MVMR
    mvmr_input <- MVMR::format_mvmr(
      BXGs = beta_matrix,
      BYG = outcome_snps_data$beta.outcome,
      seBXGs = se_matrix,
      seBYG = outcome_snps_data$se.outcome,
      RSID = common_snps
    )
    
    mvmr_result <- MVMR::ivw_mvmr(mvmr_input)
    
    # 检查结果有效性 - ivw_mvmr可能返回data.frame或list
    if (is.null(mvmr_result)) {
      return(list(error = "MVMR分析返回NULL"))
    }
    
    # 调试：打印结果类型和结构（仅在调试模式下）
    # cat(sprintf("      [调试] MVMR结果类型: %s\n", class(mvmr_result)))
    # if (is.data.frame(mvmr_result)) {
    #   cat(sprintf("      [调试] 列名: %s\n", paste(colnames(mvmr_result), collapse=", ")))
    #   cat(sprintf("      [调试] 行数: %d\n", nrow(mvmr_result)))
    # } else if (is.list(mvmr_result)) {
    #   cat(sprintf("      [调试] 元素名: %s\n", paste(names(mvmr_result), collapse=", ")))
    # }
    
    # 尝试多种方式提取结果
    beta_vals <- NULL
    se_vals <- NULL
    pval_vals <- NULL
    
    # 方法1: 如果是matrix类型
    if (is.matrix(mvmr_result) || is.array(mvmr_result)) {
      # 转换为data.frame以便统一处理
      mvmr_result_df <- as.data.frame(mvmr_result)
      
      # 尝试不同的列名变体
      estimate_col <- NULL
      se_col <- NULL
      pval_col <- NULL
      
      # 查找Estimate列
      for (col in c("Estimate", "estimate", "beta", "Beta", "coef", "Coef")) {
        if (col %in% colnames(mvmr_result_df)) {
          estimate_col <- col
          break
        }
      }
      
      # 查找SE列
      for (col in c("Std. Error", "Std Error", "SE", "se", "se.exposure", "std.error")) {
        if (col %in% colnames(mvmr_result_df)) {
          se_col <- col
          break
        }
      }
      
      # 查找P值列
      for (col in c("Pr(>|t|)", "Pr...t..", "pval", "Pvalue", "p.value", "p")) {
        if (col %in% colnames(mvmr_result_df)) {
          pval_col <- col
          break
        }
      }
      
      if (!is.null(estimate_col) && !is.null(se_col) && !is.null(pval_col)) {
        beta_vals <- mvmr_result_df[[estimate_col]]
        se_vals <- mvmr_result_df[[se_col]]
        pval_vals <- mvmr_result_df[[pval_col]]
      } else if (ncol(mvmr_result) >= 4) {
        # 如果列名匹配失败，按位置提取（标准格式：Estimate, Std.Error, t.value, Pr>|t|）
        beta_vals <- mvmr_result[, 1]
        se_vals <- mvmr_result[, 2]
        pval_vals <- mvmr_result[, 4]
      }
      
    # 方法2: 如果是data.frame
    } else if (is.data.frame(mvmr_result)) {
      # 尝试不同的列名变体
      estimate_col <- NULL
      se_col <- NULL
      pval_col <- NULL
      
      # 查找Estimate列
      for (col in c("Estimate", "estimate", "beta", "Beta", "coef", "Coef")) {
        if (col %in% colnames(mvmr_result)) {
          estimate_col <- col
          break
        }
      }
      
      # 查找SE列
      for (col in c("Std. Error", "Std Error", "SE", "se", "se.exposure", "std.error")) {
        if (col %in% colnames(mvmr_result)) {
          se_col <- col
          break
        }
      }
      
      # 查找P值列
      for (col in c("Pr(>|t|)", "Pr...t..", "pval", "Pvalue", "p.value", "p")) {
        if (col %in% colnames(mvmr_result)) {
          pval_col <- col
          break
        }
      }
      
      if (!is.null(estimate_col) && !is.null(se_col) && !is.null(pval_col)) {
        beta_vals <- mvmr_result[[estimate_col]]
        se_vals <- mvmr_result[[se_col]]
        pval_vals <- mvmr_result[[pval_col]]
      }
      
    } else if (is.list(mvmr_result) || inherits(mvmr_result, "summary")) {
      # 方法2: 如果是list或summary对象
      # 尝试提取系数
      if ("coefficients" %in% names(mvmr_result)) {
        coef_table <- mvmr_result$coefficients
        if (is.matrix(coef_table) || is.data.frame(coef_table)) {
          beta_vals <- coef_table[, 1]  # 第一列通常是系数
          if (ncol(coef_table) >= 2) se_vals <- coef_table[, 2]  # 第二列是SE
          if (ncol(coef_table) >= 4) pval_vals <- coef_table[, 4]  # 第四列是P值
        }
      } else {
        # 尝试直接访问
        if ("Estimate" %in% names(mvmr_result)) beta_vals <- mvmr_result$Estimate
        if ("Std. Error" %in% names(mvmr_result)) se_vals <- mvmr_result$`Std. Error`
        if ("Pr(>|t|)" %in% names(mvmr_result)) pval_vals <- mvmr_result$`Pr(>|t|)`
      }
    }
    
    # 如果仍未提取到，尝试使用summary
    if (is.null(beta_vals)) {
      tryCatch({
        summary_result <- summary(mvmr_result)
        if (!is.null(summary_result)) {
          if (is.data.frame(summary_result)) {
            # summary返回data.frame
            beta_vals <- summary_result$Estimate
            se_vals <- summary_result$`Std. Error`
            pval_vals <- summary_result$`Pr(>|t|)`
          } else if (is.list(summary_result) && "coefficients" %in% names(summary_result)) {
            coef_table <- summary_result$coefficients
            beta_vals <- coef_table[, 1]
            if (ncol(coef_table) >= 2) se_vals <- coef_table[, 2]
            if (ncol(coef_table) >= 4) pval_vals <- coef_table[, 4]
          }
        }
      }, error = function(e) {
        # 忽略错误
      })
    }
    
    # 如果仍然无法提取结果，返回错误
    if (is.null(beta_vals) || is.null(se_vals) || is.null(pval_vals)) {
      result_info <- sprintf("类型: %s, ", paste(class(mvmr_result), collapse=", "))
      if (is.data.frame(mvmr_result)) {
        result_info <- paste0(result_info, sprintf("列名: %s", paste(colnames(mvmr_result), collapse=", ")))
      } else if (is.matrix(mvmr_result) || is.array(mvmr_result)) {
        result_info <- paste0(result_info, sprintf("维度: %s, 列名: %s", 
                                                   paste(dim(mvmr_result), collapse="x"),
                                                   paste(colnames(mvmr_result), collapse=", ")))
      } else if (is.list(mvmr_result)) {
        result_info <- paste0(result_info, sprintf("元素名: %s", paste(names(mvmr_result), collapse=", ")))
      }
      
      return(list(error = sprintf("MVMR结果无法提取 (beta: %s, se: %s, pval: %s). %s", 
                                 ifelse(is.null(beta_vals), "NULL", "OK"),
                                 ifelse(is.null(se_vals), "NULL", "OK"),
                                 ifelse(is.null(pval_vals), "NULL", "OK"),
                                 result_info)))
    }
    
    # 检查长度是否匹配
    if (length(beta_vals) != length(valid_exposures_final) ||
        length(se_vals) != length(valid_exposures_final) ||
        length(pval_vals) != length(valid_exposures_final)) {
      return(list(error = sprintf("MVMR结果长度不匹配: beta=%d, se=%d, pval=%d, 暴露因子数=%d", 
                                 length(beta_vals), length(se_vals), length(pval_vals),
                                 length(valid_exposures_final))))
    }
    
    # 检查是否有NA值
    if (any(is.na(beta_vals)) || any(is.na(se_vals)) || any(is.na(pval_vals))) {
      return(list(error = "MVMR结果包含NA值"))
    }
    
    # 步骤6：整理结果，添加弱工具变量信息
    # 为每个暴露因子添加弱工具变量检测信息
    weak_info_list <- lapply(valid_exposures_final, function(exp_name) {
      if (exp_name %in% names(weak_instrument_check)) {
        info <- weak_instrument_check[[exp_name]]
        list(
          is_weak = info$is_weak,
          mean_f = info$mean_f,
          median_f = info$median_f,
          min_f = info$min_f,
          n_snps_f = info$n_snps,
          is_small_sample = info$is_small_sample
        )
      } else {
        list(
          is_weak = FALSE,
          mean_f = NA,
          median_f = NA,
          min_f = NA,
          n_snps_f = NA,
          is_small_sample = FALSE
        )
      }
    })
    
    results_df <- data.frame(
      outcome = outcome_name,
      exposure = valid_exposures_final,
      beta = beta_vals,
      se = se_vals,
      pval = pval_vals,
      n_snps = n_snp,
      analysis_type = analysis_type,
      is_weak_instrument = sapply(weak_info_list, function(x) x$is_weak),
      mean_f_statistic = sapply(weak_info_list, function(x) x$mean_f),
      median_f_statistic = sapply(weak_info_list, function(x) x$median_f),
      min_f_statistic = sapply(weak_info_list, function(x) x$min_f),
      is_small_sample_trait = sapply(weak_info_list, function(x) x$is_small_sample),
      stringsAsFactors = FALSE
    )
    
    # 计算OR和置信区间
    results_df <- results_df %>%
      mutate(
        or = exp(beta),
        or_lci = exp(beta - 1.96 * se),
        or_uci = exp(beta + 1.96 * se),
        or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci)
      )
    
    # 步骤6.5：执行增强的敏感性分析
    tryCatch({
      sensitivity_results <- perform_enhanced_sensitivity(
        harmonized_data_list = harmonized_data,
        outcome_name = outcome_name,
        exposure_list = valid_exposures_final,
        common_snps = common_snps
      )
      
      # 保存敏感性分析结果（作为结果的属性）
      attr(results_df, "sensitivity") <- sensitivity_results
      
    }, error = function(e) {
      cat(sprintf("      ⚠ 敏感性分析失败: %s\n", conditionMessage(e)))
    })
    
    return(results_df)
    
  }, error = function(e) {
    return(list(error = conditionMessage(e)))
  })
}

cat("✓ MVMR分析函数已定义\n\n")

# 9. 运行MVMR分析
cat("【步骤9/12】运行MVMR分析\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

all_mvmr_results <- list()
smoking_adjusted_results <- list()
failed_analyses <- list()

# 运行基础组合分析
cat("\n▶ 运行基础MVMR组合分析:\n")
for (combo_name in names(available_combinations)) {
  exposure_list <- available_combinations[[combo_name]]
  
  cat(sprintf("\n  组合: %s\n", combo_name))
  cat(sprintf("    暴露因子: %s\n", paste(exposure_list, collapse=" + ")))
  
  for (outcome_name in names(outcome_data_list)) {
    cat(sprintf("    → %s... ", outcome_name))
    
    result <- run_mvmr_analysis_robust(
      exposure_list = exposure_list,
      outcome_name = outcome_name,
      exposure_instruments = exposure_instruments,
      outcome_data = outcome_data_list[[outcome_name]],
      analysis_type = "standard"
    )
    
    if (is.data.frame(result)) {
      result$combination_name <- combo_name
      all_mvmr_results[[length(all_mvmr_results) + 1]] <- result
      cat("✓\n")
    } else if (is.list(result) && "error" %in% names(result)) {
      failed_analyses[[length(failed_analyses) + 1]] <- list(
        combination = combo_name,
        outcome = outcome_name,
        reason = result$error
      )
      cat(sprintf("✗ (%s)\n", result$error))
    }
  }
}

# 运行吸烟调整组合分析
if (length(smoking_adjusted_available) > 0) {
  cat("\n▶ 运行吸烟调整的MVMR分析:\n")
  
  for (combo_name in names(smoking_adjusted_available)) {
    exposure_list <- smoking_adjusted_available[[combo_name]]
    
    cat(sprintf("\n  组合: %s\n", combo_name))
    cat(sprintf("    暴露因子: %s\n", paste(exposure_list, collapse=" + ")))
    
    for (outcome_name in names(outcome_data_list)) {
      cat(sprintf("    → %s... ", outcome_name))
      
      result <- run_mvmr_analysis_robust(
        exposure_list = exposure_list,
        outcome_name = outcome_name,
        exposure_instruments = exposure_instruments,
        outcome_data = outcome_data_list[[outcome_name]],
        analysis_type = "smoking_adjusted"
      )
      
      if (is.data.frame(result)) {
        result$combination_name <- combo_name
        smoking_adjusted_results[[length(smoking_adjusted_results) + 1]] <- result
        cat("✓\n")
      } else if (is.list(result) && "error" %in% names(result)) {
        failed_analyses[[length(failed_analyses) + 1]] <- list(
          combination = combo_name,
          outcome = outcome_name,
          reason = result$error
        )
        cat(sprintf("✗ (%s)\n", result$error))
      }
    }
  }
}

cat(sprintf("\n分析完成：\n"))
cat(sprintf("  成功的基础分析: %d\n", length(all_mvmr_results)))
cat(sprintf("  成功的吸烟调整分析: %d\n", length(smoking_adjusted_results)))
cat(sprintf("  失败的分析: %d\n\n", length(failed_analyses)))

# 10. 合并和保存MVMR结果
cat("【步骤10/12】合并和保存MVMR结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

all_results <- c(all_mvmr_results, smoking_adjusted_results)

if (length(all_results) > 0) {
  mvmr_results_df <- bind_rows(all_results)
  
  # 计算FDR校正
  mvmr_results_df <- mvmr_results_df %>%
    mutate(
      fdr = p.adjust(pval, method = "fdr"),
      significant_fdr = fdr < 0.05,
      significant_nominal = pval < 0.05
    ) %>%
    arrange(pval)
  
  # 保存详细结果
  write.csv(mvmr_results_df, "results/tables/step06_mvmr_detailed_final.csv", row.names = FALSE)
  write.csv(mvmr_results_df, "results/tables/step06_mvmr_results_final.csv", row.names = FALSE)
  cat(sprintf("✓ 已保存MVMR结果：%d个分析\n", nrow(mvmr_results_df)))
  
  # 单独保存吸烟调整结果
  if (length(smoking_adjusted_results) > 0) {
    smoking_df <- bind_rows(smoking_adjusted_results) %>%
      mutate(
        fdr = p.adjust(pval, method = "fdr"),
        significant_fdr = fdr < 0.05,
        significant_nominal = pval < 0.05
      ) %>%
      arrange(pval)
    
    write.csv(smoking_df, "results/tables/step06_smoking_adjusted_final.csv", row.names = FALSE)
    cat(sprintf("✓ 已保存吸烟调整结果：%d个分析\n", nrow(smoking_df)))
  }
  
  # 显著性统计
  n_significant_fdr <- sum(mvmr_results_df$significant_fdr, na.rm = TRUE)
  n_significant_nominal <- sum(mvmr_results_df$significant_nominal, na.rm = TRUE)
  
  cat(sprintf("\n显著结果统计：\n"))
  cat(sprintf("  FDR < 0.05: %d个 (%.1f%%)\n", 
              n_significant_fdr, 
              100 * n_significant_fdr / nrow(mvmr_results_df)))
  cat(sprintf("  P < 0.05: %d个 (%.1f%%)\n", 
              n_significant_nominal, 
              100 * n_significant_nominal / nrow(mvmr_results_df)))
  
} else {
  warning("⚠ 没有成功的MVMR分析结果")
  write.csv(data.frame(), "results/tables/step06_mvmr_results_final.csv", row.names = FALSE)
  write.csv(data.frame(), "results/tables/step06_mvmr_detailed_final.csv", row.names = FALSE)
}

# 保存失败分析记录
if (length(failed_analyses) > 0) {
  failed_df <- bind_rows(lapply(failed_analyses, as.data.frame))
  write.csv(failed_df, "results/tables/step06_failed_analyses.csv", row.names = FALSE)
  cat(sprintf("ℹ 已保存失败分析记录：%d个\n", length(failed_analyses)))
}

cat("\n")

# 11. 创建可视化（基于单变量MR结果）
cat("【步骤11/13】创建单变量MR可视化\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 检查是否有单变量MR结果用于可视化
if (is.null(univariate_results) || nrow(univariate_results) == 0) {
  cat("⚠ 没有单变量MR结果，跳过可视化步骤\n")
  cat("  提示：请先运行第5步分析生成单变量结果\n\n")
} else {
  
  # 准备可视化数据
  mr_summary <- univariate_results
  
  # 定义完整的标签映射
  exposure_labels <- c(
    "BMI" = "Body Mass Index (BMI)",
    "HbA1c" = "Glycated Hemoglobin (HbA1c)",
    "CRP" = "C-Reactive Protein (CRP)",
    "WBC" = "White Blood Cell Count (WBC)",
    "IL6" = "Interleukin-6 (IL-6)",
    "IL6R" = "IL-6 Receptor (IL-6R)",
    "TNFR1" = "TNF Receptor 1 (TNFR1)",
    "triglycerides" = "Triglycerides",
    "HDL_cholesterol" = "HDL Cholesterol",
    "LDL_cholesterol" = "LDL Cholesterol",
    "smoking_initiation" = "Smoking Initiation",
    "alcohol_drinks" = "Alcohol Drinks per Week",
    "fasting_glucose" = "Fasting Glucose",
    "fasting_insulin" = "Fasting Insulin",
    "SBP" = "Systolic Blood Pressure",
    "DBP" = "Diastolic Blood Pressure",
    "hypertension" = "Hypertension",
    "vitamin_D" = "Vitamin D",
    "circulating_leptin" = "Circulating Leptin",
    "ApoB" = "Apolipoprotein B",
    "ApoA1" = "Apolipoprotein A1",
    "ApoB_ApoA1_ratio" = "ApoB/ApoA1 Ratio",
    "IGF1" = "Insulin-like Growth Factor 1",
    "HDL_diameter" = "HDL Particle Diameter",
    "HDL_large" = "Large HDL Particles",
    "HDL_very_large" = "Very Large HDL Particles",
    "remnant_cholesterol" = "Remnant Cholesterol",
    "LDL_small" = "Small LDL Particles",
    "BCAA" = "Branched-Chain Amino Acids",
    "GGT" = "Gamma-Glutamyl Transferase"
  )
  
  outcome_labels <- c(
    "lung_adenocarcinoma" = "Lung Adenocarcinoma",
    "squamous_cell_lung" = "Squamous Cell Lung Cancer",
    "lung_cancer_overall" = "Overall Lung Cancer"
  )
  
  # 应用标签映射
  mr_summary <- mr_summary %>%
    mutate(
      exposure_label = ifelse(exposure %in% names(exposure_labels), 
                             exposure_labels[exposure], 
                             exposure),
      outcome_label = ifelse(outcome %in% names(outcome_labels), 
                            outcome_labels[outcome], 
                            outcome)
    )
  
  # 色盲友好调色板（Okabe-Ito配色）
  okabe_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                    "#0072B2", "#D55E00", "#CC79A7")
  metabolic_color <- okabe_colors[1]  # 橙色
  inflammatory_color <- okabe_colors[2]  # 蓝色
  
  # 森林图创建函数
  create_forest_plot <- function(data, title, subtitle = NULL) {
    plot_data <- data %>% 
      filter(!is.na(or), !is.na(or_lci), !is.na(or_uci))
    
    if (nrow(plot_data) == 0) {
      cat("    ⚠ 无有效数据，跳过该图\n")
      return(NULL)
    }
    
    # 排序
    plot_data <- plot_data %>%
      arrange(exposure_order, or) %>%
      mutate(exposure_label_factor = factor(exposure_label, levels = unique(exposure_label)))
    
    # 确定使用的颜色
    present_categories <- unique(plot_data$category)
    dynamic_colors <- c(
      "Metabolic" = metabolic_color,
      "Inflammatory" = inflammatory_color
    )[present_categories]
    
    # 创建图表
    p <- ggplot(plot_data, aes(x = or, y = exposure_label_factor)) +
      geom_errorbarh(aes(xmin = or_lci, xmax = or_uci, color = category), 
                     height = 0.3, size = 0.9, alpha = 0.7) +
      geom_point(aes(color = category), size = 4, alpha = 0.9) +
      geom_point(
        data = filter(plot_data, significant_fdr == TRUE), 
        aes(x = or, y = exposure_label_factor), 
        shape = 16, size = 2.5, color = "black"
      ) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
      scale_x_continuous(
        trans = "log10",
        breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2.0),
        labels = c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00"),
        name = "Odds Ratio (95% CI, log scale)"
      ) +
      scale_y_discrete(name = NULL) +
      scale_color_manual(
        name = "Exposure Category",
        values = dynamic_colors
      ) +
      labs(
        title = title,
        subtitle = subtitle,
        caption = "Black dots indicate FDR-adjusted significance (P < 0.05)."
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 13, color = "gray30"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 13, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    return(p)
  }
  
  # 方案A：合并分面图
  cat("\n生成方案A：合并分面森林图\n")
  
  facet_data <- mr_summary %>%
    group_by(exposure_label) %>%
    mutate(
      max_sig = min(exposure_order),
      max_or = max(or, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    arrange(max_sig, desc(max_or), exposure_label) %>%
    mutate(exposure_label_factor = factor(exposure_label, levels = unique(exposure_label)))
  
  p_facet <- ggplot(facet_data, aes(x = or, y = exposure_label_factor)) +
    geom_errorbarh(aes(xmin = or_lci, xmax = or_uci, color = category), 
                   height = 0.2, size = 0.7, alpha = 0.7) +
    geom_point(aes(color = category), size = 3, alpha = 0.9) +
    geom_point(
      data = filter(facet_data, significant_fdr == TRUE), 
      aes(x = or, y = exposure_label_factor), 
      shape = 16, size = 2, color = "black"
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_x_continuous(
      trans = "log10", 
      breaks = c(0.5, 0.75, 1, 1.5, 2), 
      name = "Odds Ratio (95% CI, log scale)"
    ) +
    scale_y_discrete(name = NULL) +
    scale_color_manual(
      name = "Category", 
      values = c("Metabolic" = metabolic_color, "Inflammatory" = inflammatory_color)
    ) +
    facet_grid(. ~ outcome_label, scales = "free_x", space = "free_x") +
    labs(
      title = "Causal Effects of Metabolic and Inflammatory Traits on Lung Cancer Subtypes",
      subtitle = "A Two-Sample Mendelian Randomization Analysis",
      caption = "Black dots indicate FDR-adjusted significance (P < 0.05)."
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )
  
  ggsave("results/figures/Figure1_Facet_ForestPlot.png", 
         plot = p_facet, 
         width = 14, 
         height = 12, 
         dpi = 300)
  cat("  ✓ 已保存: Figure1_Facet_ForestPlot.png\n")
  
  # 方案B：按结局分别生成独立图
  cat("\n生成方案B：按结局分别生成独立森林图\n")
  
  unique_outcomes <- unique(mr_summary$outcome_label)
  plot_counter <- 1
  
  for (outcome in unique_outcomes) {
    outcome_data <- mr_summary %>% filter(outcome_label == outcome)
    
    plot_title <- paste("Causal Effects on", outcome)
    plot_subtitle <- "Mendelian Randomization Analysis"
    
    p <- create_forest_plot(
      data = outcome_data,
      title = plot_title,
      subtitle = plot_subtitle
    )
    
    if (!is.null(p)) {
      filename <- paste0("results/figures/Figure1_", 
                        letters[plot_counter], 
                        "_ForestPlot_", 
                        gsub(" ", "_", outcome), 
                        ".png")
      
      plot_height <- max(8, nrow(outcome_data) * 0.35)
      ggsave(filename, plot = p, width = 12, height = plot_height, dpi = 300)
      cat(sprintf("  ✓ 已保存: %s\n", basename(filename)))
      plot_counter <- plot_counter + 1
    }
  }
  
  # 汇总统计图
  cat("\n生成汇总统计图\n")
  
  summary_stats <- mr_summary %>%
    group_by(category, outcome_label) %>%
    summarise(
      total_analyses = n(),
      significant_analyses = sum(significant_fdr, na.rm = TRUE),
      prop_significant = significant_analyses / total_analyses,
      .groups = "drop"
    )
  
  p_summary <- ggplot(summary_stats, 
                     aes(x = outcome_label, y = prop_significant, fill = category)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prop_significant * 100)), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, 
              size = 3.5, 
              fontface = "bold") +
    scale_y_continuous(
      labels = percent_format(accuracy = 1), 
      limits = c(0, max(summary_stats$prop_significant, 0.1) * 1.2)
    ) +
    scale_fill_manual(
      name = "Exposure Category",
      values = c("Metabolic" = metabolic_color, "Inflammatory" = inflammatory_color)
    ) +
    labs(
      title = "Proportion of FDR-Significant Causal Associations",
      subtitle = "By Exposure Category and Lung Cancer Subtype",
      x = "Lung Cancer Subtype",
      y = "Proportion of Significant Associations"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 10, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave("results/figures/Figure2_Summary_Barplot.png", 
         plot = p_summary, 
         width = 10, 
         height = 6, 
         dpi = 300)
  cat("  ✓ 已保存: Figure2_Summary_Barplot.png\n")
  
  # 生成热图
  cat("\n生成效应热图\n")
  
  # Beta值热图
  tryCatch({
    # 去重并聚合
    heatmap_data <- mr_summary %>%
      group_by(exposure, outcome) %>%
      summarise(
        beta = mean(beta, na.rm = TRUE),
        fdr_pval = mean(fdr_pval, na.rm = TRUE),
        .groups = "drop"
      )
    
    beta_matrix <- reshape2::acast(heatmap_data, 
                                   exposure ~ outcome, 
                                   value.var = "beta", 
                                   fun.aggregate = mean)
    
    color_palette <- colorRampPalette(colors = c("blue", "white", "red"))(100)
    
    pheatmap(
      beta_matrix, 
      filename = "results/figures/step06_beta_heatmap.png", 
      main = "Beta Coefficients (Effect Sizes)", 
      color = color_palette,
      cluster_rows = TRUE, 
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      fontsize_number = 8,
      width = 8,
      height = 10
    )
    cat("  ✓ 已保存: step06_beta_heatmap.png\n")
    
    # FDR显著性热图
    fdr_matrix <- reshape2::acast(heatmap_data, 
                                  exposure ~ outcome, 
                                  value.var = "fdr_pval", 
                                  fun.aggregate = mean)
    fdr_matrix[is.na(fdr_matrix)] <- 1
    
    # 尝试使用viridis调色板
    if (requireNamespace("viridis", quietly = TRUE)) {
      color_palette_sig <- rev(viridis::magma(100))
    } else {
      color_palette_sig <- colorRampPalette(colors = c("darkred", "yellow", "white"))(100)
    }
    
    pheatmap(
      -log10(fdr_matrix), 
      filename = "results/figures/step06_fdr_significance_heatmap.png", 
      main = "-log10(FDR P-values)", 
      color = color_palette_sig,
      cluster_rows = TRUE, 
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      fontsize_number = 8,
      width = 8,
      height = 10
    )
    cat("  ✓ 已保存: step06_fdr_significance_heatmap.png\n")
    
  }, error = function(e) {
    cat(sprintf("  ⚠ 热图生成失败: %s\n", conditionMessage(e)))
  })
  
  # 生成图注模板
  cat("\n生成图注模板\n")
  
  legend_text <- c(
    "Figure 1: Mendelian Randomization estimates of the causal effects of metabolic and inflammatory traits on lung cancer subtypes.",
    "",
    "This figure presents the odds ratios (ORs) and 95% confidence intervals (CIs) for the association between genetically predicted metabolic traits (orange) and inflammatory markers (blue) with the risk of lung cancer subtypes. The analysis was conducted using a two-sample MR approach with inverse variance weighted (IVW) method. Black dots indicate FDR-adjusted significance (P < 0.05). Error bars represent 95% confidence intervals.",
    "",
    "Abbreviations: BMI, body mass index; CRP, C-reactive protein; HDL, high-density lipoprotein; LDL, low-density lipoprotein; IL-6R, interleukin-6 receptor; WBC, white blood cell count."
  )
  
  writeLines(legend_text, "results/figures/Figure1_Legend_Template.txt")
  cat("  ✓ 已保存: Figure1_Legend_Template.txt\n")
}

cat("\n")

# ===================================================================
# 12. 生成MVMR主文可视化（森林图和热图）
# ===================================================================
cat("【步骤12/13】生成MVMR主文可视化（森林图和热图）...\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 检查MVMR结果
if (exists("mvmr_results_df") && nrow(mvmr_results_df) > 0) {
  
  # 添加category字段（如果不存在）
  if (!"category" %in% colnames(mvmr_results_df)) {
    mvmr_results_df <- mvmr_results_df %>%
      mutate(
        category = case_when(
          exposure %in% names(metabolic_traits) ~ "Metabolic",
          exposure %in% names(inflammatory_traits) ~ "Inflammatory",
          TRUE ~ "Other"
        )
      )
  }
  
  # 12.1 生成MVMR森林图（Forest Plot）
  cat("\n【12.1】生成MVMR森林图...\n")
  
  # 筛选显著结果
  mvmr_forest_data <- mvmr_results_df[mvmr_results_df$significant_nominal == TRUE, ]
  
  if (nrow(mvmr_forest_data) > 0) {
    # 按p值排序
    mvmr_forest_data <- mvmr_forest_data[order(mvmr_forest_data$pval), ]
    
    # 准备森林图数据
    mvmr_forest_data$label <- paste(
      mvmr_forest_data$exposure, "→", mvmr_forest_data$outcome,
      " (", mvmr_forest_data$combination_name, ")"
    )
    mvmr_forest_data$label <- factor(mvmr_forest_data$label, levels = rev(mvmr_forest_data$label))
    
    # 创建森林图
    tryCatch({
      p_mvmr_forest <- ggplot(mvmr_forest_data, 
                               aes(x = label, y = or, ymin = or_lci, ymax = or_uci, 
                                   color = category)) +
        geom_pointrange(size = 0.5, fatten = 3) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c("Metabolic" = "#E69F00", "Inflammatory" = "#56B4E9", "Other" = "#999999")) +
        coord_flip() +
        labs(
          title = "Forest Plot: Significant MVMR Associations",
          x = "Exposure → Outcome (Combination)",
          y = "Odds Ratio (95% CI)",
          color = "Category"
        ) +
        theme_classic(base_size = 10) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 9),
          legend.position = "right"
        )
      
      ggsave("results/figures/main_figures/Figure3_MVMR_ForestPlot.pdf", 
             plot = p_mvmr_forest, width = 12, height = max(8, nrow(mvmr_forest_data) * 0.4), 
             dpi = 600, limitsize = FALSE)
      ggsave("results/figures/main_figures/Figure3_MVMR_ForestPlot.png", 
             plot = p_mvmr_forest, width = 12, height = max(8, nrow(mvmr_forest_data) * 0.4), 
             dpi = 600, limitsize = FALSE)
      
      cat("✓ MVMR森林图已生成:\n")
      cat("  - results/figures/main_figures/Figure3_MVMR_ForestPlot.pdf\n")
      cat("  - results/figures/main_figures/Figure3_MVMR_ForestPlot.png\n")
    }, error = function(e) {
      cat("⚠ MVMR森林图生成失败:", e$message, "\n")
    })
  } else {
    cat("⚠ 没有显著的MVMR结果，跳过森林图生成\n")
  }
  
  # 12.2 生成MVMR对比图（调整前vs调整后，或MVMR vs 单变量MR）
  cat("\n【12.2】生成MVMR对比图（Figure 4格式）...\n")
  
  # 确保ggplot2已加载
  if (!require("ggplot2", quietly = TRUE)) {
    cat("⚠ ggplot2未加载，正在加载...\n")
    library(ggplot2)
  }
  
  # 加载单变量MR结果用于对比
  if (file.exists("results/tables/step05_mr_results_summary.csv")) {
    tryCatch({
      univariate_mr_comp <- read.csv("results/tables/step05_mr_results_summary.csv", 
                                     stringsAsFactors = FALSE)
      
      # 确保列名匹配（转换为小写并检查）
      if (!all(c("exposure", "outcome", "or", "or_lci", "or_uci", "pval") %in% colnames(univariate_mr_comp))) {
        cat("⚠ 单变量MR结果列名不匹配，尝试重命名...\n")
        # 尝试重命名常见的列名变体
        colnames(univariate_mr_comp) <- tolower(colnames(univariate_mr_comp))
        if ("odds_ratio" %in% colnames(univariate_mr_comp)) {
          univariate_mr_comp$or <- univariate_mr_comp$odds_ratio
        }
        if ("or_95_lci" %in% colnames(univariate_mr_comp)) {
          univariate_mr_comp$or_lci <- univariate_mr_comp$or_95_lci
        }
        if ("or_95_uci" %in% colnames(univariate_mr_comp)) {
          univariate_mr_comp$or_uci <- univariate_mr_comp$or_95_uci
        }
      }
      
      # 创建对比数据：MVMR vs 单变量MR（优先选择显著结果，但如果有数据也可以显示所有）
      comparison_data <- mvmr_results_df %>%
        select(exposure, outcome, or, or_lci, or_uci, pval, analysis_type, combination_name) %>%
        filter(!is.na(or), !is.na(or_lci), !is.na(or_uci))
      
      # 确保单变量MR结果有必要的列
      if (all(c("exposure", "outcome", "or") %in% colnames(univariate_mr_comp))) {
        comparison_data <- comparison_data %>%
          left_join(
            univariate_mr_comp %>%
              select(exposure, outcome, or, or_lci, or_uci, pval) %>%
              rename(univ_or = or, univ_or_lci = or_lci, univ_or_uci = or_uci, univ_pval = pval),
            by = c("exposure", "outcome")
          ) %>%
          filter(!is.na(univ_or)) %>%
          mutate(
            label = paste(exposure, "→", outcome),
            category = case_when(
              exposure %in% names(metabolic_traits) ~ "Metabolic",
              exposure %in% names(inflammatory_traits) ~ "Inflammatory",
              TRUE ~ "Other"
            )
          )
      } else {
        cat("⚠ 单变量MR结果缺少必要列，跳过对比图生成\n")
        comparison_data <- data.frame()
      }
      
      if (nrow(comparison_data) > 0) {
        # 创建对比图（单变量MR vs MVMR）
        comparison_long <- bind_rows(
          comparison_data %>%
            select(exposure, outcome, label, category, or, or_lci, or_uci, pval) %>%
            mutate(analysis_type = "MVMR (Adjusted)"),
          comparison_data %>%
            select(exposure, outcome, label, category, or = univ_or, 
                   or_lci = univ_or_lci, or_uci = univ_or_uci, pval = univ_pval) %>%
            mutate(analysis_type = "Univariable MR")
        )
        
        tryCatch({
          p_comparison <- ggplot(comparison_long, 
                                aes(x = label, y = or, ymin = or_lci, ymax = or_uci,
                                    color = analysis_type, shape = analysis_type)) +
            geom_pointrange(position = position_dodge(width = 0.6), size = 0.6, fatten = 2.5) +
            geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
            scale_color_manual(
              values = c("Univariable MR" = "#56B4E9", "MVMR (Adjusted)" = "#E69F00"),
              name = "Analysis Type"
            ) +
            scale_shape_manual(
              values = c("Univariable MR" = 16, "MVMR (Adjusted)" = 17),
              name = "Analysis Type"
            ) +
            coord_flip() +
            labs(
              title = "Figure 4: Independent Causal Effects Adjusted for Smoking",
              subtitle = "Comparison of Univariable MR vs. Multivariable MR (MVMR) Estimates",
              x = "Exposure → Outcome",
              y = "Odds Ratio (95% CI)",
              caption = "MVMR estimates are adjusted for smoking initiation when indicated."
            ) +
            theme_classic(base_size = 10) +
            theme(
              plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
              axis.title = element_text(size = 11, face = "bold"),
              axis.text = element_text(size = 9),
              legend.position = "right",
              legend.title = element_text(size = 10, face = "bold"),
              plot.caption = element_text(size = 8, color = "gray50", hjust = 0)
            )
          
          plot_height <- max(8, nrow(comparison_data) * 0.5)
          ggsave("results/figures/main_figures/Figure4_MVMR_Comparison.pdf", 
                 plot = p_comparison, width = 12, height = plot_height, 
                 dpi = 600, limitsize = FALSE)
          ggsave("results/figures/main_figures/Figure4_MVMR_Comparison.png", 
                 plot = p_comparison, width = 12, height = plot_height, 
                 dpi = 600, limitsize = FALSE)
          
          cat("✓ MVMR对比图已生成:\n")
          cat("  - results/figures/main_figures/Figure4_MVMR_Comparison.pdf\n")
          cat("  - results/figures/main_figures/Figure4_MVMR_Comparison.png\n")
        }, error = function(e) {
          cat("⚠ MVMR对比图生成失败:", e$message, "\n")
        })
      } else {
        cat("⚠ 没有可用于对比的数据\n")
      }
    }, error = function(e) {
      cat("⚠ 无法加载单变量MR结果用于对比:", e$message, "\n")
    })
  } else {
    cat("⚠ 未找到单变量MR结果文件，跳过对比图生成\n")
  }
  
  # 12.3 生成MVMR热图（Heatmap）
  cat("\n【12.3】生成MVMR热图...\n")
  
  # 准备热图数据
  mvmr_heatmap_data <- mvmr_results_df %>%
    select(exposure, outcome, or, pval, category, combination_name) %>%
    mutate(
      log_or = log(or),
      significant = ifelse(pval < 0.05, 1, 0)
    )
  
  # 创建OR值热图矩阵（按结局分组）
  if (require("tidyr", quietly = TRUE) && require("pheatmap", quietly = TRUE)) {
    tryCatch({
      # 对于每个组合，提取主要结果（可以选择代表性结果）
      mvmr_heatmap_summary <- mvmr_heatmap_data %>%
        group_by(exposure, outcome) %>%
        summarise(
          log_or = mean(log_or, na.rm = TRUE),
          significant = ifelse(any(significant == 1), 1, 0),
          .groups = "drop"
        )
      
      # 创建OR值热图矩阵
      or_matrix_mvmr <- mvmr_heatmap_summary %>%
        select(exposure, outcome, log_or) %>%
        pivot_wider(names_from = outcome, values_from = log_or, values_fill = NA)
      
      if (ncol(or_matrix_mvmr) > 1 && nrow(or_matrix_mvmr) > 0) {
        or_mat_mvmr <- as.matrix(or_matrix_mvmr[, -1])
        rownames(or_mat_mvmr) <- or_matrix_mvmr$exposure
        
        # 生成OR值热图
        pdf("results/figures/main_figures/Figure4_MVMR_Heatmap_OR.pdf", width = 8, height = 12)
        pheatmap::pheatmap(
          or_mat_mvmr,
          color = RColorBrewer::brewer.pal(11, "RdBu"),
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          scale = "none",
          main = "MVMR Association Heatmap (log OR)",
          fontsize = 8,
          fontsize_row = 7,
          fontsize_col = 9
        )
        dev.off()
        
        png("results/figures/main_figures/Figure4_MVMR_Heatmap_OR.png", 
            width = 8, height = 12, units = "in", res = 600)
        pheatmap::pheatmap(
          or_mat_mvmr,
          color = RColorBrewer::brewer.pal(11, "RdBu"),
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          scale = "none",
          main = "MVMR Association Heatmap (log OR)",
          fontsize = 8,
          fontsize_row = 7,
          fontsize_col = 9
        )
        dev.off()
        
        # 创建显著性矩阵
        sig_matrix_mvmr <- mvmr_heatmap_summary %>%
          select(exposure, outcome, significant) %>%
          pivot_wider(names_from = outcome, values_from = significant, values_fill = 0)
        
        sig_mat_mvmr <- as.matrix(sig_matrix_mvmr[, -1])
        rownames(sig_mat_mvmr) <- sig_matrix_mvmr$exposure
        
        # 生成显著性热图
        pdf("results/figures/main_figures/Figure4_MVMR_Heatmap_Significance.pdf", width = 8, height = 12)
        pheatmap::pheatmap(
          sig_mat_mvmr,
          color = c("white", "red"),
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          scale = "none",
          main = "MVMR Association Significance (p < 0.05)",
          fontsize = 8,
          fontsize_row = 7,
          fontsize_col = 9
        )
        dev.off()
        
        png("results/figures/main_figures/Figure4_MVMR_Heatmap_Significance.png", 
            width = 8, height = 12, units = "in", res = 600)
        pheatmap::pheatmap(
          sig_mat_mvmr,
          color = c("white", "red"),
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          scale = "none",
          main = "MVMR Association Significance (p < 0.05)",
          fontsize = 8,
          fontsize_row = 7,
          fontsize_col = 9
        )
        dev.off()
        
        cat("✓ MVMR热图已生成:\n")
        cat("  - results/figures/main_figures/Figure4_MVMR_Heatmap_OR.pdf\n")
        cat("  - results/figures/main_figures/Figure4_MVMR_Heatmap_OR.png\n")
        cat("  - results/figures/main_figures/Figure4_MVMR_Heatmap_Significance.pdf\n")
        cat("  - results/figures/main_figures/Figure4_MVMR_Heatmap_Significance.png\n")
      } else {
        cat("⚠ 热图矩阵为空，跳过热图生成\n")
      }
    }, error = function(e) {
      cat("⚠ MVMR热图生成失败:", e$message, "\n")
      cat("  提示：可能需要安装tidyr或pheatmap包\n")
    })
  } else {
    cat("⚠ tidyr或pheatmap包未安装，跳过MVMR热图生成\n")
  }
  
} else {
  cat("⚠ 没有MVMR结果数据，跳过可视化生成\n")
}

cat("\n")

# ===================================================================
# 13. 生成论文级别表格
# ===================================================================
cat("【步骤13/13】生成论文级别表格...\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (exists("mvmr_results_df") && nrow(mvmr_results_df) > 0) {
  
  # 加载单变量MR结果用于对比（计算ΔOR）
  univariate_mr <- NULL
  if (file.exists("results/tables/step05_mr_results_summary.csv")) {
    tryCatch({
      univariate_mr <- read.csv("results/tables/step05_mr_results_summary.csv", 
                               stringsAsFactors = FALSE)
      cat("✓ 已加载单变量MR结果用于对比\n")
    }, error = function(e) {
      cat("⚠ 无法加载单变量MR结果:", e$message, "\n")
    })
  }
  
  # 13.1 主表：MVMR主要结果
  cat("\n【13.1】生成MVMR主结果表...\n")
  
  table_mvmr_main <- mvmr_results_df %>%
    select(outcome, exposure, combination_name, analysis_type,
           beta, se, pval, fdr, or, or_lci, or_uci, n_snps,
           significant_fdr, significant_nominal) %>%
    arrange(pval) %>%
    mutate(
      or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
      pval_formatted = sprintf("%.3e", pval),
      fdr_formatted = ifelse(is.na(fdr), "", sprintf("%.3e", fdr))
    )
  
  write.csv(table_mvmr_main, 
            "results/tables/paper_tables/Table_MVMR_Main_Results.csv", 
            row.names = FALSE)
  cat("✓ 已保存: Table_MVMR_Main_Results.csv\n")
  
  # 13.2 标准MVMR结果表
  cat("\n【13.2】生成标准MVMR结果表...\n")
  
  table_mvmr_standard <- mvmr_results_df %>%
    filter(analysis_type == "standard") %>%
    select(outcome, exposure, combination_name, beta, se, pval, fdr, or, or_lci, or_uci, n_snps) %>%
    arrange(pval)
  
  write.csv(table_mvmr_standard, 
            "results/tables/paper_tables/Table_MVMR_Standard_Results.csv", 
            row.names = FALSE)
  cat("✓ 已保存: Table_MVMR_Standard_Results.csv\n")
  
  # 13.3 Table 3格式：吸烟调整MVMR结果（包含ΔOR对比）
  cat("\n【13.3】生成Table 3格式：Multivariable MR (Adjusted for Smoking)...\n")
  
  if (sum(mvmr_results_df$analysis_type == "smoking_adjusted", na.rm = TRUE) > 0) {
    # 创建Table 3格式表格
    table3_data <- mvmr_results_df %>%
      filter(analysis_type == "smoking_adjusted") %>%
      select(outcome, exposure, combination_name, beta, se, pval, fdr, or, or_lci, or_uci, n_snps)
    
    # 添加单变量MR对比，计算ΔOR
    if (!is.null(univariate_mr) && nrow(univariate_mr) > 0) {
      # 合并单变量MR结果
      table3_data <- table3_data %>%
        left_join(
          univariate_mr %>% 
            select(exposure, outcome, or, beta) %>%
            rename(univ_or = or, univ_beta = beta),
          by = c("exposure", "outcome")
        ) %>%
        mutate(
          delta_or = or - ifelse(is.na(univ_or), 0, univ_or),
          delta_beta = beta - ifelse(is.na(univ_beta), 0, univ_beta),
          delta_or_formatted = ifelse(is.na(delta_or), "", sprintf("%.3f", delta_or))
        )
    } else {
      table3_data$delta_or <- NA
      table3_data$delta_or_formatted <- ""
    }
    
    # 格式化表格，符合Table 3要求
    table3_final <- table3_data %>%
      mutate(
        # 确定协变量调整信息
        covariate_adjusted = case_when(
          grepl("smoking", combination_name, ignore.case = TRUE) ~ "Smoking initiation",
          TRUE ~ "See combination"
        ),
        method = "MVMR",
        or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
        pval_formatted = sprintf("%.3e", pval),
        # 添加显著性标记
        significance = case_when(
          pval < 0.001 ~ "***",
          pval < 0.01 ~ "**",
          pval < 0.05 ~ "*",
          TRUE ~ ""
        ),
        # 格式化p值（带显著性标记）
        pval_with_sig = paste0(sprintf("%.3e", pval), significance)
      ) %>%
      select(
        Exposure = exposure,
        Outcome = outcome,
        Covariate_adjusted = covariate_adjusted,
        Method = method,
        OR_95CI = or_95ci,
        P_value = pval_with_sig,
        DeltaOR_vs_unadjusted = delta_or_formatted,
        N_SNPs = n_snps,
        Combination = combination_name
      ) %>%
      arrange(Outcome, P_value)
    
    write.csv(table3_final, 
              "results/tables/paper_tables/Table_3_MVMR_Adjusted_for_Smoking.csv", 
              row.names = FALSE)
    cat("✓ 已保存: Table_3_MVMR_Adjusted_for_Smoking.csv (主文Table 3格式)\n")
    
    # 同时保存详细版本（包含所有列）
    table_mvmr_smoking <- table3_data %>%
      select(outcome, exposure, combination_name, beta, se, pval, fdr, or, or_lci, or_uci, n_snps, delta_or) %>%
      arrange(pval)
    
    write.csv(table_mvmr_smoking, 
              "results/tables/paper_tables/Table_MVMR_Smoking_Adjusted.csv", 
              row.names = FALSE)
    cat("✓ 已保存: Table_MVMR_Smoking_Adjusted.csv (详细版)\n")
  } else {
    cat("ℹ 没有吸烟调整结果\n")
  }
  
  # 13.4 Supplementary Table S5: Complete MVMR Results
  cat("\n【13.4】生成Supplementary Table S5: Complete MVMR Results...\n")
  
  table_s5 <- mvmr_results_df %>%
    mutate(
      category = case_when(
        exposure %in% names(metabolic_traits) ~ "Metabolic",
        exposure %in% names(inflammatory_traits) ~ "Inflammatory",
        TRUE ~ "Other"
      ),
      method = "IVW-MVMR",
      or_95ci = sprintf("%.3f (%.3f-%.3f)", or, or_lci, or_uci),
      pval_formatted = sprintf("%.3e", pval),
      fdr_formatted = ifelse(is.na(fdr), "", sprintf("%.3e", fdr)),
      significance = case_when(
        !is.na(fdr) & fdr < 0.05 ~ "*** (FDR)",
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**",
        pval < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    select(
      Category = category,
      Exposure = exposure,
      Outcome = outcome,
      Combination = combination_name,
      `Analysis Type` = analysis_type,
      Method = method,
      `OR (95% CI)` = or_95ci,
      `Beta` = beta,
      `SE` = se,
      `P-value` = pval_formatted,
      `FDR` = fdr_formatted,
      `Significance` = significance,
      `N SNPs` = n_snps
    ) %>%
    arrange(Category, Exposure, Outcome)
  
  write.csv(table_s5, 
            "results/tables/paper_tables/Table_S5_MVMR_Complete_Results.csv", 
            row.names = FALSE)
  cat("✓ 已保存: Table_S5_MVMR_Complete_Results.csv (Supplementary Table S5)\n")
  
  # 13.5 生成Excel文件（多工作表）
  cat("\n【13.5】生成Excel文件（多工作表）...\n")
  
  if (require("openxlsx", quietly = TRUE)) {
    tryCatch({
      wb <- createWorkbook()
      
      # 工作表1：主要结果
      if (nrow(table_mvmr_main) > 0) {
        addWorksheet(wb, "MVMR_Main_Results")
        writeData(wb, "MVMR_Main_Results", table_mvmr_main)
        setColWidths(wb, "MVMR_Main_Results", cols = seq_len(ncol(table_mvmr_main)), widths = "auto")
      }
      
      # 工作表2：标准结果
      if (nrow(table_mvmr_standard) > 0) {
        addWorksheet(wb, "MVMR_Standard")
        writeData(wb, "MVMR_Standard", table_mvmr_standard)
        setColWidths(wb, "MVMR_Standard", cols = seq_len(ncol(table_mvmr_standard)), widths = "auto")
      }
      
      # 工作表3：吸烟调整结果
      if (exists("table_mvmr_smoking") && nrow(table_mvmr_smoking) > 0) {
        addWorksheet(wb, "MVMR_Smoking_Adjusted")
        writeData(wb, "MVMR_Smoking_Adjusted", table_mvmr_smoking)
        setColWidths(wb, "MVMR_Smoking_Adjusted", cols = seq_len(ncol(table_mvmr_smoking)), widths = "auto")
      }
      
      # 工作表4：Table 3 (主文格式)
      if (exists("table3_final") && nrow(table3_final) > 0) {
        addWorksheet(wb, "Table_3_Main_Text")
        writeData(wb, "Table_3_Main_Text", table3_final)
        setColWidths(wb, "Table_3_Main_Text", cols = seq_len(ncol(table3_final)), widths = "auto")
      }
      
      # 工作表5：Table S5 (补充材料)
      if (nrow(table_s5) > 0) {
        addWorksheet(wb, "Table_S5_MVMR_Complete")
        writeData(wb, "Table_S5_MVMR_Complete", table_s5)
        setColWidths(wb, "Table_S5_MVMR_Complete", cols = seq_len(ncol(table_s5)), widths = "auto")
      }
      
      saveWorkbook(wb, "results/tables/paper_tables/step06_all_paper_tables.xlsx", overwrite = TRUE)
      cat("✓ Excel文件已保存: step06_all_paper_tables.xlsx (包含所有工作表)\n")
      
    }, error = function(e) {
      cat("⚠ Excel文件生成失败:", e$message, "\n")
    })
  } else {
    cat("⚠ openxlsx包未安装，跳过Excel文件生成\n")
  }
  
  cat("\n✓ 论文级别表格生成完成\n")
  
} else {
  cat("⚠ 没有MVMR结果数据，跳过表格生成\n")
}

cat("\n")

# 14. 生成最终分析报告
cat("【步骤14/14】生成最终分析报告\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

report_lines <- c(
  paste(rep("=", 80), collapse = ""),
  "多变量孟德尔随机化（MVMR）分析报告",
  paste(rep("=", 80), collapse = ""),
  "",
  sprintf("分析日期: %s", Sys.Date()),
  sprintf("分析时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "【分析概况】",
  sprintf("  总MVMR分析数: %d", ifelse(exists("mvmr_results_df"), nrow(mvmr_results_df), 0)),
  sprintf("  标准分析: %d", ifelse(exists("mvmr_results_df"), sum(mvmr_results_df$analysis_type == "standard"), 0)),
  sprintf("  吸烟调整分析: %d", ifelse(exists("mvmr_results_df"), sum(mvmr_results_df$analysis_type == "smoking_adjusted"), 0)),
  sprintf("  失败分析: %d", length(failed_analyses)),
  "",
  "【输出文件】",
  "  MVMR结果表格:",
  "    - results/tables/step06_mvmr_results_final.csv",
  "    - results/tables/step06_mvmr_detailed_final.csv",
  "    - results/tables/step06_smoking_adjusted_final.csv",
  "    - results/tables/step06_failed_analyses.csv",
  "  论文级别表格:",
  "    - results/tables/paper_tables/step06_all_paper_tables.xlsx",
  "    - results/tables/paper_tables/Table_3_MVMR_Adjusted_for_Smoking.csv (主文Table 3)",
  "    - results/tables/paper_tables/Table_S5_MVMR_Complete_Results.csv (补充材料S5)",
  "    - results/tables/paper_tables/Table_MVMR_Main_Results.csv",
  "    - results/tables/paper_tables/Table_MVMR_Standard_Results.csv",
  "    - results/tables/paper_tables/Table_MVMR_Smoking_Adjusted.csv",
  "",
  "  可视化图表:",
  "    单变量MR图表:",
  "    - results/figures/Figure1_Facet_ForestPlot.png (推荐用于论文)",
  "    - results/figures/Figure1_a/b/c_ForestPlot_*.png (按结局独立图)",
  "    - results/figures/Figure2_Summary_Barplot.png",
  "    - results/figures/step06_beta_heatmap.png",
  "    - results/figures/step06_fdr_significance_heatmap.png",
  "    MVMR主文图表:",
  "    - results/figures/main_figures/Figure3_MVMR_ForestPlot.pdf/png",
  "    - results/figures/main_figures/Figure4_MVMR_Comparison.pdf/png (推荐用于主文)",
  "    - results/figures/main_figures/Figure4_MVMR_Heatmap_OR.pdf/png",
  "    - results/figures/main_figures/Figure4_MVMR_Heatmap_Significance.pdf/png",
  "",
  "  文档:",
  "    - results/figures/Figure1_Legend_Template.txt",
  "",
  "【数据质量】",
  sprintf("  工具变量提取成功率: %.1f%% (%d/%d)", 
          100 * length(exposure_instruments) / length(all_exposures_needed),
          length(exposure_instruments),
          length(all_exposures_needed)),
  sprintf("  结局数据提取成功率: %.1f%% (%d/%d)", 
          100 * length(outcome_data_list) / length(outcomes),
          length(outcome_data_list),
          length(outcomes)),
  "",
  paste(rep("=", 80), collapse = "")
)

writeLines(report_lines, "results/tables/step06_analysis_report.txt")
cat("✓ 已保存分析报告: step06_analysis_report.txt\n\n")

# 13. 最终总结
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("🎉 第6步多变量MR分析完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【主要输出】\n")
cat("  ✓ MVMR分析结果（标准 + 吸烟调整）\n")
cat("  ✓ 单变量MR可视化图表（参考）\n")
cat("  ✓ MVMR主文森林图和热图（Figure3和Figure4）\n")
cat("  ✓ 论文级别表格（Excel多工作表 + CSV）\n")
cat("  ✓ 汇总统计图\n")
cat("  ✓ 效应热图（Beta值 + FDR显著性）\n")
cat("  ✓ 图注模板\n")
cat("  ✓ 完整分析报告\n\n")

cat("【关键结果】\n")
if (exists("mvmr_results_df") && nrow(mvmr_results_df) > 0) {
  cat(sprintf("  MVMR分析总数: %d\n", nrow(mvmr_results_df)))
  cat(sprintf("  FDR显著关联: %d (%.1f%%)\n", 
              sum(mvmr_results_df$significant_fdr, na.rm = TRUE),
              100 * sum(mvmr_results_df$significant_fdr, na.rm = TRUE) / nrow(mvmr_results_df)))
  cat(sprintf("  名义显著关联: %d (%.1f%%)\n", 
              sum(mvmr_results_df$significant_nominal, na.rm = TRUE),
              100 * sum(mvmr_results_df$significant_nominal, na.rm = TRUE) / nrow(mvmr_results_df)))
}

cat("\n【文件位置】\n")
cat("  结果表格: results/tables/\n")
cat("  图表文件: results/figures/\n")
cat("  分析报告: results/tables/step06_analysis_report.txt\n\n")

cat("【推荐使用】\n")
cat("  单变量MR主图: Figure1_Facet_ForestPlot.png\n")
cat("  MVMR对比图: Figure4_MVMR_Comparison.png (主文Figure 4)\n")
cat("  MVMR森林图: Figure3_MVMR_ForestPlot.png (补充材料)\n")
cat("  MVMR热图: Figure4_MVMR_Heatmap_OR.png (补充材料)\n")
cat("  补充材料: Figure1_a/b/c_ForestPlot_*.png\n")
cat("  图注文字: Figure1_Legend_Template.txt\n")
cat("  论文表格: \n")
cat("    - Table 3 (主文): Table_3_MVMR_Adjusted_for_Smoking.csv\n")
cat("    - Table S5 (补充): Table_S5_MVMR_Complete_Results.csv\n")
cat("    - Excel汇总: step06_all_paper_tables.xlsx\n\n")

cat("【下一步建议】\n")
cat("  1. 查看 results/figures/ 中的所有图表\n")
cat("  2. 检查 results/tables/step06_mvmr_results_final.csv 了解详细结果\n")
cat("  3. 使用图注模板撰写论文图注\n")
cat("  4. 如有失败分析，查看 step06_failed_analyses.csv\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("分析完成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

