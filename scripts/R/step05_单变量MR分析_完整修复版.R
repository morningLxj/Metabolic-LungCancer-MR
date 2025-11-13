############################################################################
# 第5步：单变量孟德尔随机化分析（完整修复版 - 智能数据加载）
# 修复版本：
#   1. 优先从Step4和Step2数据文件加载工具变量和结局数据
#   2. 对于缺失的工具变量，自动从OpenGWAS在线数据库下载
#   3. 确保所有30个暴露 × 3个结局 = 90个分析全覆盖
############################################################################

cat("第5步：单变量MR分析（智能数据加载：本地优先，缺失时下载）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 1. 加载必要的包
required_packages <- c("TwoSampleMR", "dplyr", "openxlsx", "ggplot2", "gridExtra", 
                       "RColorBrewer", "pheatmap", "tidyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 尝试加载MR-PRESSO（可能不在CRAN上）
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
    cat("  如需安装，请运行：devtools::install_github('rondolab/MR-PRESSO')\n")
  })
}

# 2. 创建输出目录
dirs <- c("data", "results/univariable_mr", "results/tables", 
          "results/figures", "results/figures/supplementary_sensitivity",
          "results/tables/paper_tables", "results/figures/main_figures")
for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# 3. 加载Step4生成的工具变量数据（如果存在）
cat("【步骤3】加载Step4生成的工具变量数据...\n")
instruments_file <- "data/step04_all_instruments.RData"
all_instruments <- list()

if (file.exists(instruments_file)) {
  load(instruments_file)  # 加载 all_instruments 列表
  cat(sprintf("✓ 已加载工具变量数据：%d 个暴露因子\n", length(all_instruments)))
} else {
  cat("⚠ 警告：未找到Step4生成的工具变量文件，将从头开始下载\n")
  all_instruments <- list()
}

# 4. 加载Step2生成的结局数据
cat("【步骤4】加载Step2生成的结局数据...\n")
outcome_data_file <- "results/data/outcome_data_list.RData"
if (!file.exists(outcome_data_file)) {
  stop(sprintf("错误：找不到Step2生成的结局数据文件: %s\n请先运行Step2脚本", outcome_data_file))
}

load(outcome_data_file)  # 加载 outcome_data_list 列表
cat(sprintf("✓ 已加载结局数据：%d 个结局\n", length(outcome_data_list)))

# 5. 定义暴露因子名称映射和GWAS ID映射
# 名称映射：将step05中的名称映射到step04中的名称
exposure_name_mapping <- list(
  # 代谢性状
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
  # 炎症标志物
  CRP = "crp",
  WBC = "wbc",
  IL6 = "il6",
  IL6R = "il6r",
  TNFR1 = "tnfr1"
)

# GWAS ID映射：用于从网上提取工具变量
exposure_gwas_id_mapping <- list(
  # 代谢性状
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
  GGT = "ebi-a-GCST90025966",
  # 炎症标志物
  CRP = "ebi-a-GCST90029070",
  WBC = "ieu-b-30",
  IL6 = "ebi-a-GCST90012005",
  IL6R = "ebi-a-GCST90012025",
  TNFR1 = "ebi-a-GCST90012015"
)

# 定义结局名称映射
outcome_name_mapping <- list(
  lung_cancer_overall = "lung_cancer_overall",
  lung_adenocarcinoma = "lung_adenocarcinoma",
  squamous_cell_lung = "squamous_cell_lung"
)

# 获取所有暴露因子名称（从step05定义）
all_exposure_names <- names(exposure_name_mapping)
outcome_names <- names(outcome_name_mapping)

# 验证哪些暴露因子有工具变量数据
# 方法1：通过映射查找
available_exposures <- character(0)
for (expo_name in all_exposure_names) {
  mapped_name <- exposure_name_mapping[[expo_name]]
  if (mapped_name %in% names(all_instruments)) {
    available_exposures <- c(available_exposures, expo_name)
  }
}

# 方法2：如果方法1失败，尝试直接匹配（大小写不敏感）
if (length(available_exposures) == 0) {
  cat("⚠ 警告：通过映射未找到匹配，尝试直接匹配...\n")
  all_instruments_lower <- tolower(names(all_instruments))
  for (expo_name in all_exposure_names) {
    mapped_name <- tolower(exposure_name_mapping[[expo_name]])
    # 尝试精确匹配
    if (mapped_name %in% all_instruments_lower) {
      available_exposures <- c(available_exposures, expo_name)
    } else {
      # 尝试部分匹配
      matches <- grep(mapped_name, all_instruments_lower, value = FALSE)
      if (length(matches) > 0) {
        available_exposures <- c(available_exposures, expo_name)
        # 更新映射
        actual_name <- names(all_instruments)[matches[1]]
        exposure_name_mapping[[expo_name]] <- actual_name
      }
    }
  }
}

# 方法3：如果还是找不到，直接使用all_instruments中的所有键（作为备选）
if (length(available_exposures) == 0) {
  cat("⚠ 警告：直接匹配也失败，使用all_instruments中的所有键名...\n")
  available_exposures <- names(all_instruments)
  # 为这些键创建身份映射
  for (expo_name in available_exposures) {
    if (!expo_name %in% names(exposure_name_mapping)) {
      exposure_name_mapping[[expo_name]] <- expo_name
    }
  }
}

cat(sprintf("可用暴露因子（有工具变量）: %d/%d\n", 
            length(available_exposures), length(all_exposure_names)))

# 检查缺失的工具变量，并从网上下载
missing_exposures <- setdiff(all_exposure_names, available_exposures)

if (length(missing_exposures) > 0) {
  cat(sprintf("\n【步骤3.1】发现 %d 个暴露因子缺少工具变量，开始从网上下载...\n", 
              length(missing_exposures)))
  
  # 定义从网上下载工具变量的函数
  download_instruments_online <- function(exposure_name, gwas_id, mapped_name) {
    cat(sprintf("  → 下载 %s (GWAS ID: %s)...", exposure_name, gwas_id))
    
    tryCatch({
      # 尝试多个策略提取工具变量
  strategies <- list(
        list(p1 = 5e-8, r2 = 0.001, kb = 10000),
        list(p1 = 5e-7, r2 = 0.001, kb = 10000),
        list(p1 = 5e-6, r2 = 0.01, kb = 5000),
        list(p1 = 5e-5, r2 = 0.05, kb = 5000)
      )
      
      instruments <- NULL
  for (strategy in strategies) {
    tryCatch({
          instruments <- TwoSampleMR::extract_instruments(
            outcomes = gwas_id,
            p1 = strategy$p1,
        clump = TRUE,
        r2 = strategy$r2,
        kb = strategy$kb
      )
      
          if (!is.null(instruments) && nrow(instruments) >= 3) {
            # 验证必要列
            required_cols <- c("SNP", "beta.exposure", "se.exposure")
            if (all(required_cols %in% colnames(instruments))) {
              # 计算F统计量（如果缺失）
              if (!"F_statistic" %in% colnames(instruments)) {
                instruments$F_statistic <- (instruments$beta.exposure^2) / 
                                         (instruments$se.exposure^2)
              }
              cat(sprintf(" ✓ (%d个SNP, p<%.0e)\n", 
                         nrow(instruments), strategy$p1))
        return(instruments)
      }
          }
        }, error = function(e) {
          NULL
        })
      }
      
      cat(" ✗ (提取失败或SNP不足)\n")
  return(NULL)
      
    }, error = function(e) {
      cat(sprintf(" ✗ (错误: %s)\n", conditionMessage(e)))
      return(NULL)
    })
  }
  
  # 为每个缺失的暴露因子下载工具变量
  for (expo_name in missing_exposures) {
    gwas_id <- exposure_gwas_id_mapping[[expo_name]]
    mapped_name <- exposure_name_mapping[[expo_name]]
    
    if (is.null(gwas_id)) {
      cat(sprintf("  ⚠ %s: 未找到GWAS ID，跳过\n", expo_name))
      next
    }
    
    instruments <- download_instruments_online(expo_name, gwas_id, mapped_name)
    
    if (!is.null(instruments) && nrow(instruments) >= 3) {
      # 保存到all_instruments（使用映射后的名称）
      all_instruments[[mapped_name]] <- instruments
      available_exposures <- c(available_exposures, expo_name)
      cat(sprintf("    ✓ 已添加 %s 的工具变量到列表\n", mapped_name))
    }
  }
  
  # 更新可用暴露因子列表
  cat(sprintf("\n✓ 下载完成，当前可用暴露因子: %d/%d\n", 
              length(available_exposures), length(all_exposure_names)))
  
  # 保存更新后的工具变量列表
  if (length(all_instruments) > 0) {
    save(all_instruments, file = instruments_file)
    cat(sprintf("✓ 已保存更新后的工具变量到: %s\n", instruments_file))
  }
} else {
  cat("✓ 所有暴露因子都有工具变量，无需下载\n")
}

# 验证哪些结局有数据
available_outcomes <- character(0)
for (outcome_name in outcome_names) {
  mapped_name <- outcome_name_mapping[[outcome_name]]
  if (mapped_name %in% names(outcome_data_list)) {
    available_outcomes <- c(available_outcomes, outcome_name)
  }
}

# 如果映射不匹配，尝试直接匹配
if (length(available_outcomes) == 0) {
  cat("⚠ 警告：结局名称映射不匹配，尝试直接匹配...\n")
  outcome_data_lower <- tolower(names(outcome_data_list))
  for (outcome_name in outcome_names) {
    mapped_name <- tolower(outcome_name_mapping[[outcome_name]])
    if (mapped_name %in% outcome_data_lower) {
      available_outcomes <- c(available_outcomes, outcome_name)
      # 更新映射
      actual_name <- names(outcome_data_list)[which(outcome_data_lower == mapped_name)[1]]
      outcome_name_mapping[[outcome_name]] <- actual_name
    }
  }
}

# 如果还是找不到，直接使用outcome_data_list中的所有键
if (length(available_outcomes) == 0) {
  cat("⚠ 警告：直接匹配也失败，使用outcome_data_list中的所有键名...\n")
  available_outcomes <- names(outcome_data_list)
  for (outcome_name in available_outcomes) {
    if (!outcome_name %in% names(outcome_name_mapping)) {
      outcome_name_mapping[[outcome_name]] <- outcome_name
    }
  }
}

cat(sprintf("可用结局（有数据）: %d/%d\n\n", 
            length(available_outcomes), length(outcome_names)))

# 定义结局GWAS ID映射（用于从数据库提取数据）
# 注意：这里使用step05中定义的结局名称作为键
outcome_gwas_id_mapping <- list(
  lung_cancer_overall = "ebi-a-GCST004748",    # Lung cancer
  lung_adenocarcinoma = "ieu-a-984",           # Lung adenocarcinoma  
  squamous_cell_lung = "ieu-a-989"             # Squamous cell lung cancer
)

# 定义反向映射：从Step2中的结局名称映射到GWAS ID
outcome_name_to_id <- function(outcome_name) {
  # 尝试直接匹配
  if (outcome_name %in% names(outcome_gwas_id_mapping)) {
    return(outcome_gwas_id_mapping[[outcome_name]])
  }
  
  # 尝试通过Step2中的名称匹配
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

# 6. 辅助函数：计算F统计量和R²
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

# 7. 定义健壮的MR分析函数
perform_mr_analysis_robust <- function(exposure_data, outcome_data, 
                                       exposure_name, outcome_name) {
  tryCatch({
    # 确保数据格式正确并符合TwoSampleMR要求
    # 暴露数据应该已经有exposure相关列，如果没有则添加
    if (!"SNP" %in% names(exposure_data)) {
      if ("rsid" %in% names(exposure_data)) {
        exposure_data$SNP <- exposure_data$rsid
      } else if ("variant" %in% names(exposure_data)) {
        exposure_data$SNP <- exposure_data$variant
      } else {
        stop("无法找到SNP标识符列")
      }
    }
    
    # 确保暴露数据有必要的列
    if (!"id.exposure" %in% names(exposure_data)) {
      exposure_data$id.exposure <- exposure_name
    }
    if (!"exposure" %in% names(exposure_data)) {
      exposure_data$exposure <- exposure_name
    }
    
    # 获取暴露的SNP列表
    if (!"SNP" %in% names(exposure_data)) {
      if ("rsid" %in% names(exposure_data)) {
        exposure_data$SNP <- exposure_data$rsid
      } else if ("variant" %in% names(exposure_data)) {
        exposure_data$SNP <- exposure_data$variant
      } else {
        stop("无法找到SNP标识符列")
      }
    }
    
    exposure_snps <- unique(exposure_data$SNP)
    
    # 检查结局数据质量：如果SNP数量太少（<100），可能是从本地文件加载的不完整数据
    # 这种情况下，我们需要从数据库重新提取特定SNP的数据
    outcome_snp_count <- if ("SNP" %in% names(outcome_data)) {
      length(unique(outcome_data$SNP))
    } else if ("rsid" %in% names(outcome_data)) {
      length(unique(outcome_data$rsid))
    } else {
      0
    }
    
    # 如果结局数据SNP太少，尝试从数据库提取
    if (outcome_snp_count < 50) {
      # 获取结局GWAS ID（使用智能映射函数）
      outcome_id <- outcome_name_to_id(outcome_name)
      if (is.null(outcome_id)) {
        # 尝试从outcome_data中获取ID
        outcome_id <- if ("id.outcome" %in% names(outcome_data) && length(unique(outcome_data$id.outcome)) > 0) {
          unique(outcome_data$id.outcome)[1]
        } else {
          NULL
        }
      }
      
      if (!is.null(outcome_id) && length(exposure_snps) > 0) {
        cat(sprintf("    [提取] 从数据库提取结局数据（暴露SNP: %d个，结局ID: %s）...", 
                   length(exposure_snps), outcome_id))
        tryCatch({
          outcome_data_extracted <- TwoSampleMR::extract_outcome_data(
            snps = exposure_snps,
            outcomes = outcome_id
          )
          
          if (!is.null(outcome_data_extracted) && nrow(outcome_data_extracted) > 0) {
            outcome_data <- outcome_data_extracted
            cat(sprintf(" ✓ (%d个SNP)\n", nrow(outcome_data)))
          } else {
            cat(" ✗ (未提取到数据)\n")
            # 继续使用原始数据尝试
          }
        }, error = function(e) {
          cat(sprintf(" ✗ (错误: %s)\n", conditionMessage(e)))
          # 继续使用原始数据尝试
        })
      }
    }
    
    # 确保结局数据有必要的列
    if (!"SNP" %in% names(outcome_data)) {
      if ("rsid" %in% names(outcome_data)) {
        outcome_data$SNP <- outcome_data$rsid
      } else if ("variant" %in% names(outcome_data)) {
        outcome_data$SNP <- outcome_data$variant
      }
    }
    
    if (!"id.outcome" %in% names(outcome_data)) {
      outcome_data$id.outcome <- outcome_name
    }
    if (!"outcome" %in% names(outcome_data)) {
      outcome_data$outcome <- outcome_name
    }
    
    # 检查数据格式
    if (!"SNP" %in% names(exposure_data) || !"SNP" %in% names(outcome_data)) {
      cat("    ⚠ SNP列缺失\n")
      return(NULL)
    }
    
    # 直接使用harmonise_data，它会自动处理SNP匹配、格式转换和等位基因一致性
    harmonized <- harmonise_data(exposure_data, outcome_data)
    if (is.null(harmonized) || nrow(harmonized) == 0) {
      # 如果harmonise_data失败，可能是数据格式问题，尝试诊断
      exp_snp_count <- ifelse("SNP" %in% names(exposure_data), length(unique(exposure_data$SNP)), 0)
      out_snp_count <- ifelse("SNP" %in% names(outcome_data), length(unique(outcome_data$SNP)), 0)
      cat(sprintf("    ⚠ 数据协调失败（暴露SNP: %d, 结局SNP: %d）\n", 
                 exp_snp_count, out_snp_count))
      return(NULL)
    }
    
    # 根据SNP数量选择分析方法
    n_harmonized <- nrow(harmonized)
    if (n_harmonized >= 3) {
      # 使用所有可用的MR方法
      mr_results <- mr(harmonized, method_list = c("mr_ivw", "mr_egger_regression", 
                                                    "mr_weighted_median", "mr_weighted_mode",
                                                    "mr_simple_mode"))
      heterogeneity <- mr_heterogeneity(harmonized)
      pleiotropy <- mr_pleiotropy_test(harmonized)
      single_snp <- mr_singlesnp(harmonized)
      loo <- mr_leaveoneout(harmonized)
      
      # MR-PRESSO分析（如果包可用且SNP数量足够）
      presso_result <- NULL
      if (require("MRPRESSO", quietly = TRUE)) {
        tryCatch({
          if (n_harmonized >= 10) {  # PRESSO需要至少10个SNP
            # 准备PRESSO所需的数据格式
            presso_data <- data.frame(
              BetaOutcome = harmonized$beta.outcome,
              BetaExposure = harmonized$beta.exposure,
              SdOutcome = harmonized$se.outcome,
              SdExposure = harmonized$se.exposure
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
          }
        }, error = function(e) {
          # 静默失败，不影响主分析
        })
      }
    } else if (n_harmonized == 2) {
      mr_results <- mr(harmonized, method_list = c("mr_ivw", "mr_wald_ratio"))
      heterogeneity <- NULL
      pleiotropy <- NULL
      single_snp <- mr_singlesnp(harmonized)
      loo <- NULL
      presso_result <- NULL
    } else if (n_harmonized == 1) {
      mr_results <- mr(harmonized, method_list = "mr_wald_ratio")
      heterogeneity <- NULL
      pleiotropy <- NULL
      single_snp <- mr_singlesnp(harmonized)
      loo <- NULL
      presso_result <- NULL
    } else {
      return(NULL)
    }
    
    # 计算工具变量强度
    iv_strength <- calculate_instrument_strength(harmonized)
    
    return(list(
      exposure = exposure_name,
      outcome = outcome_name,
      harmonized_data = harmonized,
      mr_results = mr_results,
      heterogeneity = heterogeneity,
      pleiotropy = pleiotropy,
      single_snp = single_snp,
      loo = loo,
      presso = presso_result,
      iv_strength = iv_strength,
      n_snps = nrow(exposure_data),
      n_harmonized = n_harmonized
    ))
  }, error = function(e) {
    cat(sprintf("    ✗ 错误: %s\n", e$message))
    return(NULL)
  })
}

# 7. 定义暴露因子分类
metabolic_exposures <- c("circulating_leptin", "vitamin_D", "HbA1c", "ApoB", "ApoA1", 
                        "IGF1", "ApoB_ApoA1_ratio", "HDL_diameter", "HDL_large", 
                        "remnant_cholesterol", "LDL_small", "BCAA", "HDL_very_large", 
                        "BMI", "HDL_cholesterol", "LDL_cholesterol", "fasting_glucose", 
                        "fasting_insulin", "triglycerides", "GGT", "SBP", "DBP", 
                        "hypertension", "smoking_initiation", "alcohol_drinks")
inflammatory_exposures <- c("CRP", "WBC", "IL6", "IL6R", "TNFR1")

# 8. 执行所有MR分析
cat("【步骤8】开始执行所有MR分析...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

all_mr_results <- list()
extraction_log <- data.frame(
  exposure = character(),
  category = character(),
  n_snps = integer(),
  strategy = character(),
  stringsAsFactors = FALSE
)

analysis_counter <- 0
success_counter <- 0
expected_total <- length(available_exposures) * length(available_outcomes)

for (exposure_name in available_exposures) {
  # 获取映射后的名称
  mapped_exposure_name <- if (exposure_name %in% names(exposure_name_mapping)) {
    exposure_name_mapping[[exposure_name]]
  } else {
    exposure_name
  }
  
  # 检查是否有工具变量数据
  if (!mapped_exposure_name %in% names(all_instruments)) {
    cat(sprintf("⚠ 跳过 %s：未找到工具变量数据\n", exposure_name))
    analysis_counter <- analysis_counter + length(available_outcomes)
    next
  }
  
  exposure_data <- all_instruments[[mapped_exposure_name]]
  if (is.null(exposure_data) || nrow(exposure_data) == 0) {
    cat(sprintf("⚠ 跳过 %s：工具变量数据为空\n", exposure_name))
    analysis_counter <- analysis_counter + length(available_outcomes)
    next
  }
  
  category <- ifelse(exposure_name %in% metabolic_exposures, "Metabolic", "Inflammatory")
  
  cat(sprintf("\n【%s】(%s) - %d SNPs\n", exposure_name, category, nrow(exposure_data)))
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # 记录提取日志
  extraction_strategy <- if("extraction_strategy" %in% names(exposure_data)) {
    exposure_data$extraction_strategy[1]
  } else {
    "Standard"
  }
  
  extraction_log <- rbind(extraction_log, data.frame(
    exposure = exposure_name,
    category = category,
    n_snps = nrow(exposure_data),
    strategy = extraction_strategy,
    stringsAsFactors = FALSE
  ))
  
  # 对每个结局进行MR分析
  for (outcome_name in available_outcomes) {
    # 获取映射后的结局名称
    mapped_outcome_name <- if (outcome_name %in% names(outcome_name_mapping)) {
      outcome_name_mapping[[outcome_name]]
    } else {
      outcome_name
    }
    
    # 检查是否有结局数据
    if (!mapped_outcome_name %in% names(outcome_data_list)) {
      cat(sprintf("  ⚠ 跳过结局 %s：未找到数据\n", outcome_name))
      analysis_counter <- analysis_counter + 1
      next
    }
    
    outcome_data <- outcome_data_list[[mapped_outcome_name]]
    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      cat(sprintf("  ⚠ 跳过结局 %s：数据为空\n", outcome_name))
      analysis_counter <- analysis_counter + 1
      next
    }
    
    analysis_counter <- analysis_counter + 1
    
    cat(sprintf("  [%d/%d] %s -> %s\n", 
               analysis_counter, expected_total, exposure_name, outcome_name))
    
    result <- perform_mr_analysis_robust(
      exposure_data = exposure_data,
      outcome_data = outcome_data,
      exposure_name = exposure_name,
      outcome_name = outcome_name
    )
    
    if (!is.null(result)) {
      key <- paste(exposure_name, outcome_name, sep = "_")
      all_mr_results[[key]] <- result
      success_counter <- success_counter + 1
      
      # 保存单个结果
      save(result, file = paste0("results/univariable_mr/", category, "_", key, ".RData"))
      
      # 生成图表
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
            filename = file.path("results/figures/supplementary_sensitivity", 
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
            filename = file.path("results/figures/supplementary_sensitivity", 
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
              filename = file.path("results/figures/supplementary_sensitivity", 
                                   paste0("loo_", key, ".png")),
              plot = p_loo,
              width = 10, height = plot_height, dpi = 300, limitsize = FALSE
            )
          }
        }
        
        cat("    ✓ 分析和图表生成成功\n")
      }, error = function(e) {
        cat(sprintf("    ⚠ 图表生成失败: %s\n", e$message))
      })
    } else {
      cat("    ✗ 分析失败\n")
    }
  }
}

# 9. 保存完整结果
save(all_mr_results, extraction_log, file = "data/step05_all_results.RData")
cat("\n✓ 完整结果已保存: data/step05_all_results.RData\n")

# 10. 生成汇总表
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
  
  # 添加工具变量强度信息
  for (i in seq_len(nrow(results_summary))) {
    key <- paste(results_summary$exposure[i], results_summary$outcome[i], sep = "_")
    if (key %in% names(all_mr_results)) {
      result <- all_mr_results[[key]]
      if (!is.null(result$iv_strength)) {
        results_summary$mean_f_statistic[i] <- result$iv_strength$mean_f
        results_summary$r_squared[i] <- result$iv_strength$r_squared
      }
    }
  }
}

# 保存汇总表
write.xlsx(results_summary, "results/tables/step05_mr_results_summary.xlsx", row.names = FALSE)
write.csv(results_summary, "results/tables/step05_mr_results_summary.csv", row.names = FALSE)
write.csv(extraction_log, "results/tables/step05_extraction_log.csv", row.names = FALSE)

cat("✓ 汇总表已保存:\n")
cat("  - results/tables/step05_mr_results_summary.xlsx\n")
cat("  - results/tables/step05_mr_results_summary.csv\n")
cat("  - results/tables/step05_extraction_log.csv\n\n")

# ===================================================================
# 10.1 生成增强版表格（论文级别）
# ===================================================================
cat("【步骤10.1】生成增强版论文表格...\n")

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

for (key in names(all_mr_results)) {
  result <- all_mr_results[[key]]
  if (is.null(result) || is.null(result$mr_results)) next
  
  exposure <- result$exposure
  outcome <- result$outcome
  category <- ifelse(exposure %in% metabolic_exposures, "Metabolic", "Inflammatory")
  
  extraction_info <- extraction_log[extraction_log$exposure == exposure, ]
  strategy <- ifelse(nrow(extraction_info) > 0, extraction_info$strategy[1], "Unknown")
  
  # 提取所有MR方法结果
  if (nrow(result$mr_results) > 0) {
    for (i in seq_len(nrow(result$mr_results))) {
      method_row <- result$mr_results[i, ]
      method_name <- method_row$method
      
      beta <- method_row$b
      se <- method_row$se
      pval <- method_row$pval
      nsnp <- ifelse(is.null(method_row$nsnp) || is.na(method_row$nsnp), 
                     result$n_harmonized, method_row$nsnp)
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
        extraction_strategy = strategy,
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
    ivw_beta <- ivw_result$b[1]
    ivw_se <- ivw_result$se[1]
    ivw_pval <- ivw_result$pval[1]
    ivw_nsnp <- ifelse(is.null(ivw_result$nsnp[1]) || is.na(ivw_result$nsnp[1]), 
                       result$n_harmonized, ivw_result$nsnp[1])
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
        het_p <- ivw_het$Q_pval[1]
      }
    }
    
    # 多效性检验
    pleo_p <- NA
    pleo_intercept <- NA
    if (!is.null(result$pleiotropy) && nrow(result$pleiotropy) > 0) {
      pleo_p <- result$pleiotropy$pval[1]
      pleo_intercept <- result$pleiotropy$egger_intercept[1]
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
    
    # Table 2增强版
    table2_enhanced_data[[length(table2_enhanced_data) + 1]] <- data.frame(
      category = category,
      exposure = exposure,
      outcome = outcome,
      method = "IVW",
      n_snps = ivw_nsnp,
      or_95ci = sprintf("%.3f (%.3f-%.3f)", ivw_or, ivw_or_lci, ivw_or_uci),
      pval = ivw_pval,
      egger_or_95ci = ifelse(nrow(egger_result) > 0,
                            sprintf("%.3f (%.3f-%.3f)", 
                                   exp(egger_result$b[1]), 
                                   exp(egger_result$b[1] - 1.96 * egger_result$se[1]),
                                   exp(egger_result$b[1] + 1.96 * egger_result$se[1])),
                            NA),
      egger_pval = ifelse(nrow(egger_result) > 0, egger_result$pval[1], NA),
      wm_or_95ci = ifelse(nrow(wm_result) > 0,
                         sprintf("%.3f (%.3f-%.3f)", 
                                exp(wm_result$b[1]), 
                                exp(wm_result$b[1] - 1.96 * wm_result$se[1]),
                                exp(wm_result$b[1] + 1.96 * wm_result$se[1])),
                         NA),
      wm_pval = ifelse(nrow(wm_result) > 0, wm_result$pval[1], NA),
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
      egger_beta = ifelse(nrow(egger_result) > 0, egger_result$b[1], NA),
      egger_se = ifelse(nrow(egger_result) > 0, egger_result$se[1], NA),
      egger_pval = ifelse(nrow(egger_result) > 0, egger_result$pval[1], NA),
      wm_beta = ifelse(nrow(wm_result) > 0, wm_result$b[1], NA),
      wm_se = ifelse(nrow(wm_result) > 0, wm_result$se[1], NA),
      wm_pval = ifelse(nrow(wm_result) > 0, wm_result$pval[1], NA),
      pleiotropy_intercept = pleo_intercept,
      pleiotropy_p = pleo_p,
      heterogeneity_q = ifelse(!is.null(result$heterogeneity) && 
                               nrow(result$heterogeneity) > 0,
                               result$heterogeneity[result$heterogeneity$method == "Inverse variance weighted", "Q"][1],
                               NA),
      heterogeneity_q_p = het_p,
      presso_global_p = presso_global_p,
      presso_outlier_corrected_or_95ci = ifelse(!is.na(presso_outlier_corrected_or),
                                                sprintf("%.3f (%.3f-%.3f)", 
                                                       presso_outlier_corrected_or,
                                                       presso_outlier_corrected_or_lci,
                                                       presso_outlier_corrected_or_uci),
                                                NA),
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
      extraction_strategy = strategy,
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
  write.csv(table_s1_full, "results/tables/paper_tables/Table_S1_Full_MR_Results.csv", row.names = FALSE)
}
if (nrow(table2_enhanced_full) > 0) {
  write.csv(table2_enhanced_full, "results/tables/paper_tables/Table_2_Enhanced_Main_MR_Results.csv", row.names = FALSE)
}
if (nrow(table5_full) > 0) {
  write.csv(table5_full, "results/tables/paper_tables/Table_5_Sensitivity_Summary.csv", row.names = FALSE)
}
if (nrow(table_s2_full) > 0) {
  write.csv(table_s2_full, "results/tables/paper_tables/Table_S2_Sensitivity_Details.csv", row.names = FALSE)
}
if (nrow(table_s3_full) > 0) {
  write.csv(table_s3_full, "results/tables/paper_tables/Table_S3_Leave_One_Out.csv", row.names = FALSE)
}
if (nrow(table_s6_full) > 0) {
  write.csv(table_s6_full, "results/tables/paper_tables/Table_S6_IV_Strength.csv", row.names = FALSE)
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

saveWorkbook(wb, "results/tables/paper_tables/step05_all_paper_tables.xlsx", overwrite = TRUE)

cat("✓ 论文表格已保存到 results/tables/paper_tables/\n")
cat("  - step05_all_paper_tables.xlsx (多工作表)\n")
cat("  - Table_S1_Full_MR_Results.csv\n")
cat("  - Table_2_Enhanced_Main_MR_Results.csv\n")
cat("  - Table_5_Sensitivity_Summary.csv\n")
cat("  - Table_S2_Sensitivity_Details.csv\n")
cat("  - Table_S3_Leave_One_Out.csv\n")
cat("  - Table_S6_IV_Strength.csv\n\n")

# 11. 生成分析报告
capture.output({
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("单变量MR分析完整报告\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  cat("【分析覆盖度】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat(sprintf("预期总分析数:   %d 个 (%d暴露 × %d结局)\n", 
              expected_total, length(available_exposures), length(available_outcomes)))
  cat(sprintf("实际尝试分析:   %d 个\n", analysis_counter))
  cat(sprintf("成功完成分析:   %d 个\n", success_counter))
  cat(sprintf("成功率:         %.1f%%\n\n", 100 * success_counter / max(1, analysis_counter)))
  
  cat("【工具变量提取】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat(sprintf("总暴露变量:     %d 个\n", length(available_exposures)))
  cat(sprintf("成功提取工具变量: %d 个\n", nrow(extraction_log)))
  cat(sprintf("提取成功率:     %.1f%%\n\n", 100 * nrow(extraction_log) / max(1, length(available_exposures))))
  
  if (nrow(extraction_log) > 0) {
    cat("按提取策略分类:\n")
    print(table(extraction_log$strategy))
    cat("\n")
  }
  
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
      for (i in seq_len(min(5, nrow(sig_fdr)))) {
        cat(sprintf("%d. %s -> %s\n", i, sig_fdr$exposure[i], sig_fdr$outcome[i]))
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
}, file = "results/step05_analysis_report.txt")

cat("✓ 分析报告已保存: results/step05_analysis_report.txt\n\n")

# 12. 生成策略分布图（如果适用）
if (nrow(extraction_log) > 0 && "strategy" %in% names(extraction_log)) {
  cat("生成策略分布图...\n")

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

# ===================================================================
# 13. 生成主文可视化（森林图和热图）
# ===================================================================
cat("\n【步骤13】生成主文可视化（森林图和热图）...\n")

# 13.1 生成森林图（Forest Plot）
if (nrow(results_summary) > 0) {
  # 筛选显著结果（可选：只显示显著结果）
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
          title = "Forest Plot: Significant MR Associations",
          x = "Exposure → Outcome",
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
      
      ggsave("results/figures/main_figures/Figure1_ForestPlot.pdf", 
             plot = p_forest, width = 10, height = max(8, nrow(forest_data) * 0.4), 
             dpi = 600, limitsize = FALSE)
      ggsave("results/figures/main_figures/Figure1_ForestPlot.png", 
             plot = p_forest, width = 10, height = max(8, nrow(forest_data) * 0.4), 
             dpi = 600, limitsize = FALSE)
      
      cat("✓ 森林图已生成:\n")
      cat("  - results/figures/main_figures/Figure1_ForestPlot.pdf\n")
      cat("  - results/figures/main_figures/Figure1_ForestPlot.png\n")
    }, error = function(e) {
      cat("⚠ 森林图生成失败:", e$message, "\n")
    })
  }
}

# 13.2 生成热图（Heatmap）
if (nrow(results_summary) > 0) {
  # 准备热图数据：按结局分组
  heatmap_data <- results_summary %>%
    select(exposure, outcome, or, pval, category) %>%
    mutate(
      log_or = log(or),
      significant = ifelse(pval < 0.05, 1, 0)
    )
  
  # 创建OR值热图矩阵
  if (require("tidyr", quietly = TRUE) && require("pheatmap", quietly = TRUE)) {
    tryCatch({
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
      pdf("results/figures/main_figures/Figure2_Heatmap_OR.pdf", width = 8, height = 12)
      pheatmap::pheatmap(
        or_mat,
        color = RColorBrewer::brewer.pal(11, "RdBu"),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        scale = "none",
        main = "MR Association Heatmap (log OR)",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 9
      )
      dev.off()
      
      png("results/figures/main_figures/Figure2_Heatmap_OR.png", 
          width = 8, height = 12, units = "in", res = 600)
      pheatmap::pheatmap(
        or_mat,
        color = RColorBrewer::brewer.pal(11, "RdBu"),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        scale = "none",
        main = "MR Association Heatmap (log OR)",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 9
      )
      dev.off()
      
      # 生成显著性热图
      pdf("results/figures/main_figures/Figure2_Heatmap_Significance.pdf", width = 8, height = 12)
      pheatmap::pheatmap(
        sig_mat,
        color = c("white", "red"),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        scale = "none",
        main = "MR Association Significance (p < 0.05)",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 9
      )
      dev.off()
      
      png("results/figures/main_figures/Figure2_Heatmap_Significance.png", 
          width = 8, height = 12, units = "in", res = 600)
      pheatmap::pheatmap(
        sig_mat,
        color = c("white", "red"),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        scale = "none",
        main = "MR Association Significance (p < 0.05)",
        fontsize = 8,
        fontsize_row = 7,
        fontsize_col = 9
      )
      dev.off()
      
      cat("✓ 热图已生成:\n")
      cat("  - results/figures/main_figures/Figure2_Heatmap_OR.pdf\n")
      cat("  - results/figures/main_figures/Figure2_Heatmap_OR.png\n")
      cat("  - results/figures/main_figures/Figure2_Heatmap_Significance.pdf\n")
      cat("  - results/figures/main_figures/Figure2_Heatmap_Significance.png\n")
    }, error = function(e) {
      cat("⚠ 热图生成失败:", e$message, "\n")
      cat("  提示：可能需要安装tidyr或pheatmap包\n")
    })
} else {
    cat("⚠ tidyr或pheatmap包未安装，跳过热图生成\n")
  }
}

# 14. 最终总结
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("第5步分析完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【最终统计】\n")
cat(sprintf("预期分析:       %d 个\n", expected_total))
cat(sprintf("实际尝试:       %d 个\n", analysis_counter))
cat(sprintf("成功完成:       %d 个\n", success_counter))
cat(sprintf("成功率:         %.1f%%\n", 100 * success_counter / max(1, analysis_counter)))
if (nrow(results_summary) > 0) {
cat(sprintf("FDR显著关联:    %d 个\n", sum(results_summary$significant_fdr, na.rm = TRUE)))
cat(sprintf("名义显著关联:   %d 个\n\n", sum(results_summary$significant_nominal, na.rm = TRUE)))
}

cat("【保存的文件】\n")
cat("主要结果:\n")
cat("  - data/step05_all_results.RData\n")
cat("汇总表格:\n")
cat("  - results/tables/step05_mr_results_summary.xlsx\n")
cat("  - results/tables/step05_mr_results_summary.csv\n")
cat("  - results/tables/step05_extraction_log.csv\n")
cat("详细报告:\n")
cat("  - results/step05_analysis_report.txt\n")
cat("图表文件:\n")
cat("  - results/figures/strategy_comparison_plot.pdf\n")
cat("  - results/figures/strategy_comparison_plot.png\n")
cat("主文图表:\n")
cat("  - results/figures/main_figures/Figure1_ForestPlot.pdf/png\n")
cat("  - results/figures/main_figures/Figure2_Heatmap_OR.pdf/png\n")
cat("  - results/figures/main_figures/Figure2_Heatmap_Significance.pdf/png\n")
cat("敏感性分析图表:\n")
cat("  - results/figures/supplementary_sensitivity/*.png\n")
cat("论文表格:\n")
cat("  - results/tables/paper_tables/step05_all_paper_tables.xlsx\n")
cat("  - results/tables/paper_tables/*.csv\n\n")

cat("分析完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n提示: 可以使用以下命令查看结果:\n")
cat("  - load('data/step05_all_results.RData')\n")
cat("  - View(results_summary)\n")
cat("  - readLines('results/step05_analysis_report.txt') %>% cat(sep='\\n')\n")
