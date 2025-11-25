############################################################################
# Step05 诊断脚本：找出失败的MR分析
# 用途：识别哪些暴露-结局对的分析失败了
############################################################################

cat("Step05 诊断：识别失败的MR分析\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 1. 加载必要的包
required_packages <- c("TwoSampleMR", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 2. 加载已有的结果
cat("正在检查结果文件...\n")
results_file <- "data/step05_all_results.RData"
cat(sprintf("结果文件路径: %s\n", results_file))

if (!file.exists(results_file)) {
  stop("错误：找不到结果文件 data/step05_all_results.RData\n请先运行 step05_单变量MR分析_完整修复版.R")
}

cat("正在加载结果文件...\n")
load(results_file)  # 加载 all_mr_results 和 extraction_log

# 检查加载的对象是否存在
if (!exists("all_mr_results")) {
  stop("错误：结果文件中未找到 all_mr_results 对象。请确保已正确运行 step05_单变量MR分析_完整修复版.R")
}
cat(sprintf("✓ 成功加载结果文件，包含 %d 个分析结果\n", length(all_mr_results)))

# 3. 从结果文件中提取已完成的分析信息
# 注意：我们不需要运行整个Step5脚本，只需要从数据文件中加载信息

# 从all_mr_results中提取已成功的分析
if (length(all_mr_results) == 0) {
  cat("警告：all_mr_results 为空，可能没有完成任何分析\n")
  completed_analyses <- character(0)
} else {
  completed_analyses <- names(all_mr_results)
}

# 4. 获取所有预期的暴露-结局组合
# 需要加载Step5中的可用暴露和结局列表
instruments_file <- "data/step04_all_instruments.RData"
outcome_file <- "results/data/outcome_data_list.RData"

if (!file.exists(instruments_file) || !file.exists(outcome_file)) {
  stop("错误：找不到必要的数据文件\n")
}

load(instruments_file)
load(outcome_file)

# 定义暴露和结局名称（与Step5保持一致）
all_exposure_names <- c(
  "circulating_leptin", "vitamin_D", "HbA1c", "ApoB", "ApoA1", 
  "IGF1", "ApoB_ApoA1_ratio", "HDL_diameter", "HDL_large", 
  "remnant_cholesterol", "LDL_small", "BCAA", "HDL_very_large", 
  "BMI", "HDL_cholesterol", "LDL_cholesterol", "fasting_glucose", 
  "fasting_insulin", "triglycerides", "GGT", "SBP", "DBP", 
  "hypertension", "smoking_initiation", "alcohol_drinks",
  "CRP", "WBC", "IL6", "IL6R", "TNFR1"
)

outcome_names <- c(
  "lung_cancer_overall",
  "lung_cancer_adenocarcinoma", 
  "lung_cancer_squamous"
)

# 定义GWAS ID映射（与Step5保持一致）
outcome_gwas_id_mapping <- list(
  "Lung cancer overall" = "ebi-a-GCST004748",
  "Lung cancer (Adenocarcinoma)" = "ieu-a-984",
  "Lung cancer (Squamous cell lung cancer)" = "ieu-a-989",
  "lung_cancer_overall" = "ebi-a-GCST004748",
  "lung_cancer_adenocarcinoma" = "ieu-a-984",
  "lung_cancer_squamous" = "ieu-a-989"
)

# 定义反向映射函数
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

# 名称映射
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

# 获取实际可用的暴露和结局
available_exposures <- character(0)
for (expo_name in all_exposure_names) {
  mapped_name <- exposure_name_mapping[[expo_name]]
  if (mapped_name %in% names(all_instruments)) {
    available_exposures <- c(available_exposures, expo_name)
  }
}

available_outcomes <- character(0)
for (outcome_name in outcome_names) {
  mapped_name <- outcome_name_mapping[[outcome_name]]
  if (mapped_name %in% names(outcome_data_list)) {
    available_outcomes <- c(available_outcomes, outcome_name)
  }
}

# 如果还是空的，尝试直接匹配
if (length(available_exposures) == 0) {
  available_exposures <- names(all_instruments)
}
if (length(available_outcomes) == 0) {
  available_outcomes <- names(outcome_data_list)
}

# 5. 生成所有预期的分析组合
expected_analyses <- character(0)
for (expo in available_exposures) {
  for (outcome in available_outcomes) {
    key <- paste(expo, outcome, sep = "_")
    expected_analyses <- c(expected_analyses, key)
  }
}

# 6. 找出缺失的分析
missing_analyses <- setdiff(expected_analyses, completed_analyses)

cat("【诊断结果】\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("预期分析总数:   %d 个\n", length(expected_analyses)))
cat(sprintf("已完成分析:      %d 个\n", length(completed_analyses)))
cat(sprintf("缺失分析:        %d 个\n\n", length(missing_analyses)))

if (length(missing_analyses) > 0) {
  cat("【缺失的分析详情】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  missing_df <- data.frame(
    exposure = character(),
    outcome = character(),
    stringsAsFactors = FALSE
  )
  
  # 更智能地解析键名：使用已知的暴露和结局列表
  for (key in missing_analyses) {
    exposure_name <- NULL
    outcome_name <- NULL
    
    # 首先尝试匹配完整的结局名称
    for (outcome in available_outcomes) {
      if (grepl(paste0("_", outcome, "$"), key) || grepl(paste0("^", outcome, "_"), key)) {
        outcome_name <- outcome
        # 移除结局名称部分得到暴露名称
        exposure_pattern <- paste0("_", outcome, "$")
        exposure_name <- gsub(exposure_pattern, "", key)
        break
      }
    }
    
    # 如果还是没找到，尝试通过匹配结局的关键词
    if (is.null(outcome_name)) {
      if (grepl("_overall$", key)) {
        outcome_name <- "lung_cancer_overall"
        exposure_name <- gsub("_overall$", "", key)
      } else if (grepl("_adenocarcinoma$", key)) {
        outcome_name <- "lung_cancer_adenocarcinoma"
        exposure_name <- gsub("_adenocarcinoma$", "", key)
      } else if (grepl("_squamous$", key)) {
        outcome_name <- "lung_cancer_squamous"
        exposure_name <- gsub("_squamous$", "", key)
      } else {
        # 最后尝试：匹配已知的暴露名称
        for (expo in available_exposures) {
          if (grepl(paste0("^", expo, "_"), key)) {
            exposure_name <- expo
            outcome_name <- gsub(paste0("^", expo, "_"), "", key)
            break
          }
        }
      }
    }
    
    # 如果还是找不到，使用简单的分割方式
    if (is.null(exposure_name) || is.null(outcome_name)) {
      parts <- strsplit(key, "_")[[1]]
      outcome_patterns <- c("overall", "adenocarcinoma", "squamous")
      outcome_start <- 0
      for (i in seq_along(parts)) {
        if (parts[i] %in% outcome_patterns || 
            (i < length(parts) && paste(parts[i], parts[i+1], sep="_") %in% c("lung_cancer"))) {
          outcome_start <- i
          break
        }
      }
      
      if (outcome_start > 1) {
        exposure_name <- paste(parts[1:(outcome_start-1)], collapse = "_")
        outcome_name <- paste(parts[outcome_start:length(parts)], collapse = "_")
      } else if (length(parts) >= 2) {
        # 假设最后2-3部分是结局
        if (length(parts) >= 4) {
          exposure_name <- paste(parts[1:(length(parts)-3)], collapse = "_")
          outcome_name <- paste(parts[(length(parts)-2):length(parts)], collapse = "_")
        } else {
          exposure_name <- parts[1]
          outcome_name <- paste(parts[2:length(parts)], collapse = "_")
        }
      } else {
        exposure_name <- key
        outcome_name <- "unknown"
      }
    }
    
    missing_df <- rbind(missing_df, data.frame(
      exposure = exposure_name,
      outcome = outcome_name,
      stringsAsFactors = FALSE
    ))
  }
  
  print(missing_df)
  
  # 保存缺失分析列表
  write.csv(missing_df, "results/tables/step05_missing_analyses.csv", row.names = FALSE)
  cat(sprintf("\n✓ 缺失分析列表已保存: results/tables/step05_missing_analyses.csv\n\n"))
  
  # 生成修复脚本的命令
  cat("【修复建议】\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("运行以下命令修复缺失的分析：\n")
  cat("\n  source('step05_fix_missing.R')\n\n")
  
  # 保存缺失分析的键名供修复脚本使用
  save(missing_analyses, missing_df, file = "data/step05_missing_analyses.RData")
  cat("✓ 缺失分析数据已保存: data/step05_missing_analyses.RData\n")
  
} else {
  cat("✓ 所有分析均已完成，无需修复！\n\n")
}

cat("\n诊断完成！\n")

