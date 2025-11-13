############################################################################
# 步骤7：中介分析（完整整合版）
# 孟德尔随机化研究 - 代谢性状、炎症标志物与肺癌亚型
# 
# 整合内容：
#   1. 参考step05的数据集提取方式，确保无偏差
#   2. 整合混合策略和穷举法
#   3. 统一数据加载和工具变量提取逻辑
############################################################################

cat("步骤7：中介分析（完整整合版）\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("注意：这是高级分析，将检验炎症标志物的中介作用\n\n")

# 加载必要的包
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(tidyr)
})

# 声明全局变量（消除dplyr NSE的linter警告）
utils::globalVariables(c("exposure", "outcome", "method"))

# 创建输出目录
dirs <- c("results/mediation_analysis", "results/figures", "results/tables", 
          "results/figures/step07_publication",  # 用于发表级图表
          "data", "data/mediation_cache")
for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# 设置超时和选项
options(timeout = 300)

# ============================================================================
# 步骤1：加载前序步骤的数据（参考step05的方式）
# ============================================================================

cat("【步骤1】加载数据\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 1.1 加载Step4的工具变量（优先使用本地文件，与step05保持一致）
instruments_file <- "data/step04_all_instruments.RData"
all_instruments <- list()

if (file.exists(instruments_file)) {
  load(instruments_file)
  # 检查加载后的变量是否存在
  if (!exists("all_instruments") || is.null(all_instruments)) {
    all_instruments <- list()
    cat("⚠ Step4工具变量文件存在但数据为空，将在分析中按需提取\n")
  } else {
    cat(sprintf("✓ 已加载Step4工具变量：%d 个暴露因子\n", length(all_instruments)))
  }
} else {
  cat("⚠ 未找到Step4工具变量文件，将在分析中按需提取\n")
  all_instruments <- list()
}

# 1.2 加载Step2的结局数据
outcome_data_file <- "results/data/outcome_data_list.RData"
if (!file.exists(outcome_data_file)) {
  stop(sprintf("错误：找不到Step2结局数据文件：%s\n请先运行Step2", outcome_data_file))
}
load(outcome_data_file)
cat(sprintf("✓ 已加载Step2结局数据：%d 个结局\n", length(outcome_data_list)))

# 1.3 加载Step5的单变量MR结果
if (!file.exists("results/tables/step05_mr_results_summary.csv")) {
  stop("错误：找不到 step05_mr_results_summary.csv，请先运行步骤5")
}
univariable_results <- read.csv("results/tables/step05_mr_results_summary.csv",
                                stringsAsFactors = FALSE)
cat(sprintf("✓ 已加载Step5单变量MR结果：%d 个分析结果\n\n", nrow(univariable_results)))

# ============================================================================
# 步骤2：定义数据集映射（与step05保持一致，确保无偏差）
# ============================================================================

cat("【步骤2】定义数据集映射\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 暴露名称映射（与step05完全一致）
exposure_name_mapping <- list(
  circulating_leptin = "circulating_leptin",
  vitamin_D = "vitamin_d",
  HbA1c = "hba1c",
  ApoB = "apob",
  ApoA1 = "apoa1",
  IGF1 = "igf1",
  ApoB_ApoA1_ratio = "apob_apoa1_ratio",
  HDL_diameter = "hdl_diameter",
  HDL_large = "hdl_large",
  remnant_cholesterol = "remnant_cholesterol",
  LDL_small = "ldl_small",
  BCAA = "bcaa",
  HDL_very_large = "hdl_very_large",
  BMI = "bmi",
  HDL_cholesterol = "hdl_cholesterol",
  LDL_cholesterol = "ldl_cholesterol",
  smoking_initiation = "smoking_initiation",
  alcohol_drinks = "alcohol_drinks",
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

# GWAS ID映射（与step05完全一致）
exposure_gwas_id_mapping <- list(
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
  CRP = "ebi-a-GCST90029070",
  WBC = "ieu-b-30",
  IL6 = "ebi-a-GCST90012005",
  IL6R = "ebi-a-GCST90012025",
  TNFR1 = "ebi-a-GCST90012015"
)

# 代谢性状（25个）
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

# 炎症标志物（5个）
inflammatory_traits <- list(
  CRP = "ebi-a-GCST90029070",
  WBC = "ieu-b-30",
  IL6 = "ebi-a-GCST90012005",
  IL6R = "ebi-a-GCST90012025",
  TNFR1 = "ebi-a-GCST90012015"
)

# 结局变量（3个）
outcomes <- list(
  lung_cancer_overall = "ebi-a-GCST004748",
  lung_adenocarcinoma = "ieu-a-984",
  squamous_cell_lung = "ieu-a-989"
)

cat(sprintf("✓ 代谢性状: %d 个\n", length(metabolic_traits)))
cat(sprintf("✓ 炎症标志物: %d 个\n", length(inflammatory_traits)))
cat(sprintf("✓ 结局变量: %d 个\n\n", length(outcomes)))

# ============================================================================
# 步骤3：工具变量提取函数（参考step05的方式，确保无偏差）
# ============================================================================

cat("【步骤3】定义工具变量提取函数\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 工具变量提取函数（与step05完全一致的策略）
get_instruments_robust <- function(exposure_name, gwas_id, mapped_name) {
  # 步骤1：优先从本地all_instruments加载
  if (mapped_name %in% names(all_instruments)) {
    instruments <- all_instruments[[mapped_name]]
    if (!is.null(instruments) && is.data.frame(instruments) && nrow(instruments) >= 3) {
      cat(sprintf("    [本地] %s: %d个SNP\n", exposure_name, nrow(instruments)))
      return(instruments)
    }
  }
  
  # 步骤2：如果本地没有，尝试从网上下载（使用与step05相同的策略）
  cat(sprintf("    [提取] %s (ID: %s)...", exposure_name, gwas_id))
  
  tryCatch({
    # 多个策略尝试（与step05完全一致）
    strategies <- list(
      list(p1 = 5e-8, r2 = 0.001, kb = 10000),
      list(p1 = 5e-7, r2 = 0.001, kb = 10000),
      list(p1 = 5e-6, r2 = 0.01, kb = 5000),
      list(p1 = 5e-5, r2 = 0.05, kb = 5000)
    )
    
    for (strategy in strategies) {
      tryCatch({
        instruments <- TwoSampleMR::extract_instruments(
          exposures = gwas_id,
          p1 = strategy$p1,
          clump = TRUE,
          r2 = strategy$r2,
          kb = strategy$kb
        )
        
        if (!is.null(instruments) && is.data.frame(instruments) && nrow(instruments) >= 3) {
          # 验证必要列（与step05一致）
          required_cols <- c("SNP", "beta.exposure", "se.exposure")
          if (all(required_cols %in% colnames(instruments))) {
            # 计算F统计量（如果缺失）
            if (!"F_statistic" %in% colnames(instruments)) {
              instruments$F_statistic <- (instruments$beta.exposure^2) / 
                                       (instruments$se.exposure^2)
            }
            # 保存到all_instruments以便后续使用
            all_instruments[[mapped_name]] <<- instruments
            cat(sprintf(" ✓ (%d个SNP, p<%.0e)\n", 
                       nrow(instruments), strategy$p1))
            return(instruments)
          }
        }
      }, error = function(e) {
        NULL
      })
      Sys.sleep(0.5)  # 避免API限流
    }
    
    cat(sprintf(" ✗ (提取失败: 无法从GWAS ID %s获取足够SNP)\n", gwas_id))
    return(NULL)
    
  }, error = function(e) {
    cat(sprintf(" ✗ (错误: %s)\n", conditionMessage(e)))
    return(NULL)
  })
}

# 带缓存的工具变量获取函数
instrument_cache <- list()
get_instruments_cached <- function(exposure_name, gwas_id, mapped_name, category) {
  cache_key <- paste(mapped_name, category, sep = "_")
  
  # 检查内存缓存
  if (cache_key %in% names(instrument_cache)) {
    return(instrument_cache[[cache_key]])
  }
  
  # 调用提取函数
  instruments <- get_instruments_robust(exposure_name, gwas_id, mapped_name)
  
  if (!is.null(instruments) && nrow(instruments) >= 3) {
    instrument_cache[[cache_key]] <<- instruments
  }
  
  return(instruments)
}

cat("✓ 工具变量提取函数已定义\n\n")

# ============================================================================
# 步骤4：定义中介分析路径
# ============================================================================

cat("【步骤4】定义中介分析路径\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 4.1 扩展智能选择策略（26条核心策略，基于文献证据和生物学合理性）
smart_selection_paths <- list(
  # 1-4. 肥胖相关炎症（BMI → 4种炎症标志物）
  list(exp = "BMI", med = c("CRP", "IL6", "WBC", "TNFR1"), 
       rationale = "Obesity-induced systemic inflammation and cytokine elevation"),
  
  # 5-7. 维生素D抗炎效应
  list(exp = "vitamin_D", med = c("CRP", "IL6", "TNFR1"), 
       rationale = "Vitamin D deficiency associated with increased inflammatory markers"),
  
  # 8-10. 胰岛素抵抗和炎症（3条路径）
  list(exp = "fasting_insulin", med = c("CRP", "IL6", "IL6R", "WBC"), 
       rationale = "Insulin resistance triggers inflammatory cascade"),
  list(exp = "HbA1c", med = c("CRP", "IL6"), 
       rationale = "Chronic hyperglycemia promotes inflammation"),
  list(exp = "fasting_glucose", med = c("CRP", "IL6"), 
       rationale = "Glucose dysregulation and inflammation"),
  
  # 11-17. 脂蛋白代谢（7条路径）
  list(exp = "HDL_cholesterol", med = c("CRP", "IL6R"), 
       rationale = "HDL anti-inflammatory properties and reverse cholesterol transport"),
  list(exp = "ApoA1", med = c("CRP", "TNFR1"), 
       rationale = "ApoA1 deficiency associated with inflammation"),
  list(exp = "triglycerides", med = c("CRP", "IL6", "WBC"), 
       rationale = "Hypertriglyceridemia linked to inflammatory response"),
  list(exp = "LDL_cholesterol", med = c("CRP", "IL6"), 
       rationale = "Oxidized LDL triggers inflammatory pathways"),
  list(exp = "ApoB", med = c("CRP", "IL6"), 
       rationale = "ApoB-containing lipoproteins promote inflammation"),
  list(exp = "HDL_large", med = c("CRP", "IL6"), 
       rationale = "Large HDL particles have superior anti-inflammatory function"),
  list(exp = "HDL_diameter", med = c("CRP"), 
       rationale = "HDL size inversely correlates with inflammatory risk"),
  
  # 18-19. 更多脂蛋白亚型
  list(exp = "HDL_very_large", med = c("CRP"), 
       rationale = "Very large HDL and anti-inflammatory capacity"),
  list(exp = "remnant_cholesterol", med = c("CRP", "IL6"), 
       rationale = "Remnant cholesterol promotes arterial inflammation"),
  list(exp = "LDL_small", med = c("CRP", "IL6"), 
       rationale = "Small dense LDL particles highly atherogenic and pro-inflammatory"),
  list(exp = "ApoB_ApoA1_ratio", med = c("CRP", "IL6"), 
       rationale = "ApoB/ApoA1 ratio reflects lipid-driven inflammation"),
  
  # 20. 生长因子相关炎症
  list(exp = "IGF1", med = c("IL6", "TNFR1", "CRP"), 
       rationale = "IGF-1 signaling pathways interact with inflammatory cytokines"),
  
  # 21-22. 血压和炎症
  list(exp = "SBP", med = c("CRP", "IL6", "TNFR1"), 
       rationale = "Hypertension associated with vascular inflammation"),
  list(exp = "DBP", med = c("CRP", "IL6"), 
       rationale = "Diastolic pressure and endothelial inflammation"),
  
  # 23. 肝功能和炎症
  list(exp = "GGT", med = c("CRP", "IL6", "WBC"), 
       rationale = "Liver enzymes reflect metabolic stress and inflammation"),
  
  # 24. 吸烟相关炎症（重点）
  list(exp = "smoking_initiation", med = c("CRP", "WBC", "IL6R", "IL6"), 
       rationale = "Smoking-induced chronic inflammation and immune activation"),
  
  # 25. 酒精相关炎症
  list(exp = "alcohol_drinks", med = c("CRP", "GGT", "WBC"), 
       rationale = "Excessive alcohol consumption triggers inflammatory response"),
  
  # 26. 其他：代谢综合征、氨基酸代谢、脂肪因子
  list(exp = "hypertension", med = c("CRP", "IL6", "WBC"), 
       rationale = "Hypertension as component of metabolic inflammation"),
  list(exp = "BCAA", med = c("CRP", "IL6"), 
       rationale = "Branched-chain amino acids linked to inflammatory pathways"),
  list(exp = "circulating_leptin", med = c("CRP", "IL6", "TNFR1"), 
       rationale = "Leptin as pro-inflammatory adipokine")
)

# 4.2 穷举+预筛策略（基于单变量MR结果）
exhaustive_screening_paths <- list()

# 结局名称映射（step05的CSV使用不同的格式）
outcome_name_to_csv <- list(
  lung_cancer_overall = "Lung cancer",
  lung_adenocarcinoma = "Lung adenocarcinoma",
  squamous_cell_lung = "Squamous cell lung cancer"
)

# 预筛函数
pre_screen_pathway <- function(exp_name, med_name, outcome_name) {
  # 转换结局名称为CSV格式
  outcome_csv_name <- if (outcome_name %in% names(outcome_name_to_csv)) {
    outcome_name_to_csv[[outcome_name]]
  } else {
    outcome_name  # 如果不在映射中，直接使用
  }
  
  # 检查暴露->结局
  exp_out_check <- univariable_results %>%
    filter(exposure == exp_name, outcome == outcome_csv_name)
  exp_out_sig <- if(nrow(exp_out_check) > 0) exp_out_check$pval[1] < 0.1 else FALSE
  
  # 检查暴露->中介（中介作为"结局"，需要在univariable_results中查找）
  # 注意：中介名称需要匹配CSV中的格式
  # 在step05中，炎症标志物可能作为暴露，所以需要检查exposure列
  exp_med_check <- univariable_results %>%
    filter(exposure == exp_name, outcome == med_name)
  # 如果没找到，尝试将中介名称转换为可能的格式
  if (nrow(exp_med_check) == 0) {
    # 尝试常见的中介名称格式转换
    med_name_variants <- c(
      med_name,
      toupper(med_name),
      tools::toTitleCase(med_name),
      gsub("_", " ", tools::toTitleCase(med_name))
    )
    for (med_var in med_name_variants) {
      exp_med_check <- univariable_results %>%
        filter(exposure == exp_name, outcome == med_var)
      if (nrow(exp_med_check) > 0) break
    }
  }
  exp_med_sig <- if(nrow(exp_med_check) > 0) exp_med_check$pval[1] < 0.1 else FALSE
  
  # 检查中介->结局（中介作为"暴露"）
  med_out_check <- univariable_results %>%
    filter(exposure == med_name, outcome == outcome_csv_name)
  # 如果没找到，尝试中介名称的变体
  if (nrow(med_out_check) == 0) {
    med_name_variants <- c(
      med_name,
      toupper(med_name),
      tools::toTitleCase(med_name),
      gsub("_", " ", tools::toTitleCase(med_name))
    )
    for (med_var in med_name_variants) {
      med_out_check <- univariable_results %>%
        filter(exposure == med_var, outcome == outcome_csv_name)
      if (nrow(med_out_check) > 0) break
    }
  }
  med_out_sig <- if(nrow(med_out_check) > 0) med_out_check$pval[1] < 0.1 else FALSE
  
  # 至少两个路径显著（放宽条件：至少一个路径显著即可）
  sig_count <- sum(c(exp_out_sig, exp_med_sig, med_out_sig), na.rm = TRUE)
  return(sig_count >= 1)  # 改为至少1个路径显著，更宽松的筛选
}

# 生成穷举路径
if (nrow(univariable_results) > 0) {
  # 首先打印一些调试信息
  cat("  调试信息：检查单变量MR结果格式...\n")
  cat(sprintf("    - 结果总数: %d\n", nrow(univariable_results)))
  cat(sprintf("    - 暴露列名: %s\n", paste(names(univariable_results), collapse = ", ")))
  if (nrow(univariable_results) > 0) {
    cat(sprintf("    - 暴露示例: %s\n", paste(unique(univariable_results$exposure[seq_len(min(5, nrow(univariable_results)))]), collapse = ", ")))
    cat(sprintf("    - 结局示例: %s\n", paste(unique(univariable_results$outcome[seq_len(min(5, nrow(univariable_results)))]), collapse = ", ")))
  }
  
  total_checked <- 0
  passed_count <- 0
  
  for (exp_name in names(metabolic_traits)) {
    for (med_name in names(inflammatory_traits)) {
      # 预筛检查
      for (outcome_name in names(outcomes)) {
        total_checked <- total_checked + 1
        
        # 检查是否有相关的MR结果
        outcome_csv_name <- outcome_name_to_csv[[outcome_name]]
        
        # 检查暴露->结局
        exp_out_exists <- nrow(univariable_results %>%
          filter(exposure == exp_name, outcome == outcome_csv_name)) > 0
        
        # 检查暴露->中介（可能不存在）
        exp_med_exists <- nrow(univariable_results %>%
          filter(exposure == exp_name, outcome == med_name)) > 0
        
        # 检查中介->结局
        med_out_exists <- nrow(univariable_results %>%
          filter(exposure == med_name, outcome == outcome_csv_name)) > 0
        
        if (total_checked <= 5) {  # 前5个打印详细信息
          cat(sprintf("    检查路径 %d: %s -> %s -> %s\n", total_checked, exp_name, med_name, outcome_name))
          cat(sprintf("      暴露->结局: %s\n", ifelse(exp_out_exists, "存在", "不存在")))
          cat(sprintf("      暴露->中介: %s\n", ifelse(exp_med_exists, "存在", "不存在")))
          cat(sprintf("      中介->结局: %s\n", ifelse(med_out_exists, "存在", "不存在")))
        }
        
        if (!pre_screen_pathway(exp_name, med_name, outcome_name)) {
          next
        }
        
        passed_count <- passed_count + 1
        exhaustive_screening_paths[[length(exhaustive_screening_paths) + 1]] <- list(
          source = "exhaustive_screening",
          exposure_name = exp_name,
          exposure_id = metabolic_traits[[exp_name]],
          mediator_name = med_name,
          mediator_id = inflammatory_traits[[med_name]],
          outcome_name = outcome_name,
          outcome_id = outcomes[[outcome_name]],
          rationale = "Statistical significance based screening"
        )
      }
    }
  }
  
  cat(sprintf("    总共检查: %d 条路径, 通过预筛: %d 条\n", total_checked, passed_count))
} else {
  cat("  ⚠ 警告: 单变量MR结果为空，无法进行穷举+预筛\n")
}

cat(sprintf("✓ 智能选择路径: %d 条\n", length(smart_selection_paths)))
cat(sprintf("✓ 穷举+预筛路径: %d 条\n\n", length(exhaustive_screening_paths)))

# ============================================================================
# 步骤5：合并所有路径
# ============================================================================

cat("【步骤5】合并所有中介路径\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

all_pathways <- list()
pathway_counter <- 0

# 添加智能选择路径
for (core in smart_selection_paths) {
  exp_name <- core$exp
  if (!exp_name %in% names(metabolic_traits)) next
  
  exp_id <- metabolic_traits[[exp_name]]
  
  for (med_name in core$med) {
    if (!med_name %in% names(inflammatory_traits)) next
    
    med_id <- inflammatory_traits[[med_name]]
    
    for (outcome_name in names(outcomes)) {
      outcome_id <- outcomes[[outcome_name]]
      
      pathway_counter <- pathway_counter + 1
      
      all_pathways[[pathway_counter]] <- list(
        pathway_id = pathway_counter,
        source = "smart_selection",
        exposure_name = exp_name,
        exposure_id = exp_id,
        mediator_name = med_name,
        mediator_id = med_id,
        outcome_name = outcome_name,
        outcome_id = outcome_id,
        rationale = core$rationale
      )
    }
  }
}

# 添加穷举+预筛路径
for (pathway in exhaustive_screening_paths) {
  pathway_counter <- pathway_counter + 1
  pathway$pathway_id <- pathway_counter
  all_pathways[[pathway_counter]] <- pathway
}

cat(sprintf("✓ 总共定义了 %d 条中介路径\n", length(all_pathways)))
cat(sprintf("  - 智能选择: %d 条\n", 
           sum(sapply(all_pathways, function(x) x$source == "smart_selection"))))
cat(sprintf("  - 穷举+预筛: %d 条\n\n", 
           sum(sapply(all_pathways, function(x) x$source == "exhaustive_screening"))))

# ============================================================================
# 步骤6：定义中介分析函数（增强版，带重试机制）
# ============================================================================

cat("【步骤6】定义中介分析函数\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# MR分析函数（带重试）
perform_single_mr_robust <- function(exposure_data, outcome_id, outcome_name, 
                                     max_retries = 3, batch_size = 50) {
  
  for (attempt in 1:max_retries) {
    result <- tryCatch({
      # 提取结局数据
      outcome_data <- tryCatch({
        extract_outcome_data(snps = exposure_data$SNP, outcomes = outcome_id)
      }, error = function(e) {
        # 如果一次性提取失败，尝试分批提取
        snp_list <- exposure_data$SNP
        n_batches <- ceiling(length(snp_list) / batch_size)
        all_data <- list()
        
        for (i in seq_len(n_batches)) {
          start_idx <- (i - 1) * batch_size + 1
          end_idx <- min(i * batch_size, length(snp_list))
          batch_snps <- snp_list[start_idx:end_idx]
          
          batch_data <- tryCatch({
            extract_outcome_data(snps = batch_snps, outcomes = outcome_id)
          }, error = function(e2) {
            NULL
          })
          
          if (!is.null(batch_data) && nrow(batch_data) > 0) {
            all_data[[i]] <- batch_data
          }
          
          Sys.sleep(1)
        }
        
        if (length(all_data) > 0) {
          do.call(rbind, all_data)
        } else {
          NULL
        }
      })
      
      if (is.null(outcome_data) || nrow(outcome_data) == 0) {
        if (attempt < max_retries) {
          Sys.sleep(3)
          return(NULL)
        } else {
          return(list(error = "无法提取结局数据"))
        }
      }
      
      # 数据协调
      harmonized <- tryCatch({
        harmonise_data(exposure_data, outcome_data, action = 2)
      }, error = function(e) {
        common_snps <- intersect(exposure_data$SNP, outcome_data$SNP)
        if (length(common_snps) == 0) {
          stop("没有共同的SNP")
        }
        exp_filtered <- exposure_data[exposure_data$SNP %in% common_snps, ]
        out_filtered <- outcome_data[outcome_data$SNP %in% common_snps, ]
        harmonise_data(exp_filtered, out_filtered, action = 2)
      })
      
      if (is.null(harmonized) || nrow(harmonized) == 0) {
        if (attempt < max_retries) {
          Sys.sleep(3)
          return(NULL)
        } else {
          return(list(error = "数据协调后无可用SNP"))
        }
      }
      
      # 移除回文SNP
      harmonized <- harmonized[harmonized$palindromic == FALSE | 
                              harmonized$ambiguous == FALSE, ]
      
      if (nrow(harmonized) < 3) {
        return(list(error = sprintf("协调后SNP数量不足 (n=%d)", nrow(harmonized))))
      }
      
      # MR分析
      mr_results <- mr(harmonized, method_list = c("mr_ivw", "mr_egger_regression", 
                                                   "mr_weighted_median"))
      
      heterogeneity <- tryCatch({
        mr_heterogeneity(harmonized)
      }, error = function(e) {
        data.frame(method = "IVW", Q = NA, Q_df = NA, Q_pval = NA)
      })
      
      pleiotropy <- tryCatch({
        mr_pleiotropy_test(harmonized)
      }, error = function(e) {
        data.frame(egger_intercept = NA, se = NA, pval = NA)
      })
      
      return(list(
        harmonized = harmonized,
        mr_results = mr_results,
        heterogeneity = heterogeneity,
        pleiotropy = pleiotropy,
        n_snps = nrow(harmonized)
      ))
      
    }, error = function(e) {
      if (attempt < max_retries) {
        Sys.sleep(3)
        return(NULL)
      } else {
        return(list(error = sprintf("重试%d次后仍失败: %s", 
                                   max_retries, conditionMessage(e))))
      }
    })
    
    if (!is.null(result) && !"error" %in% names(result)) {
      return(result)
    }
  }
  
  return(list(error = "所有重试尝试均失败"))
}

# 完整的中介分析函数
perform_mediation_analysis_robust <- function(pathway) {
  source_label <- ifelse(is.null(pathway$source), "未知", pathway$source)
  
  cat(sprintf("\n  路径 %d: %s -> %s -> %s (%s)\n",
             pathway$pathway_id,
             pathway$exposure_name,
             pathway$mediator_name,
             pathway$outcome_name,
             source_label))
  
  result <- list(
    pathway = pathway,
    exp_to_med = NULL,
    med_to_out = NULL,
    exp_to_out = NULL,
    indirect_effect = NA,
    direct_effect = NA,
    total_effect = NA,
    mediation_proportion = NA,
    partial_success = FALSE,
    error_message = NA
  )
  
  # 获取映射名称
  mapped_exp_name <- exposure_name_mapping[[pathway$exposure_name]]
  mapped_med_name <- exposure_name_mapping[[pathway$mediator_name]]
  
  # 步骤1: 暴露 -> 中介
  cat("    (1) 暴露 -> 中介...")
  exp_instruments <- get_instruments_cached(
    pathway$exposure_name,
    pathway$exposure_id,
    mapped_exp_name,
    "metabolic"
  )
  
  if (is.null(exp_instruments) || nrow(exp_instruments) < 3) {
    n_snps <- ifelse(is.null(exp_instruments), 0, nrow(exp_instruments))
    cat(sprintf(" 失败 (SNP数量: %d, 需要≥3)\n", n_snps))
    result$error_message <- sprintf("无法获取暴露工具变量 (仅获取%d个SNP)", n_snps)
    return(result)
  }
  
  exp_to_med <- perform_single_mr_robust(exp_instruments, pathway$mediator_id, 
                                        pathway$mediator_name)
  
  if (is.null(exp_to_med) || "error" %in% names(exp_to_med)) {
    error_msg <- if("error" %in% names(exp_to_med)) {
      paste("暴露->中介:", exp_to_med$error)
    } else {
      "暴露->中介: 分析失败"
    }
    cat(sprintf(" 失败 (%s)\n", error_msg))
    result$error_message <- error_msg
    return(result)
  }
  
  result$exp_to_med <- exp_to_med
  cat(" 成功\n")
  
  # 步骤2: 中介 -> 结局
  cat("    (2) 中介 -> 结局...")
  med_instruments <- get_instruments_cached(
    pathway$mediator_name,
    pathway$mediator_id,
    mapped_med_name,
    "inflammatory"
  )
  
  if (is.null(med_instruments) || nrow(med_instruments) < 3) {
    n_snps <- ifelse(is.null(med_instruments), 0, nrow(med_instruments))
    cat(sprintf(" 失败 (SNP数量: %d, 需要≥3)\n", n_snps))
    result$error_message <- sprintf("无法获取中介工具变量 (仅获取%d个SNP)", n_snps)
    return(result)
  }
  
  med_to_out <- perform_single_mr_robust(med_instruments, pathway$outcome_id, 
                                        pathway$outcome_name)
  
  if (is.null(med_to_out) || "error" %in% names(med_to_out)) {
    error_msg <- if("error" %in% names(med_to_out)) {
      paste("中介->结局:", med_to_out$error)
    } else {
      "中介->结局: 分析失败"
    }
    cat(sprintf(" 失败 (%s)\n", error_msg))
    result$error_message <- error_msg
    return(result)
  }
  
  result$med_to_out <- med_to_out
  cat(" 成功\n")
  
  # 步骤3: 暴露 -> 结局（总效应）
  cat("    (3) 暴露 -> 结局...")
  exp_to_out <- perform_single_mr_robust(exp_instruments, pathway$outcome_id, 
                                        pathway$outcome_name)
  
  if (is.null(exp_to_out) || "error" %in% names(exp_to_out)) {
    cat(" 失败(部分成功)\n")
    result$error_message <- if("error" %in% names(exp_to_out)) {
      paste("暴露->结局:", exp_to_out$error)
    } else {
      "暴露->结局: 分析失败"
    }
    result$partial_success <- TRUE
    
    # 即使第3步失败，仍计算间接效应
    exp_to_med_ivw <- result$exp_to_med$mr_results %>% 
      filter(method == "Inverse variance weighted")
    med_to_out_ivw <- result$med_to_out$mr_results %>% 
      filter(method == "Inverse variance weighted")
    
    if (nrow(exp_to_med_ivw) > 0 && nrow(med_to_out_ivw) > 0) {
      alpha <- exp_to_med_ivw$b[1]
      beta <- med_to_out_ivw$b[1]
      result$indirect_effect <- alpha * beta
      cat(sprintf("    ⊕ 部分成功: 间接效应 = %.4f (缺少总效应)\n",
                 result$indirect_effect))
    }
    
    return(result)
  }
  
  result$exp_to_out <- exp_to_out
  cat(" 成功\n")
  
  # 计算中介效应
  exp_to_med_ivw <- result$exp_to_med$mr_results %>% 
    filter(method == "Inverse variance weighted")
  med_to_out_ivw <- result$med_to_out$mr_results %>% 
    filter(method == "Inverse variance weighted")
  exp_to_out_ivw <- result$exp_to_out$mr_results %>% 
    filter(method == "Inverse variance weighted")
  
  if (nrow(exp_to_med_ivw) > 0 && nrow(med_to_out_ivw) > 0 && 
      nrow(exp_to_out_ivw) > 0) {
    alpha <- exp_to_med_ivw$b[1]
    beta <- med_to_out_ivw$b[1]
    total <- exp_to_out_ivw$b[1]
    
    result$indirect_effect <- alpha * beta
    result$total_effect <- total
    result$direct_effect <- total - result$indirect_effect
    
    # 中介比例
    if (abs(total) > 1e-10) {
      result$mediation_proportion <- result$indirect_effect / total
    }
    
    cat(sprintf("    ✓ 间接效应 = %.4f, 中介比例 = %.1f%%\n",
               result$indirect_effect,
               ifelse(is.na(result$mediation_proportion), 0, 
                     result$mediation_proportion * 100)))
  }
  
  Sys.sleep(2)  # 避免API限流
  
  return(result)
}

cat("✓ 中介分析函数已定义\n\n")

# ============================================================================
# 步骤7：执行中介分析
# ============================================================================

cat("【步骤7】执行中介分析\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

mediation_results <- list()
mediation_summary <- data.frame()

total_pathways <- length(all_pathways)
completed <- 0
success <- 0
partial_success <- 0
start_time <- Sys.time()

if (total_pathways > 0) {
  pb <- txtProgressBar(min = 0, max = total_pathways, style = 3)
  
  for (i in seq_len(total_pathways)) {
    pathway <- all_pathways[[i]]
    
    setTxtProgressBar(pb, i)
    
    result <- perform_mediation_analysis_robust(pathway)
    
    key <- paste(pathway$exposure_name, pathway$mediator_name, 
                pathway$outcome_name, sep = "_")
    mediation_results[[key]] <- result
    
    # 提取结果用于汇总
    if (!is.null(result$exp_to_med) && !is.null(result$med_to_out)) {
      
      if (result$partial_success) {
        partial_success <- partial_success + 1
      } else if (!is.null(result$exp_to_out)) {
        success <- success + 1
      }
      
      exp_to_med_ivw <- result$exp_to_med$mr_results %>% 
        filter(method == "Inverse variance weighted")
      med_to_out_ivw <- result$med_to_out$mr_results %>% 
        filter(method == "Inverse variance weighted")
      
      if (nrow(exp_to_med_ivw) > 0 && nrow(med_to_out_ivw) > 0) {
        
        alpha <- exp_to_med_ivw$b[1]
        beta <- med_to_out_ivw$b[1]
        se_alpha <- exp_to_med_ivw$se[1]
        se_beta <- med_to_out_ivw$se[1]
        
        # 计算间接效应（确保与result中的值一致）
        indirect_effect_calc <- alpha * beta
        se_indirect <- sqrt((alpha^2 * se_beta^2) + (beta^2 * se_alpha^2))
        z_indirect <- indirect_effect_calc / se_indirect
        p_indirect <- 2 * (1 - pnorm(abs(z_indirect)))
        
        summary_row <- data.frame(
          pathway_id = pathway$pathway_id,
          source = pathway$source,
          exposure = pathway$exposure_name,
          mediator = pathway$mediator_name,
          outcome = pathway$outcome_name,
          rationale = pathway$rationale,
          
          exp_to_med_beta = exp_to_med_ivw$b[1],
          exp_to_med_se = exp_to_med_ivw$se[1],
          exp_to_med_pval = exp_to_med_ivw$pval[1],
          exp_to_med_n_snps = result$exp_to_med$n_snps,
          
          med_to_out_beta = med_to_out_ivw$b[1],
          med_to_out_se = med_to_out_ivw$se[1],
          med_to_out_pval = med_to_out_ivw$pval[1],
          med_to_out_n_snps = result$med_to_out$n_snps,
          
          total_effect = result$total_effect,
          total_effect_beta = if (!is.null(result$exp_to_out)) {
            exp_to_out_ivw <- result$exp_to_out$mr_results %>% 
              filter(method == "Inverse variance weighted")
            if (nrow(exp_to_out_ivw) > 0) exp_to_out_ivw$b[1] else NA
          } else NA,
          total_effect_se = if (!is.null(result$exp_to_out)) {
            exp_to_out_ivw <- result$exp_to_out$mr_results %>% 
              filter(method == "Inverse variance weighted")
            if (nrow(exp_to_out_ivw) > 0) exp_to_out_ivw$se[1] else NA
          } else NA,
          total_effect_pval = if (!is.null(result$exp_to_out)) {
            exp_to_out_ivw <- result$exp_to_out$mr_results %>% 
              filter(method == "Inverse variance weighted")
            if (nrow(exp_to_out_ivw) > 0) exp_to_out_ivw$pval[1] else NA
          } else NA,
          total_effect_n_snps = if (!is.null(result$exp_to_out)) {
            result$exp_to_out$n_snps
          } else NA,
          
          indirect_effect = indirect_effect_calc,
          indirect_effect_se = se_indirect,
          indirect_effect_pval = p_indirect,
          direct_effect = result$direct_effect,
          mediation_proportion = result$mediation_proportion,
          mediation_proportion_percent = ifelse(is.na(result$mediation_proportion), 
                                               NA, 
                                               result$mediation_proportion * 100),
          
          all_paths_significant = exp_to_med_ivw$pval[1] < 0.05 & 
                                 med_to_out_ivw$pval[1] < 0.05,
          indirect_significant = p_indirect < 0.05,
          partial_success = result$partial_success,
          error_message = ifelse(is.na(result$error_message), "", result$error_message),
          
          stringsAsFactors = FALSE
        )
        
        mediation_summary <- rbind(mediation_summary, summary_row)
      }
    }
    
    completed <- completed + 1
    
    # 定期保存进度
    if (completed %% 10 == 0) {
      save(mediation_results, mediation_summary, all_pathways, instrument_cache,
           file = "data/mediation_cache/step07_progress.RData")
      
      elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
      remaining_time <- (elapsed_time / completed) * (total_pathways - completed)
      
      cat(sprintf("\n  💾 进度已保存: %d/%d (%.1f%%), 完全成功: %d, 部分成功: %d, 预计剩余: %.1f分钟\n",
                 completed, total_pathways,
                 100 * completed / total_pathways, success, partial_success,
                 as.numeric(remaining_time)))
    }
    
    # 每20个路径长暂停
    if (completed %% 20 == 0 && completed < total_pathways) {
      cat("\n  ⏸ 暂停30秒以避免API限流...\n")
      Sys.sleep(30)
    }
  }
  
  close(pb)
  cat(sprintf("\n分析完成: 完全成功 %d, 部分成功 %d, 总计 %d/%d\n", 
             success, partial_success, success + partial_success, total_pathways))
} else {
  cat("⚠ 没有定义中介路径，跳过分析\n")
}

# ============================================================================
# 步骤8：多重检验校正和结果保存
# ============================================================================

cat("\n【步骤8】多重检验校正和保存结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (nrow(mediation_summary) > 0) {
  # FDR校正
  mediation_summary$fdr_pval_indirect <- p.adjust(
    mediation_summary$indirect_effect_pval, method = "fdr")
  
  mediation_summary$significant_fdr <- mediation_summary$fdr_pval_indirect < 0.05
  
  # 按显著性和效应量排序
  mediation_summary <- mediation_summary %>%
    arrange(fdr_pval_indirect, desc(abs(indirect_effect)))
  
  # 保存结果
  save(mediation_results, mediation_summary, all_pathways, instrument_cache,
       file = "data/step07_all_mediation_results.RData")
  cat("✓ 已保存: data/step07_all_mediation_results.RData\n")
  
  write.csv(mediation_summary, 
           "results/tables/step07_mediation_results.csv",
           row.names = FALSE)
  cat("✓ 已保存: results/tables/step07_mediation_results.csv\n")
  
  write.xlsx(mediation_summary,
            "results/tables/step07_mediation_results.xlsx",
            rowNames = FALSE)
  cat("✓ 已保存: results/tables/step07_mediation_results.xlsx\n\n")
  
  # 保存更新后的工具变量列表（包含新下载的）
  if (length(all_instruments) > 0) {
    save(all_instruments, file = instruments_file)
    cat("✓ 已保存更新后的工具变量列表\n\n")
  }
} else {
  cat("⚠ 没有成功的中介分析结果\n\n")
}

# ============================================================================
# 步骤9：创建可视化（简化版）
# ============================================================================

cat("【步骤9】创建可视化\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (nrow(mediation_summary) > 0) {
  
  # 9.1 中介效应森林图
  sig_mediation <- mediation_summary %>%
    filter(all_paths_significant == TRUE | indirect_significant == TRUE)
  
  if (nrow(sig_mediation) > 0) {
    sig_mediation <- sig_mediation %>%
      mutate(
        pathway_label = paste0(exposure, " → ", mediator, " → ", outcome),
        or_indirect = exp(indirect_effect),
        or_lci = exp(indirect_effect - 1.96 * indirect_effect_se),
        or_uci = exp(indirect_effect + 1.96 * indirect_effect_se),
        significance = ifelse(significant_fdr, "FDR < 0.05", "P < 0.05"),
        source_label = case_when(
          source == "smart_selection" ~ "智能选择",
          source == "exhaustive_screening" ~ "穷举+预筛",
          TRUE ~ source
        )
      )
    
    p_forest <- ggplot(sig_mediation, 
                      aes(x = or_indirect, 
                          y = reorder(pathway_label, or_indirect),
                          fill = source_label)) +
      geom_point(aes(shape = significance), size = 3) +
      geom_errorbarh(aes(xmin = or_lci, xmax = or_uci),
                    height = 0.3, alpha = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      scale_x_continuous(trans = "log10") +
      scale_shape_manual(values = c("FDR < 0.05" = 16, "P < 0.05" = 1)) +
      scale_fill_manual(values = c("智能选择" = "#4E79A7", "穷举+预筛" = "#F28E2B")) +
      labs(
        title = "Mediation Analysis: Indirect Effects via Inflammatory Markers",
        subtitle = sprintf("%d significant mediation pathways", nrow(sig_mediation)),
        x = "Odds Ratio (95% CI, Indirect Effect)",
        y = "Mediation Pathway",
        shape = "Significance",
        fill = "Strategy"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom"
      )
    
    ggsave("results/figures/step07_mediation_forest_plot.png",
           p_forest, width = 14, height = max(8, nrow(sig_mediation) * 0.3), 
           dpi = 300)
    cat("✓ 已保存: results/figures/step07_mediation_forest_plot.png\n")
  }
  
  # 9.2 中介比例柱状图
  sig_prop <- mediation_summary %>%
    filter(all_paths_significant == TRUE,
           !is.na(mediation_proportion_percent),
           partial_success == FALSE) %>%
    arrange(desc(abs(mediation_proportion_percent))) %>%
    head(20)
  
  if (nrow(sig_prop) > 0) {
    sig_prop$pathway_label <- paste0(
      sig_prop$exposure, " → ", sig_prop$mediator, " → ", sig_prop$outcome
    )
    
    sig_prop <- sig_prop %>%
      mutate(
        source_label = case_when(
          source == "smart_selection" ~ "智能选择",
          source == "exhaustive_screening" ~ "穷举+预筛",
          TRUE ~ source
        )
      )
    
    p_proportion <- ggplot(sig_prop, 
                          aes(x = reorder(pathway_label, mediation_proportion_percent),
                              y = mediation_proportion_percent,
                              fill = source_label)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      scale_fill_manual(values = c("智能选择" = "#4E79A7", "穷举+预筛" = "#F28E2B")) +
      labs(
        title = "Proportion of Effect Mediated by Inflammatory Markers",
        subtitle = "Top 20 pathways with significant indirect effects",
        x = "Mediation Pathway",
        y = "Mediation Proportion (%)",
        fill = "Strategy"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 8)
      )
    
    ggsave("results/figures/step07_mediation_proportion.png",
           p_proportion, width = 12, height = max(8, nrow(sig_prop) * 0.4),
           dpi = 300)
    cat("✓ 已保存: results/figures/step07_mediation_proportion.png\n")
  }
  
  # 9.3 Figure 5 - 整合炎症介质的因果模型图（SCI 10分期刊标准）
  cat("【步骤9.3】创建 Figure 5 - 整合炎症介质的因果模型图\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # 筛选FDR显著的中介路径
  sig_data_fig5 <- mediation_summary %>%
    filter(significant_fdr == TRUE) %>%
    mutate(
      Outcome_label = case_when(
        outcome == "lung_cancer_overall" ~ "Overall\nLung Cancer",
        outcome == "lung_adenocarcinoma" ~ "Lung\nAdenocarcinoma",
        outcome == "squamous_cell_lung" ~ "Squamous Cell\nLung Cancer",
        TRUE ~ outcome
      )
    )
  
  if (nrow(sig_data_fig5) > 0) {
    
    # 提取唯一节点
    unique_exposures <- unique(sig_data_fig5$exposure)
    unique_mediators <- unique(sig_data_fig5$mediator)
    unique_outcomes <- unique(sig_data_fig5$outcome)
    
    # 节点数据框
    node_data <- data.frame(
      name = c(unique_exposures, unique_mediators, unique_outcomes),
      type = c(
        rep("Metabolic", length(unique_exposures)),
        rep("Inflammatory", length(unique_mediators)),
        rep("Outcome", length(unique_outcomes))
      ),
      layer = c(
        rep(1, length(unique_exposures)),      # 暴露层
        rep(2, length(unique_mediators)),       # 中介层
        rep(3, length(unique_outcomes))         # 结局层
      ),
      stringsAsFactors = FALSE
    )
    
    # 计算每层节点数量，用于垂直布局
    n_exposures <- length(unique_exposures)
    n_mediators <- length(unique_mediators)
    n_outcomes <- length(unique_outcomes)
    
    # 为每层节点分配y坐标（垂直均匀分布）
    max_nodes <- max(n_exposures, n_mediators, n_outcomes, na.rm = TRUE)
    y_spacing <- ifelse(max_nodes > 1, 1.0 / (max_nodes - 1), 0)
    
    # 分配y坐标
    node_data <- node_data %>%
      group_by(layer) %>%
      mutate(
        idx = row_number(),
        y = ifelse(n() > 1, 
                   1 - (idx - 1) * y_spacing * (max_nodes - 1) / (n() - 1),
                   0.5)
      ) %>%
      ungroup() %>%
      mutate(
        x = case_when(
          layer == 1 ~ 0.2,    # 暴露层（左侧）
          layer == 2 ~ 0.5,     # 中介层（中间）
          layer == 3 ~ 0.8      # 结局层（右侧）
        )
      )
    
    # 创建边数据
    # 1. 间接路径：暴露 → 中介 → 结局
    # 确保 mediation_proportion_percent 存在，如果不存在则从 mediation_proportion 计算
    if (!"mediation_proportion_percent" %in% names(sig_data_fig5)) {
      sig_data_fig5$mediation_proportion_percent <- ifelse(
        !is.na(sig_data_fig5$mediation_proportion),
        sig_data_fig5$mediation_proportion * 100,
        NA_real_
      )
    }
    
    indirect_edges <- sig_data_fig5 %>%
      select(exposure, mediator, outcome, mediation_proportion_percent, 
             indirect_effect, indirect_effect_pval) %>%
      mutate(
        path_type = "indirect"
      ) %>%
      rename(
        from = exposure,
        to = mediator
      ) %>%
      select(from, to, mediation_proportion_percent, path_type, 
             indirect_effect, indirect_effect_pval)
    
    indirect_edges_part2 <- sig_data_fig5 %>%
      select(mediator, outcome, mediation_proportion_percent,
             indirect_effect, indirect_effect_pval) %>%
      mutate(
        path_type = "indirect"
      ) %>%
      rename(
        from = mediator,
        to = outcome
      ) %>%
      select(from, to, mediation_proportion_percent, path_type,
             indirect_effect, indirect_effect_pval)
    
    # 2. 直接路径：暴露 → 结局（如果存在direct_effect列）
    direct_edges <- data.frame()
    if ("direct_effect" %in% names(sig_data_fig5) && 
        "direct_effect_pval" %in% names(sig_data_fig5)) {
      
      # 检查是否有显著的直接效应
      direct_edges <- sig_data_fig5 %>%
        filter(!is.na(direct_effect), 
               !is.na(direct_effect_pval),
               direct_effect_pval < 0.05) %>%
        select(exposure, outcome, direct_effect, direct_effect_pval) %>%
        mutate(
          path_type = "direct",
          mediation_proportion_percent = 0  # 直接路径没有中介比例
        ) %>%
        rename(
          from = exposure,
          to = outcome
        ) %>%
        select(from, to, mediation_proportion_percent, path_type,
               indirect_effect = direct_effect, 
               indirect_effect_pval = direct_effect_pval)
    }
    
    # 合并所有边
    all_edges <- bind_rows(indirect_edges, indirect_edges_part2, direct_edges)
    
    # 连接节点坐标
    edge_data <- all_edges %>%
      left_join(node_data %>% select(name, from_x = x, from_y = y), 
                by = c("from" = "name")) %>%
      left_join(node_data %>% select(name, to_x = x, to_y = y), 
                by = c("to" = "name"))
    
    # 根据路径类型设置颜色和线型
    edge_data <- edge_data %>%
      mutate(
        edge_color = case_when(
          path_type == "indirect" ~ "#E69F00",  # 橙色 - 间接路径
          path_type == "direct" ~ "#0072B2",     # 蓝色 - 直接路径
          TRUE ~ "#999999"
        ),
        edge_linetype = case_when(
          path_type == "indirect" ~ "solid",
          path_type == "direct" ~ "dashed",
          TRUE ~ "solid"
        ),
        edge_width = case_when(
          path_type == "indirect" ~ pmax(0.8, log10(pmax(0.1, abs(mediation_proportion_percent)) + 1) * 1.2),
          path_type == "direct" ~ 0.6,
          TRUE ~ 0.5
        )
      )
    
    # Okabe-Ito配色方案（色盲友好）
    node_colors <- c(
      "Metabolic" = "#E69F00",      # 橙色 - 代谢性状
      "Inflammatory" = "#56B4E9",    # 天蓝色 - 炎症标志物
      "Outcome" = "#009E73"          # 绿色 - 结局
    )
    
    # 分离间接路径和直接路径数据
    indirect_edge_data <- edge_data %>% filter(path_type == "indirect")
    direct_edge_data <- edge_data %>% filter(path_type == "direct")
    
    # 创建图表（符合SCI 10分期刊标准）
    p_figure5 <- ggplot() +
      # 绘制间接路径（橙色实线）
      {if (nrow(indirect_edge_data) > 0)
        geom_segment(
          data = indirect_edge_data,
          aes(x = from_x, y = from_y, xend = to_x, yend = to_y, linewidth = edge_width),
          color = "#E69F00",
          alpha = 0.6,
          lineend = "round",
          arrow = arrow(length = unit(2.5, "mm"), type = "closed", 
                       ends = "last")
        )
      } +
      # 绘制直接路径（蓝色虚线）
      {if (nrow(direct_edge_data) > 0)
        geom_segment(
          data = direct_edge_data,
          aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
          color = "#0072B2",
          linewidth = 0.6,
          alpha = 0.7,
          linetype = "dashed",
          lineend = "round",
          arrow = arrow(length = unit(2.5, "mm"), type = "closed", 
                       ends = "last")
        )
      } +
      scale_linewidth_continuous(range = c(0.8, 2.5), guide = "none") +
      # 绘制节点
      geom_point(
        data = node_data,
        aes(x = x, y = y, fill = type, size = layer),
        shape = 21,
        color = "white",
        stroke = 0.8
      ) +
      # 节点标签
      geom_text(
        data = node_data,
        aes(x = x, y = y, label = name),
        size = 2.5,  # 约8 pt
        fontface = "bold",
        vjust = ifelse(node_data$layer == 2, 2, ifelse(node_data$layer == 1, -1.5, -1.5)),
        hjust = 0.5,
        color = "black"
      ) +
      # 颜色映射
      scale_fill_manual(
        name = "Variable Type",
        values = node_colors,
        labels = c("Metabolic Trait", "Inflammatory Marker", "Lung Cancer")
      ) +
      # 节点大小（根据层级）
      scale_size_continuous(
        name = "Layer",
        range = c(8, 12),
        guide = "none"
      ) +
      # 坐标轴和主题
      xlim(0, 1) +
      ylim(-0.1, 1.1) +
      labs(
        title = "Figure 5. Causal Model Integrating Inflammatory Mediators",
        subtitle = "Mendelian Randomization Pathways (FDR < 0.05)",
        caption = "Solid orange lines: indirect effects through inflammatory mediators; Dashed blue lines: direct effects"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(
          size = 11,  # 11 pt
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 5, unit = "mm")
        ),
        plot.subtitle = element_text(
          size = 9,  # 9 pt
          hjust = 0.5,
          margin = margin(b = 10, unit = "mm")
        ),
        plot.caption = element_text(
          size = 7,  # 7 pt
          hjust = 0.5,
          color = "gray40",
          margin = margin(t = 10, unit = "mm")
        ),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),  # 9 pt
        legend.text = element_text(size = 8),  # 8 pt
        legend.margin = margin(t = 10, unit = "mm"),
        plot.margin = margin(15, 15, 20, 15, unit = "mm")
      ) +
      coord_fixed()
    
    # 保存图形（SCI 10分期刊标准：双栏宽度174mm）
    ggsave(
      "results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.png",
      p_figure5,
      width = 174,  # 双栏宽度（mm）
      height = 120,  # 高度（mm）
      units = "mm",
      dpi = 600,  # 高分辨率
      bg = "white"
    )
    
    ggsave(
      "results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.pdf",
      p_figure5,
      width = 174,  # 双栏宽度（mm）
      height = 120,  # 高度（mm）
      units = "mm",
      device = "pdf",
      useDingbats = FALSE,
      bg = "white"
    )
    
    cat("✓ Figure 5已保存\n")
    cat("  - PNG: results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.png (174×120 mm, 600 DPI)\n")
    cat("  - PDF: results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.pdf (174×120 mm, 矢量)\n")
    n_indirect <- nrow(filter(all_edges, path_type == "indirect")) / 2  # 除以2因为是两段路径
    cat(sprintf("  - 显示路径: %d 条间接路径", n_indirect))
    if (nrow(filter(all_edges, path_type == "direct")) > 0) {
      cat(sprintf(", %d 条直接路径", nrow(filter(all_edges, path_type == "direct"))))
    }
    cat("\n")
    
  } else {
    cat("⚠ 没有FDR显著的中介路径，跳过Figure 5的生成\n")
  }
}

cat("\n")

# ============================================================================
# 步骤10：生成分析报告
# ============================================================================

cat("【步骤10】生成分析报告\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (nrow(mediation_summary) > 0) {
  
  # 统计汇总
  stats_summary <- list(
    total_pathways = total_pathways,
    successful_analyses = success,
    partial_success_analyses = partial_success,
    success_rate = round(100 * (success + partial_success) / total_pathways, 1),
    
    by_strategy = mediation_summary %>%
      group_by(source) %>%
      summarise(
        n_total = n(),
        n_sig_fdr = sum(significant_fdr, na.rm = TRUE),
        n_sig_nominal = sum(indirect_significant, na.rm = TRUE),
        prop_sig_fdr = round(100 * n_sig_fdr / n_total, 1),
        mean_indirect_effect = mean(indirect_effect, na.rm = TRUE),
        .groups = "drop"
      ),
    
    by_mediator = mediation_summary %>%
      group_by(mediator) %>%
      summarise(
        n_pathways = n(),
        n_significant = sum(significant_fdr, na.rm = TRUE),
        mean_indirect_effect = mean(indirect_effect, na.rm = TRUE),
        .groups = "drop"
      ),
    
    by_exposure = mediation_summary %>%
      group_by(exposure) %>%
      summarise(
        n_pathways = n(),
        n_significant = sum(significant_fdr, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_significant)),
    
    by_outcome = mediation_summary %>%
      group_by(outcome) %>%
      summarise(
        n_pathways = n(),
        n_significant = sum(significant_fdr, na.rm = TRUE),
        .groups = "drop"
      )
  )
  
  # 生成报告
  capture.output({
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("步骤7：中介分析详细报告（完整整合版）\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    
    cat("【分析概况】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    cat(sprintf("定义的中介路径数:           %d\n", stats_summary$total_pathways))
    cat(sprintf("完全成功的分析:             %d\n", stats_summary$successful_analyses))
    cat(sprintf("部分成功的分析:             %d\n", stats_summary$partial_success_analyses))
    cat(sprintf("总成功率:                   %.1f%%\n\n", stats_summary$success_rate))
    
    cat("【策略对比统计】\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    for (i in seq_len(nrow(stats_summary$by_strategy))) {
      row <- stats_summary$by_strategy[i, ]
      strategy_label <- ifelse(row$source == "smart_selection", "智能选择", "穷举+预筛")
      cat(sprintf("%-12s: %3d 路径, %2d FDR显著, %2d 名义显著 (%.1f%%)\n",
                 strategy_label,
                 row$n_total,
                 row$n_sig_fdr,
                 row$n_sig_nominal,
                 row$prop_sig_fdr))
    }
    cat("\n")
    
    # FDR显著的中介路径详情
    n_sig_fdr <- sum(mediation_summary$significant_fdr, na.rm = TRUE)
    if (n_sig_fdr > 0) {
      cat("【FDR显著的中介路径（FDR<0.05）】\n")
      cat(paste(rep("=", 80), collapse = ""), "\n\n")
      
      sig_pathways <- mediation_summary %>%
        filter(significant_fdr == TRUE) %>%
        arrange(fdr_pval_indirect)
      
      for (i in seq_len(min(20, nrow(sig_pathways)))) {
        pathway <- sig_pathways[i, ]
        strategy_label <- ifelse(pathway$source == "smart_selection", "智能选择", "穷举+预筛")
        
        cat(sprintf("%d. %s → %s → %s (%s)\n",
                   i,
                   pathway$exposure,
                   pathway$mediator,
                   pathway$outcome,
                   strategy_label))
        cat(sprintf("   生物学依据: %s\n", pathway$rationale))
        cat(sprintf("   暴露 → 中介:  β=%.4f, P=%.2e (%d SNPs)\n",
                   pathway$exp_to_med_beta,
                   pathway$exp_to_med_pval,
                   pathway$exp_to_med_n_snps))
        cat(sprintf("   中介 → 结局:  β=%.4f, P=%.2e (%d SNPs)\n",
                   pathway$med_to_out_beta,
                   pathway$med_to_out_pval,
                   pathway$med_to_out_n_snps))
        if (!is.na(pathway$total_effect_beta)) {
          cat(sprintf("   总效应:       β=%.4f, P=%.2e (%d SNPs)\n",
                     pathway$total_effect_beta,
                     pathway$total_effect_pval,
                     pathway$total_effect_n_snps))
        }
        cat(sprintf("   间接效应:     β=%.4f, P=%.2e, FDR=%.2e\n",
                   pathway$indirect_effect,
                   pathway$indirect_effect_pval,
                   pathway$fdr_pval_indirect))
        if (!is.na(pathway$mediation_proportion_percent)) {
          cat(sprintf("   中介比例:     %.1f%%\n\n",
                     pathway$mediation_proportion_percent))
        } else {
          cat("\n")
        }
      }
    }
    
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("报告生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    
  }, file = "results/step07_mediation_report.txt")
  
  cat("✓ 已保存: results/step07_mediation_report.txt\n\n")
  
  save(stats_summary, file = "data/step07_stats_summary.RData")
  cat("✓ 已保存: data/step07_stats_summary.RData\n\n")
}

# ============================================================================
# 步骤11：精确亚型分析（腺癌 vs 鳞癌特异性比较）
# ============================================================================

cat("\n【步骤11】精确亚型分析（腺癌 vs 鳞癌特异性比较）\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

if (nrow(mediation_summary) > 0) {
  
  # 创建亚型特异性输出目录
  dir.create("results/tables/step07_subtype_specific", showWarnings = FALSE, recursive = TRUE)
  dir.create("results/figures/step07_subtype_specific", showWarnings = FALSE, recursive = TRUE)
  
  # 亚型名称映射
  subtype_mapping <- data.frame(
    outcome = c("lung_cancer_overall", "lung_adenocarcinoma", "squamous_cell_lung"),
    subtype_label = c("Overall Lung Cancer", "Lung Adenocarcinoma", "Squamous Cell Lung Cancer"),
    subtype_short = c("Overall", "Adenocarcinoma", "Squamous"),
    stringsAsFactors = FALSE
  )
  
  # 按亚型分组汇总
  mediation_by_subtype <- mediation_summary %>%
    left_join(subtype_mapping, by = "outcome") %>%
    group_by(subtype_short)
  
  # 生成每个亚型的独立汇总表
  subtype_summaries <- list()
  
  for (subtype in subtype_mapping$subtype_short) {
    subtype_data <- mediation_summary %>%
      filter(outcome == subtype_mapping$outcome[subtype_mapping$subtype_short == subtype])
    
    if (nrow(subtype_data) > 0) {
      # 亚型内FDR校正
      subtype_data$fdr_pval_indirect_subtype <- p.adjust(
        subtype_data$indirect_effect_pval, method = "fdr"
      )
      subtype_data$significant_fdr_subtype <- subtype_data$fdr_pval_indirect_subtype < 0.05
      
      # 按显著性排序
      subtype_data <- subtype_data %>%
        arrange(fdr_pval_indirect_subtype, desc(abs(indirect_effect)))
      
      subtype_summaries[[subtype]] <- subtype_data
      
      # 保存亚型特异性结果
      write.xlsx(
        subtype_data,
        sprintf("results/tables/step07_subtype_specific/mediation_results_%s.xlsx", subtype),
        rowNames = FALSE
      )
      
      cat(sprintf("✓ %s: %d 条路径, %d 个FDR显著\n",
                 subtype,
                 nrow(subtype_data),
                 sum(subtype_data$significant_fdr_subtype, na.rm = TRUE)))
    }
  }
  
  # 生成亚型比较表（比较相同路径在不同亚型中的效应）
  cat("\n  生成亚型比较表...\n")
  
  # 提取所有唯一的路径（暴露-中介组合，不考虑结局）
  unique_pathways <- mediation_summary %>%
    select(exposure, mediator, rationale) %>%
    distinct()
  
  # 创建亚型比较矩阵
  subtype_comparison <- data.frame()
  
  for (i in seq_len(nrow(unique_pathways))) {
    pathway <- unique_pathways[i, ]
    
    # 提取该路径在三个亚型中的结果
    pathway_results <- mediation_summary %>%
      filter(exposure == pathway$exposure,
             mediator == pathway$mediator)
    
    # 总体肺癌结果
    overall_result <- pathway_results %>%
      filter(outcome == "lung_cancer_overall")
    
    # 腺癌结果
    adeno_result <- pathway_results %>%
      filter(outcome == "lung_adenocarcinoma")
    
    # 鳞癌结果
    squamous_result <- pathway_results %>%
      filter(outcome == "squamous_cell_lung")
    
    # 构建比较行
    comparison_row <- data.frame(
      exposure = pathway$exposure,
      mediator = pathway$mediator,
      rationale = pathway$rationale,
      
      # 总体肺癌
      overall_indirect = if(nrow(overall_result) > 0) overall_result$indirect_effect[1] else NA,
      overall_pval = if(nrow(overall_result) > 0) overall_result$indirect_effect_pval[1] else NA,
      overall_fdr = if(nrow(overall_result) > 0) overall_result$fdr_pval_indirect[1] else NA,
      overall_mediation_pct = if(nrow(overall_result) > 0) overall_result$mediation_proportion_percent[1] else NA,
      
      # 腺癌
      adeno_indirect = if(nrow(adeno_result) > 0) adeno_result$indirect_effect[1] else NA,
      adeno_pval = if(nrow(adeno_result) > 0) adeno_result$indirect_effect_pval[1] else NA,
      adeno_fdr = if(nrow(adeno_result) > 0) {
        if("fdr_pval_indirect" %in% names(adeno_result)) {
          adeno_result$fdr_pval_indirect[1]
        } else {
          p.adjust(adeno_result$indirect_effect_pval, method = "fdr")[1]
        }
      } else NA,
      adeno_mediation_pct = if(nrow(adeno_result) > 0) adeno_result$mediation_proportion_percent[1] else NA,
      
      # 鳞癌
      squamous_indirect = if(nrow(squamous_result) > 0) squamous_result$indirect_effect[1] else NA,
      squamous_pval = if(nrow(squamous_result) > 0) squamous_result$indirect_effect_pval[1] else NA,
      squamous_fdr = if(nrow(squamous_result) > 0) {
        if("fdr_pval_indirect" %in% names(squamous_result)) {
          squamous_result$fdr_pval_indirect[1]
        } else {
          p.adjust(squamous_result$indirect_effect_pval, method = "fdr")[1]
        }
      } else NA,
      squamous_mediation_pct = if(nrow(squamous_result) > 0) squamous_result$mediation_proportion_percent[1] else NA,
      
      stringsAsFactors = FALSE
    )
    
    # 判断亚型特异性（仅在某一亚型显著）
    comparison_row$adeno_specific <- !is.na(comparison_row$adeno_fdr) && 
                                    comparison_row$adeno_fdr < 0.05 &&
                                    (is.na(comparison_row$squamous_fdr) || comparison_row$squamous_fdr >= 0.05)
    
    comparison_row$squamous_specific <- !is.na(comparison_row$squamous_fdr) && 
                                       comparison_row$squamous_fdr < 0.05 &&
                                       (is.na(comparison_row$adeno_fdr) || comparison_row$adeno_fdr >= 0.05)
    
    comparison_row$both_significant <- !is.na(comparison_row$adeno_fdr) && 
                                       !is.na(comparison_row$squamous_fdr) &&
                                       comparison_row$adeno_fdr < 0.05 && 
                                       comparison_row$squamous_fdr < 0.05
    
    subtype_comparison <- rbind(subtype_comparison, comparison_row)
  }
  
  # 保存亚型比较表
  if (nrow(subtype_comparison) > 0) {
    write.xlsx(
      subtype_comparison,
      "results/tables/step07_subtype_specific/subtype_comparison_table.xlsx",
      rowNames = FALSE
    )
    
    write.csv(
      subtype_comparison,
      "results/tables/step07_subtype_specific/subtype_comparison_table.csv",
      row.names = FALSE
    )
    
    cat(sprintf("✓ 亚型比较表已保存: %d 条路径\n", nrow(subtype_comparison)))
    
    # 统计亚型特异性路径
    n_adeno_specific <- sum(subtype_comparison$adeno_specific, na.rm = TRUE)
    n_squamous_specific <- sum(subtype_comparison$squamous_specific, na.rm = TRUE)
    n_both_significant <- sum(subtype_comparison$both_significant, na.rm = TRUE)
    
    cat(sprintf("  - 腺癌特异性路径: %d 条\n", n_adeno_specific))
    cat(sprintf("  - 鳞癌特异性路径: %d 条\n", n_squamous_specific))
    cat(sprintf("  - 两亚型均显著: %d 条\n", n_both_significant))
  }
  
  # 生成亚型特异性可视化（如果ggplot2可用）
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    cat("\n  生成亚型特异性可视化...\n")
    
    # 筛选显著路径用于可视化
    sig_comparison <- subtype_comparison %>%
      filter((!is.na(adeno_fdr) & adeno_fdr < 0.05) | 
             (!is.na(squamous_fdr) & squamous_fdr < 0.05)) %>%
      arrange(desc(abs(adeno_indirect) + abs(squamous_indirect)))
    
    if (nrow(sig_comparison) > 0) {
      # 限制显示数量（最多20条）
      sig_comparison <- head(sig_comparison, 20)
      
      # 准备绘图数据
      plot_data <- sig_comparison %>%
        mutate(
          pathway_label = paste0(exposure, " → ", mediator),
          pathway_label = factor(pathway_label, levels = rev(unique(pathway_label)))
        ) %>%
        select(pathway_label, adeno_indirect, squamous_indirect, 
               adeno_fdr, squamous_fdr) %>%
        tidyr::pivot_longer(
          cols = c(adeno_indirect, squamous_indirect),
          names_to = "subtype",
          values_to = "indirect_effect"
        ) %>%
        mutate(
          subtype = ifelse(subtype == "adeno_indirect", "Adenocarcinoma", "Squamous"),
          fdr_pval = ifelse(subtype == "Adenocarcinoma", adeno_fdr, squamous_fdr),
          significant = !is.na(fdr_pval) & fdr_pval < 0.05
        )
      
      # 创建亚型比较图
      p_subtype_comparison <- ggplot(plot_data, 
                                     aes(x = indirect_effect, 
                                         y = pathway_label, 
                                         fill = subtype)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        scale_fill_manual(
          values = c("Adenocarcinoma" = "#E69F00", "Squamous" = "#56B4E9"),
          name = "Lung Cancer Subtype"
        ) +
        labs(
          title = "Subtype-Specific Mediation Effects",
          subtitle = "Comparison of indirect effects between adenocarcinoma and squamous cell carcinoma",
          x = "Indirect Effect (β)",
          y = "Mediation Pathway"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 7),
          legend.position = "bottom"
        )
      
      ggsave(
        "results/figures/step07_subtype_specific/subtype_comparison_plot.png",
        p_subtype_comparison,
        width = 12,
        height = max(8, nrow(sig_comparison) * 0.3),
        dpi = 300
      )
      
      cat("✓ 亚型比较图已保存\n")
    }
  }
  
  # 保存亚型分析结果
  save(subtype_summaries, subtype_comparison,
       file = "data/step07_subtype_analysis.RData")
  cat("✓ 已保存: data/step07_subtype_analysis.RData\n\n")
  
} else {
  cat("⚠ 没有中介分析结果，跳过亚型分析\n\n")
}

# ============================================================================
# 最终总结
# ============================================================================

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤7完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【最终统计】\n")
cat(sprintf("定义的中介路径:     %d 条\n", total_pathways))

if (total_pathways > 0) {
  cat(sprintf("完全成功分析:       %d 条 (%.1f%%)\n", 
             success, 100 * success / total_pathways))
  cat(sprintf("部分成功分析:       %d 条 (%.1f%%)\n",
             partial_success, 100 * partial_success / total_pathways))
  cat(sprintf("总体成功率:         %.1f%%\n",
             100 * (success + partial_success) / total_pathways))
}

if (nrow(mediation_summary) > 0) {
  cat(sprintf("FDR显著中介:        %d 条\n", 
             sum(mediation_summary$significant_fdr, na.rm = TRUE)))
  cat(sprintf("名义显著中介:       %d 条\n\n",
             sum(mediation_summary$indirect_significant, na.rm = TRUE)))
}

cat("【关键改进】\n")
cat("  ✓ 参考step05的数据提取方式，确保数据无偏差\n")
cat("  ✓ 优先使用本地工具变量文件，缺失时自动下载\n")
cat("  ✓ 使用与step05完全一致的GWAS ID映射和提取策略\n")
cat("  ✓ 扩展智能选择路径：26条核心策略（包含所有生物学合理路径）\n")
cat("  ✓ 整合智能选择和穷举+预筛两种策略\n")
cat("  ✓ 增强的错误处理和重试机制\n")
cat("  ✓ 自动保存进度，支持断点续传\n")
cat("  ✓ 精确亚型分析：腺癌 vs 鳞癌特异性比较\n\n")

cat("【保存的文件】\n")
cat("  数据文件:\n")
cat("    - data/step07_all_mediation_results.RData\n")
cat("    - data/step07_stats_summary.RData\n")
cat("    - data/mediation_cache/step07_progress.RData (检查点)\n\n")
cat("  结果表格:\n")
cat("    - results/tables/step07_mediation_results.xlsx\n")
cat("    - results/tables/step07_mediation_results.csv\n\n")
cat("  报告文件:\n")
cat("    - results/step07_mediation_report.txt\n\n")
cat("  可视化图表:\n")
if (file.exists("results/figures/step07_mediation_forest_plot.png")) {
  cat("    - results/figures/step07_mediation_forest_plot.png\n")
}
if (file.exists("results/figures/step07_mediation_proportion.png")) {
  cat("    - results/figures/step07_mediation_proportion.png\n")
}
if (file.exists("results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.png")) {
  cat("    - results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.png (SCI 10分期刊标准)\n")
  cat("    - results/figures/step07_publication/Figure5_Causal_Model_Integrating_Mediators.pdf (SCI 10分期刊标准)\n")
}
if (file.exists("results/figures/step07_subtype_specific/subtype_comparison_plot.png")) {
  cat("    - results/figures/step07_subtype_specific/subtype_comparison_plot.png (亚型特异性比较图)\n")
}
if (file.exists("results/tables/step07_subtype_specific/subtype_comparison_table.xlsx")) {
  cat("    - results/tables/step07_subtype_specific/subtype_comparison_table.xlsx (亚型比较表)\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("中介分析完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

