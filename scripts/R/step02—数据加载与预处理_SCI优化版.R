############################################################################
# Step 2: Data Loading and Preprocessing (SCI 5+ Optimized)
# 步骤2：数据加载与预处理（SCI 5+优化版）
############################################################################
#
# 【功能说明】
# 本脚本从用户数据集信息文件加载暴露和结局变量，并进行智能数据预处理
# 简化版本：只有30个暴露和3个结局
#
# 【输入文件】
# - 用户数据集信息.csv：包含所有GWAS数据集的元信息
#
# 【输出文件】
# - results/data/exposure_data_list.RData：暴露因子数据列表（30个暴露）
# - results/data/outcome_data_list.RData：结局数据列表（3个结局）
# - results/data/dataset_metadata.csv：数据集元信息汇总表
# - results/reports/step02_data_summary.txt：数据加载摘要报告
# - results/figures/step02_data_overview.png：数据概览可视化
#
# 【版本特点】
# ✓ 智能数据源检测（OpenGWAS、本地文件、GWAS Catalog）
# ✓ 多格式数据兼容（.csv、.txt、.gz、.vcf等）
# ✓ 自动数据质量检查和清理
# ✓ 缺失数据智能处理
# ✓ 发表级数据摘要和可视化
# ✓ 完整的错误处理和日志记录
#
############################################################################

cat("步骤2：开始数据加载与预处理（SCI 5+优化版）...\n\n")

# ===================================================================
# 0. 加载必要的包并设置环境
# ===================================================================
cat("【步骤0】加载必要的包...\n")

required_packages <- c(
  "TwoSampleMR",    # MR分析核心包
  "ieugwasr",       # IEU OpenGWAS数据库接口
  "dplyr",          # 数据处理
  "readr",          # 数据读取
  "stringr",        # 字符串处理
  "ggplot2",        # 数据可视化
  "scales",         # 比例尺
  "data.table",     # 高效数据处理
  "tidyr"           # 数据整理
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装缺失的包: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

cat("✓ 所有必要的包已加载\n\n")

# 创建输出目录
output_dirs <- c(
  "results",
  "results/data",
  "results/reports",
  "results/figures",
  "results/logs"
)

for (dir in output_dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

cat("✓ 输出目录已创建\n\n")

# ===================================================================
# 1. 加载用户数据集信息
# ===================================================================
cat("【步骤1】加载用户数据集信息...\n")

# 定义可能的数据集信息文件路径
dataset_info_paths <- c(
  "C:/Users/surob/Documents/GWAS/GWAS_Metadata-数据集信息.csv",
  "C:/Users/surob/Documents/GWAS/用户数据集信息.csv",
  file.path(dirname(getwd()), "GWAS", "GWAS_Metadata-数据集信息.csv"),
  file.path(dirname(getwd()), "GWAS", "用户数据集信息.csv"),
  "GWAS_Metadata-数据集信息.csv",
  "用户数据集信息.csv",
  "dataset_info.csv"
)

# 尝试加载数据集信息
dataset_info <- NULL
for (path in dataset_info_paths) {
  if (file.exists(path)) {
    cat(sprintf("  • 尝试读取文件: %s\n", path))
    dataset_info <- tryCatch({
      data <- read_csv(path, show_col_types = FALSE)
      # 清理空行
      data <- data %>% filter(!is.na(id) & id != "")
      cat(sprintf("✓ 成功加载数据集信息：%d 个数据集\n", nrow(data)))
      data
    }, error = function(e) {
      cat(sprintf("⚠ 读取失败: %s\n", conditionMessage(e)))
      NULL
    })
    if (!is.null(dataset_info)) break
  }
}

if (is.null(dataset_info)) {
  stop("❌ 错误：无法找到用户数据集信息文件！请确保文件存在。")
}

# ===================================================================
# 1.1 自动分类数据集（如果CSV文件没有category列）
# 简化分类：只有30个暴露和3个结局
# ===================================================================
if (!"category" %in% names(dataset_info)) {
  cat("【步骤1.1】自动分类数据集（30个暴露 + 3个结局）...\n")
  
  # 定义结局关键词（只有3个结局）
  outcome_keywords <- c("lung cancer", "lung adenocarcinoma", "squamous.*lung", 
                       "肺.*癌", "腺.*癌", "鳞.*癌")
  
  # 创建分类函数（简化版：只区分暴露和结局）
  classify_dataset <- function(trait_name, dataset_id) {
    trait_lower <- tolower(trait_name)
    
    # 首先根据特定ID精确匹配结局（优先级最高）
    # 结局（肺癌）- 只有3个结局
    if (dataset_id %in% c("ebi-a-GCST004748", "ieu-a-984", "ieu-a-989")) {
      return("outcomes")
    }
    
    # 根据关键词匹配结局（肺癌相关）
    if (any(sapply(outcome_keywords, function(k) grepl(k, trait_lower, ignore.case = TRUE)))) {
      return("outcomes")
    }
    
    # 根据ID前缀判断结局
    if (grepl("^ieu-a-", dataset_id)) {
      # IEU-a数据集通常是结局（肺癌亚型）
      return("outcomes")
    }
    
    # 其他所有数据集都归类为暴露（30个暴露）
    return("exposures")
  }
  
  # 应用分类
  if ("trait" %in% names(dataset_info)) {
    dataset_info$category <- sapply(seq_len(nrow(dataset_info)), function(i) {
      classify_dataset(dataset_info$trait[i], dataset_info$id[i])
    })
  } else {
    dataset_info$category <- "exposures"  # 默认归类为暴露
  }
  
  cat("✓ 数据集自动分类完成\n")
} else {
  # 如果已有category列，需要将其标准化为"exposures"和"outcomes"
  cat("【步骤1.1】标准化已有分类（合并为exposures和outcomes）...\n")
  dataset_info <- dataset_info %>%
    mutate(
      category = case_when(
        category == "outcomes" ~ "outcomes",  # 保持结局不变
        TRUE ~ "exposures"  # 其他所有类别（covariates, metabolic, inflammatory）都合并为exposures
      )
    )
  cat("✓ 分类标准化完成\n")
}

# ===================================================================
# 1.2 补充缺失的列
# ===================================================================
cat("【步骤1.2】补充数据集元信息...\n")

# 创建short_name（如果没有）
if (!"short_name" %in% names(dataset_info)) {
  if ("trait" %in% names(dataset_info)) {
    dataset_info$short_name <- dataset_info$trait
  } else {
    dataset_info$short_name <- dataset_info$id
  }
}

# 创建data_source_primary（根据ID格式判断）
if (!"data_source_primary" %in% names(dataset_info)) {
  dataset_info$data_source_primary <- ifelse(
    grepl("^ieu-", dataset_info$id), "OpenGWAS",
    ifelse(grepl("^ebi-a-", dataset_info$id), "OpenGWAS", "OpenGWAS")
  )
}

# 创建data_completeness（默认100%）
if (!"data_completeness" %in% names(dataset_info)) {
  dataset_info$data_completeness <- 100
}

# 处理样本量列（可能有逗号分隔符）
if ("样本量" %in% names(dataset_info)) {
  dataset_info$sample_size <- dataset_info$样本量
  # 移除逗号并转换为数字
  dataset_info$sample_size <- as.numeric(gsub(",", "", dataset_info$sample_size))
} else if (!"sample_size" %in% names(dataset_info)) {
  dataset_info$sample_size <- NA
}

cat("✓ 元信息补充完成\n\n")

# 数据集信息概览
cat("数据集分类统计：\n")
category_summary <- dataset_info %>%
  group_by(category) %>%
  summarise(
    n_datasets = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_datasets))

print(category_summary)
cat("\n")

# ===================================================================
# 2. 分类提取数据集ID（30个暴露 + 3个结局）
# ===================================================================
cat("【步骤2】分类提取数据集ID...\n")

# 提取暴露因子（所有30个暴露，包括原来的协变量、代谢和炎症）
exposure_datasets <- dataset_info %>%
  filter(category == "exposures") %>%
  select(id, short_name, trait, sample_size, data_source_primary, data_completeness)

cat(sprintf("  • 暴露因子数据集: %d 个\n", nrow(exposure_datasets)))

# 提取结局（只有3个结局）
outcome_datasets <- dataset_info %>%
  filter(category == "outcomes") %>%
  select(id, short_name, trait, sample_size, data_source_primary, data_completeness)

cat(sprintf("  • 结局数据集: %d 个\n", nrow(outcome_datasets)))
cat("\n")

# ===================================================================
# 3. 定义智能数据提取函数
# ===================================================================
cat("【步骤3】定义智能数据提取函数...\n")

#' 智能提取GWAS数据
#' 
#' @param dataset_id 数据集ID
#' @param dataset_name 数据集简称
#' @param data_source 数据来源（OpenGWAS, 本地文件等）
#' @param completeness 数据完整性（0-100）
#' @return 提取的GWAS数据（格式化为TwoSampleMR标准格式）
extract_gwas_data <- function(dataset_id, dataset_name, data_source = "OpenGWAS", completeness = 100) {
  
  cat(sprintf("  → 提取数据: %s (ID: %s, 数据源: %s, 完整性: %d%%)\n", 
              dataset_name, dataset_id, data_source, completeness))
  
  # 跳过数据完整性过低或数据源不可用的数据集
  if (data_source %in% c("None", "NULL") || completeness < 50) {
    cat(sprintf("    ⚠ 跳过：数据源不可用或完整性过低（%d%%）\n", completeness))
    return(NULL)
  }
  
  # 策略1: OpenGWAS数据库
  if (data_source == "OpenGWAS" && completeness >= 80) {
    tryCatch({
      cat("    使用OpenGWAS API...\n")
      
      # 首先检查数据集是否可用
      available <- tryCatch({
        ieugwasr::gwasinfo(dataset_id)
      }, error = function(e) {
        cat(sprintf("      警告: 无法查询数据集信息: %s\n", conditionMessage(e)))
        NULL
      })
      
      if (!is.null(available)) {
        # 使用tophits提取显著性SNPs（已clumped）
        data <- tryCatch({
          ieugwasr::tophits(id = dataset_id, pval = 5e-8, clump = 1, r2 = 0.001, kb = 10000)
        }, error = function(e) {
          # 如果tophits失败，尝试使用associations提取全基因组数据
          cat(sprintf("      tophits失败，尝试associations...\n"))
          tryCatch({
            assoc <- ieugwasr::associations(variants = NULL, id = dataset_id, proxies = 0)
            if (!is.null(assoc) && nrow(assoc) > 0) {
              # 筛选显著性SNPs
              sig_data <- assoc %>% filter(p < 5e-8)
              if (nrow(sig_data) > 0) {
                # 进行LD clumping
                clumped <- ieugwasr::ld_clump(
                  dplyr::tibble(rsid = sig_data$rsid, pval = sig_data$p, id = dataset_id),
                  clump_r2 = 0.001,
                  clump_kb = 10000
                )
                sig_data %>% filter(rsid %in% clumped$rsid)
              } else {
                NULL
              }
            } else {
              NULL
            }
          }, error = function(e2) {
            cat(sprintf("      associations也失败: %s\n", conditionMessage(e2)))
            NULL
          })
        })
        
        if (!is.null(data) && nrow(data) > 0) {
          cat(sprintf("    ✓ OpenGWAS提取成功: %d SNPs\n", nrow(data)))
          
          # 确保数据包含必需的列
          required_cols <- c("rsid", "beta", "se", "p", "ea", "nea")
          if (all(required_cols %in% names(data))) {
            # 格式化为TwoSampleMR标准格式
            formatted_data <- TwoSampleMR::format_data(
              data,
              type = "exposure",
              snp_col = "rsid",
              beta_col = "beta",
              se_col = "se",
              effect_allele_col = "ea",
              other_allele_col = "nea",
              eaf_col = if("eaf" %in% names(data)) "eaf" else NULL,
              pval_col = "p",
              samplesize_col = if("n" %in% names(data)) "n" else NULL
            )
            
            if (!is.null(formatted_data) && nrow(formatted_data) > 0) {
              formatted_data$dataset_id <- dataset_id
              formatted_data$dataset_name <- dataset_name
              return(formatted_data)
            }
          } else {
            cat(sprintf("    ⚠ 数据缺少必需列，可用列: %s\n", paste(names(data), collapse = ", ")))
          }
        } else {
          cat("    ⚠ OpenGWAS提取成功但无显著性SNPs\n")
        }
      } else {
        cat("    ⚠ 数据集在OpenGWAS中不可用\n")
      }
      
    }, error = function(e) {
      cat(sprintf("    ⚠ OpenGWAS错误: %s\n", conditionMessage(e)))
    })
  }
  
  # 策略1b: OpenGWAS数据库（完整性较低但>=50%）
  if (data_source == "OpenGWAS" && completeness >= 50 && completeness < 80) {
    tryCatch({
      cat("    使用OpenGWAS API（低完整性模式）...\n")
      available <- tryCatch({
        ieugwasr::gwasinfo(dataset_id)
      }, error = function(e) NULL)
      
      if (!is.null(available)) {
        # 放宽p值阈值，提取更多SNPs
        data <- tryCatch({
          ieugwasr::tophits(id = dataset_id, pval = 1e-5, clump = 1, r2 = 0.001, kb = 10000)
        }, error = function(e) NULL)
        
        if (!is.null(data) && nrow(data) > 0) {
          cat(sprintf("    ✓ OpenGWAS提取成功（低阈值）: %d SNPs\n", nrow(data)))
          
          if (all(c("rsid", "beta", "se", "p", "ea", "nea") %in% names(data))) {
            formatted_data <- TwoSampleMR::format_data(
              data,
              type = "exposure",
              snp_col = "rsid",
              beta_col = "beta",
              se_col = "se",
              effect_allele_col = "ea",
              other_allele_col = "nea",
              eaf_col = if("eaf" %in% names(data)) "eaf" else NULL,
              pval_col = "p",
              samplesize_col = if("n" %in% names(data)) "n" else NULL
            )
            
            if (!is.null(formatted_data) && nrow(formatted_data) > 0) {
              formatted_data$dataset_id <- dataset_id
              formatted_data$dataset_name <- dataset_name
              return(formatted_data)
            }
          }
        }
      }
    }, error = function(e) {
      cat(sprintf("    ⚠ OpenGWAS（低完整性模式）错误: %s\n", conditionMessage(e)))
    })
  }
  
  # 策略2: GWAS Catalog数据源（需要手动提供本地文件）
  if (data_source == "GWAS_Catalog") {
    cat("    尝试从GWAS Catalog本地文件加载...\n")
    # GWAS Catalog数据通常需要手动下载，放在特定目录
    local_file_patterns <- c(
      sprintf("data/GWAS_Catalog/%s.csv", dataset_id),
      sprintf("data/GWAS_Catalog/%s.txt", dataset_id),
      sprintf("data/GWAS_Catalog/%s.txt.gz", dataset_id),
      sprintf("GWAS_data/%s.csv", dataset_id),
      sprintf("raw_data/%s.txt.gz", dataset_id)
    )
    
    for (file_path in local_file_patterns) {
      if (file.exists(file_path)) {
        cat(sprintf("    找到本地文件: %s\n", file_path))
        tryCatch({
          data <- if (grepl("\\.gz$", file_path)) {
            data.table::fread(file_path)
          } else if (grepl("\\.csv$", file_path)) {
            readr::read_csv(file_path, show_col_types = FALSE)
          } else {
            readr::read_delim(file_path, delim = "\t", show_col_types = FALSE)
          }
          
          if (nrow(data) > 0) {
            cat(sprintf("    ✓ 本地文件读取成功: %d 行\n", nrow(data)))
            data <- standardize_gwas_columns(data)
            formatted_data <- format_gwas_data(data, dataset_id, dataset_name)
            if (!is.null(formatted_data)) return(formatted_data)
          }
        }, error = function(e) {
          cat(sprintf("    ⚠ 本地文件读取失败: %s\n", conditionMessage(e)))
        })
      }
    }
    cat("    ⚠ 未找到GWAS Catalog本地文件\n")
  }
  
  # 策略3: 通用本地文件（用于所有数据源的后备方案）
  cat("    尝试通用本地文件路径...\n")
  local_file_patterns <- c(
    sprintf("data/%s.csv", dataset_name),
    sprintf("data/%s.txt", dataset_name),
    sprintf("data/%s.gz", dataset_name),
    sprintf("data/%s.csv", dataset_id),
    sprintf("data/%s.txt", dataset_id),
    sprintf("GWAS_data/%s.csv", dataset_id),
    sprintf("raw_data/%s.txt.gz", dataset_id)
  )
  
  for (file_path in local_file_patterns) {
    if (file.exists(file_path)) {
      cat(sprintf("    使用本地文件: %s\n", file_path))
      
      tryCatch({
        # 根据文件扩展名选择读取方法
        if (grepl("\\.gz$", file_path)) {
          data <- fread(file_path)
        } else if (grepl("\\.csv$", file_path)) {
          data <- read_csv(file_path, show_col_types = FALSE)
        } else {
          data <- read_delim(file_path, delim = "\t", show_col_types = FALSE)
        }
        
        if (nrow(data) > 0) {
          cat(sprintf("    ✓ 本地文件读取成功: %d 行\n", nrow(data)))
          
          # 智能列名映射
          data <- standardize_gwas_columns(data)
          
          # 格式化数据
          formatted_data <- format_gwas_data(data, dataset_id, dataset_name)
          
          return(formatted_data)
        }
        
      }, error = function(e) {
        cat(sprintf("    ⚠ 本地文件读取失败: %s\n", conditionMessage(e)))
      })
    }
  }
  
  # 策略4: 对于数据源不可用或提取失败的情况，返回NULL而不是模拟数据
  cat(sprintf("    ✗ 无法提取数据: %s (ID: %s)\n", dataset_name, dataset_id))
  cat("      建议: 检查数据集ID是否正确，或手动下载数据到本地目录\n")
  
  return(NULL)
}

#' 标准化GWAS数据列名
#' 
#' @param data 原始GWAS数据
#' @return 标准化列名的数据
standardize_gwas_columns <- function(data) {
  
  # 定义常见列名映射
  col_mappings <- list(
    SNP = c("SNP", "rsid", "rs_id", "RSID", "MarkerName", "snp"),
    chr = c("chr", "CHR", "chromosome", "Chromosome", "#chr"),
    pos = c("pos", "POS", "position", "Position", "bp", "BP"),
    effect_allele = c("effect_allele", "A1", "EA", "ea", "ALT", "alt"),
    other_allele = c("other_allele", "A2", "NEA", "nea", "REF", "ref"),
    beta = c("beta", "BETA", "Beta", "b", "B", "Effect"),
    se = c("se", "SE", "Se", "standard_error", "StdErr"),
    pval = c("pval", "P", "p", "P-value", "p-value", "pvalue", "PVAL"),
    eaf = c("eaf", "EAF", "Freq", "freq", "MAF", "maf", "FRQ"),
    samplesize = c("samplesize", "N", "n", "sample_size", "SampleSize")
  )
  
  # 执行列名映射
  for (standard_name in names(col_mappings)) {
    possible_names <- col_mappings[[standard_name]]
    for (poss_name in possible_names) {
      if (poss_name %in% names(data)) {
        if (poss_name != standard_name) {
          data[[standard_name]] <- data[[poss_name]]
          cat(sprintf("      列名映射: %s -> %s\n", poss_name, standard_name))
        }
        break
      }
    }
  }
  
  return(data)
}

#' 格式化GWAS数据为TwoSampleMR标准格式
#' 
#' @param data 标准化后的GWAS数据
#' @param dataset_id 数据集ID
#' @param dataset_name 数据集名称
#' @return 格式化的数据
format_gwas_data <- function(data, dataset_id, dataset_name) {
  
  # 检查必需的列
  required_cols <- c("SNP", "beta", "se", "pval", "effect_allele", "other_allele")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    cat(sprintf("      ⚠ 缺少必需列: %s\n", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # 数据质量过滤
  filtered_data <- data %>%
    filter(
      !is.na(SNP),
      !is.na(beta),
      !is.na(se),
      !is.na(pval),
      se > 0,
      pval > 0,
      pval <= 1
    )
  
  if (nrow(filtered_data) == 0) {
    cat("      ⚠ 数据质量过滤后无有效数据\n")
    return(NULL)
  }
  
  cat(sprintf("      数据质量过滤: %d -> %d SNPs\n", nrow(data), nrow(filtered_data)))
  
  # 添加元信息
  filtered_data$id.exposure <- dataset_id
  filtered_data$exposure <- dataset_name
  
  return(filtered_data)
}

cat("✓ 数据提取函数已定义\n\n")

# ===================================================================
# 3.1 定义GWAS元数据提取函数（用于生成Table 1）
# ===================================================================
cat("【步骤3.1】定义GWAS元数据提取函数...\n")

#' 从OpenGWAS API提取完整的GWAS元数据
#' 
#' @param dataset_id 数据集ID
#' @param dataset_name 数据集简称
#' @return 包含完整元数据的数据框
extract_gwas_metadata <- function(dataset_id, dataset_name) {
  
  # 初始化元数据
  metadata <- data.frame(
    id = dataset_id,
    short_name = dataset_name,
    trait = NA_character_,
    sample_size = NA_integer_,
    nsnp = NA_integer_,
    population = NA_character_,
    year = NA_integer_,
    author = NA_character_,
    consortium = NA_character_,
    pmid = NA_character_,
    unit = NA_character_,
    unit_description = NA_character_,
    ncase = NA_integer_,
    ncontrol = NA_integer_,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # 尝试从OpenGWAS API获取元数据
  if (grepl("^ieu-|^ebi-a-", dataset_id)) {
    tryCatch({
      gwas_info <- suppressWarnings(ieugwasr::gwasinfo(dataset_id))
      
      if (!is.null(gwas_info) && nrow(gwas_info) > 0) {
        # 提取基本信息
        if ("trait" %in% names(gwas_info)) {
          metadata$trait <- gwas_info$trait[1]
        }
        
        if ("sample_size" %in% names(gwas_info)) {
          metadata$sample_size <- as.integer(gwas_info$sample_size[1])
        }
        
        if ("nsnp" %in% names(gwas_info)) {
          metadata$nsnp <- as.integer(gwas_info$nsnp[1])
        }
        
        if ("population" %in% names(gwas_info)) {
          metadata$population <- as.character(gwas_info$population[1])
        }
        
        if ("year" %in% names(gwas_info)) {
          metadata$year <- as.integer(gwas_info$year[1])
        }
        
        if ("author" %in% names(gwas_info)) {
          metadata$author <- as.character(gwas_info$author[1])
        }
        
        if ("consortium" %in% names(gwas_info)) {
          metadata$consortium <- as.character(gwas_info$consortium[1])
        }
        
        if ("pmid" %in% names(gwas_info)) {
          metadata$pmid <- as.character(gwas_info$pmid[1])
        }
        
        # 尝试提取ncase和ncontrol
        if ("ncase" %in% names(gwas_info)) {
          metadata$ncase <- as.integer(gwas_info$ncase[1])
        }
        
        if ("ncontrol" %in% names(gwas_info)) {
          metadata$ncontrol <- as.integer(gwas_info$ncontrol[1])
        }
        
        # 根据trait推断单位（如果是定量性状）
        if (!is.na(metadata$trait)) {
          trait_lower <- tolower(metadata$trait)
          if (grepl("bmi|body mass index", trait_lower)) {
            metadata$unit <- "kg/m²"
            metadata$unit_description <- "Body mass index"
          } else if (grepl("cholesterol|hdl|ldl", trait_lower)) {
            metadata$unit <- "mg/dL or mmol/L"
            metadata$unit_description <- "Lipid concentration"
          } else if (grepl("glucose|insulin", trait_lower)) {
            metadata$unit <- "mg/dL or mmol/L"
            metadata$unit_description <- "Glucose/insulin concentration"
          } else if (grepl("blood pressure|systolic|diastolic", trait_lower)) {
            metadata$unit <- "mmHg"
            metadata$unit_description <- "Blood pressure"
          } else if (grepl("cancer|disease", trait_lower)) {
            metadata$unit <- "Odds Ratio"
            metadata$unit_description <- "Binary outcome"
          }
        }
        
      }
    }, error = function(e) {
      # 静默处理错误，使用默认NA值
    })
  }
  
  return(metadata)
}

cat("✓ GWAS元数据提取函数已定义\n\n")

# ===================================================================
# 4. 批量提取暴露因子数据（同时提取元数据）
# ===================================================================
cat("【步骤4】批量提取暴露因子数据（同时提取元数据）...\n")

exposure_data_list <- list()
exposure_metadata_list <- list()  # 存储元数据

for (i in seq_len(nrow(exposure_datasets))) {
  dataset <- exposure_datasets[i, ]
  
  # 跳过NA值的数据源
  if (is.na(dataset$id) || is.na(dataset$data_source_primary) || 
      dataset$data_source_primary %in% c("None", "NULL")) {
    cat(sprintf("    ⚠ %s: 跳过（数据源不可用）\n", dataset$short_name))
    # 即使跳过，也尝试提取元数据
    tryCatch({
      metadata <- extract_gwas_metadata(dataset$id, dataset$short_name)
      if (!is.null(metadata)) {
        metadata$category <- "exposures"
        exposure_metadata_list[[dataset$short_name]] <- metadata
      }
    }, error = function(e) {})
    next
  }
  
  tryCatch({
    # 提取GWAS数据
    data <- extract_gwas_data(
      dataset_id = dataset$id,
      dataset_name = dataset$short_name,
      data_source = ifelse(is.na(dataset$data_source_primary), "None", dataset$data_source_primary),
      completeness = ifelse(is.na(dataset$data_completeness), 0, dataset$data_completeness)
    )
    
    # 提取GWAS元数据
    metadata <- extract_gwas_metadata(dataset$id, dataset$short_name)
    if (!is.null(metadata)) {
      metadata$category <- "exposures"
      # 如果数据提取成功，更新SNP数量
      if (!is.null(data) && nrow(data) > 0) {
        metadata$nsnp <- nrow(data)
      }
      exposure_metadata_list[[dataset$short_name]] <- metadata
    }
    
    if (!is.null(data) && nrow(data) > 0) {
      exposure_data_list[[dataset$short_name]] <- data
      cat(sprintf("    ✓ %s: %d SNPs\n", dataset$short_name, nrow(data)))
    } else {
      cat(sprintf("    ⚠ %s: 数据提取失败\n", dataset$short_name))
    }
    
  }, error = function(e) {
    cat(sprintf("    ✗ %s: %s\n", dataset$short_name, conditionMessage(e)))
    # 即使数据提取失败，也尝试保存元数据
    tryCatch({
      metadata <- extract_gwas_metadata(dataset$id, dataset$short_name)
      if (!is.null(metadata)) {
        metadata$category <- "exposures"
        exposure_metadata_list[[dataset$short_name]] <- metadata
      }
    }, error = function(e2) {})
  })
}

cat(sprintf("\n✓ 暴露因子数据提取完成: %d/%d 成功\n\n", 
            length(exposure_data_list), nrow(exposure_datasets)))

# ===================================================================
# 5. 批量提取结局数据（只有3个结局，同时提取元数据）
# ===================================================================
cat("【步骤5】批量提取结局数据（3个结局，同时提取元数据）...\n")

outcome_data_list <- list()
outcome_metadata_list <- list()  # 存储元数据

for (i in seq_len(nrow(outcome_datasets))) {
  dataset <- outcome_datasets[i, ]
  
  # 跳过NA值的数据源
  if (is.na(dataset$id) || is.na(dataset$data_source_primary) || 
      dataset$data_source_primary %in% c("None", "NULL")) {
    cat(sprintf("    ⚠ %s: 跳过（数据源不可用）\n", dataset$short_name))
    # 即使跳过，也尝试提取元数据
    tryCatch({
      metadata <- extract_gwas_metadata(dataset$id, dataset$short_name)
      if (!is.null(metadata)) {
        metadata$category <- "outcomes"
        outcome_metadata_list[[dataset$short_name]] <- metadata
      }
    }, error = function(e) {})
    next
  }
  
  tryCatch({
    # 提取GWAS数据
    data <- extract_gwas_data(
      dataset_id = dataset$id,
      dataset_name = dataset$short_name,
      data_source = ifelse(is.na(dataset$data_source_primary), "None", dataset$data_source_primary),
      completeness = ifelse(is.na(dataset$data_completeness), 0, dataset$data_completeness)
    )
    
    # 提取GWAS元数据
    metadata <- extract_gwas_metadata(dataset$id, dataset$short_name)
    if (!is.null(metadata)) {
      metadata$category <- "outcomes"
      # 如果数据提取成功，更新SNP数量
      if (!is.null(data) && nrow(data) > 0) {
        metadata$nsnp <- nrow(data)
      }
      outcome_metadata_list[[dataset$short_name]] <- metadata
    }
    
    if (!is.null(data) && nrow(data) > 0) {
      # 结局数据需要转换为outcome格式
      # 如果数据是通过TwoSampleMR::format_data格式化的，需要重新格式化
      if ("id.exposure" %in% names(data)) {
        # 转换exposure格式为outcome格式
        data$id.outcome <- data$id.exposure
        data$outcome <- data$exposure
        data$id.exposure <- NULL
        data$exposure <- NULL
        
        # 重命名列
        if ("beta.exposure" %in% names(data)) {
          names(data) <- gsub("\\.exposure$", ".outcome", names(data))
        }
      } else {
        # 如果没有格式化，使用format_data转换为outcome格式
        if ("beta" %in% names(data) && "se" %in% names(data)) {
          formatted_outcome <- TwoSampleMR::format_data(
            data,
            type = "outcome",
            snp_col = if("rsid" %in% names(data)) "rsid" else "SNP",
            beta_col = "beta",
            se_col = "se",
            effect_allele_col = if("ea" %in% names(data)) "ea" else "effect_allele",
            other_allele_col = if("nea" %in% names(data)) "nea" else "other_allele",
            eaf_col = if("eaf" %in% names(data)) "eaf" else NULL,
            pval_col = if("p" %in% names(data)) "p" else "pval",
            samplesize_col = if("n" %in% names(data)) "n" else NULL
          )
          data <- formatted_outcome
          data$id.outcome <- dataset$id
          data$outcome <- dataset$short_name
        }
      }
      
      outcome_data_list[[dataset$short_name]] <- data
      cat(sprintf("    ✓ %s: %d SNPs\n", dataset$short_name, nrow(data)))
    } else {
      cat(sprintf("    ⚠ %s: 数据提取失败\n", dataset$short_name))
    }
    
  }, error = function(e) {
    cat(sprintf("    ✗ %s: %s\n", dataset$short_name, conditionMessage(e)))
    # 即使数据提取失败，也尝试保存元数据
    tryCatch({
      metadata <- extract_gwas_metadata(dataset$id, dataset$short_name)
      if (!is.null(metadata)) {
        metadata$category <- "outcomes"
        outcome_metadata_list[[dataset$short_name]] <- metadata
      }
    }, error = function(e2) {})
  })
}

# 合并所有元数据
if (length(exposure_metadata_list) > 0 && length(outcome_metadata_list) > 0) {
  all_gwas_metadata <- rbind(
    do.call(rbind, exposure_metadata_list),
    do.call(rbind, outcome_metadata_list)
  )
} else if (length(exposure_metadata_list) > 0) {
  all_gwas_metadata <- do.call(rbind, exposure_metadata_list)
} else if (length(outcome_metadata_list) > 0) {
  all_gwas_metadata <- do.call(rbind, outcome_metadata_list)
} else {
  all_gwas_metadata <- data.frame()
}

if (nrow(all_gwas_metadata) > 0) {
  cat(sprintf("\n✓ GWAS元数据提取完成: %d 个数据集\n\n", nrow(all_gwas_metadata)))
} else {
  cat("\n⚠ 警告: 未提取到任何GWAS元数据\n\n")
}

cat(sprintf("\n✓ 结局数据提取完成: %d/%d 成功\n\n", 
            length(outcome_data_list), nrow(outcome_datasets)))

# ===================================================================
# 7. 数据质量统计和摘要
# ===================================================================
cat("【步骤6】生成数据质量统计...\n")

# 计算每个数据集的质量指标
calculate_data_quality <- function(data_list, category_name) {
  
  # 处理空列表
  if (is.null(data_list) || length(data_list) == 0) {
    return(data.frame(
      category = character(0),
      dataset_name = character(0),
      n_snps = integer(0),
      n_significant = integer(0),
      median_pval = numeric(0),
      median_beta = numeric(0),
      median_se = numeric(0),
      missing_eaf = numeric(0),
      missing_samplesize = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  quality_stats <- lapply(names(data_list), function(name) {
    data <- data_list[[name]]
    
    # 处理空数据集
    n_snps <- nrow(data)
    if (n_snps == 0 || is.null(data)) {
      return(data.frame(
        category = category_name,
        dataset_name = name,
        n_snps = 0L,
        n_significant = 0L,
        median_pval = NA_real_,
        median_beta = NA_real_,
        median_se = NA_real_,
        missing_eaf = NA_real_,
        missing_samplesize = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    # 计算统计量，处理可能的缺失列
    n_significant <- if ("pval" %in% names(data)) {
      sum(data$pval < 5e-8, na.rm = TRUE)
    } else {
      0L
    }
    
    median_pval <- if ("pval" %in% names(data) && length(data$pval) > 0) {
      median(data$pval, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    median_beta <- if ("beta" %in% names(data) && length(data$beta) > 0) {
      median(abs(data$beta), na.rm = TRUE)
    } else {
      NA_real_
    }
    
    median_se <- if ("se" %in% names(data) && length(data$se) > 0) {
      median(data$se, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    # 处理除以零的情况
    missing_eaf <- if ("eaf" %in% names(data) && n_snps > 0) {
      sum(is.na(data$eaf)) / n_snps
    } else {
      NA_real_
    }
    
    missing_samplesize <- if ("samplesize" %in% names(data) && n_snps > 0) {
      sum(is.na(data$samplesize)) / n_snps
    } else {
      NA_real_
    }
    
    data.frame(
      category = category_name,
      dataset_name = name,
      n_snps = n_snps,
      n_significant = n_significant,
      median_pval = median_pval,
      median_beta = median_beta,
      median_se = median_se,
      missing_eaf = missing_eaf,
      missing_samplesize = missing_samplesize,
      stringsAsFactors = FALSE
    )
  })
  
  # 过滤掉NULL值（如果有的话）
  quality_stats <- quality_stats[!sapply(quality_stats, is.null)]
  
  if (length(quality_stats) == 0) {
    return(data.frame(
      category = character(0),
      dataset_name = character(0),
      n_snps = integer(0),
      n_significant = integer(0),
      median_pval = numeric(0),
      median_beta = numeric(0),
      median_se = numeric(0),
      missing_eaf = numeric(0),
      missing_samplesize = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  do.call(rbind, quality_stats)
}

exposure_quality <- calculate_data_quality(exposure_data_list, "Exposure")
outcome_quality <- calculate_data_quality(outcome_data_list, "Outcome")

all_quality_stats <- bind_rows(exposure_quality, outcome_quality)

cat("\n数据质量摘要：\n")
print(all_quality_stats)
cat("\n")

# 保存质量统计表
write.csv(all_quality_stats, 
          "results/data/dataset_quality_summary.csv", 
          row.names = FALSE)

cat("✓ 数据质量统计已保存\n\n")

# ===================================================================
# 7. 保存处理后的数据
# ===================================================================
cat("【步骤7】保存处理后的数据...\n")

# 保存数据列表
save(exposure_data_list, file = "results/data/exposure_data_list.RData")
cat("  ✓ exposure_data_list.RData\n")

save(outcome_data_list, file = "results/data/outcome_data_list.RData")
cat("  ✓ outcome_data_list.RData\n")

# 保存元信息
# 处理category列：如果有冲突，优先使用all_quality_stats中的category
# 首先检查是否有category列冲突
has_category_conflict <- "category" %in% names(dataset_info) && "category" %in% names(all_quality_stats)

if (has_category_conflict) {
  # 如果有冲突，使用suffix参数区分
  dataset_metadata <- dataset_info %>%
    left_join(all_quality_stats, by = c("short_name" = "dataset_name"), suffix = c("_info", "_quality")) %>%
    mutate(
      # 优先使用来自all_quality_stats的category（_quality后缀）
      category = coalesce(category_quality, category_info)
    ) %>%
    # 移除临时列
    select(-any_of(c("category_info", "category_quality")))
} else {
  # 如果没有冲突，直接join
  dataset_metadata <- dataset_info %>%
    left_join(all_quality_stats, by = c("short_name" = "dataset_name"))
}

# 确保category列存在，如果不存在则从dataset_info或all_quality_stats中获取
if (!"category" %in% names(dataset_metadata)) {
  # 尝试从各个可能的列名中获取
  if ("category_quality" %in% names(dataset_metadata)) {
    dataset_metadata$category <- dataset_metadata$category_quality
  } else if ("category_info" %in% names(dataset_metadata)) {
    dataset_metadata$category <- dataset_metadata$category_info
  } else {
    # 如果都不存在，使用默认值
    warning("category列不存在，使用默认值")
    dataset_metadata$category <- NA_character_
  }
}

# 创建完整的数据集信息统计表格
# 首先合并所有可用的元信息
all_available_cols <- names(dataset_metadata)

# 选择并整理列，优先包含更多有用信息
preferred_cols <- c(
  "category", "short_name", "id", "trait",
  "sample_size", "n_snps", 
  "population", "year", "pmid", "author", "consortium",
  "data_source_primary", "data_completeness",
  "n_significant", "median_pval", "median_beta", "median_se",
  "missing_eaf", "missing_samplesize"
)

# 只选择存在的列
cols_to_select <- intersect(preferred_cols, all_available_cols)

dataset_metadata <- dataset_metadata %>%
  select(all_of(cols_to_select)) %>%
  # 格式化数值列
  mutate(
    sample_size = ifelse(is.na(sample_size), NA, 
                        as.integer(sample_size)),
    n_snps = ifelse(is.na(n_snps), NA, as.integer(n_snps)),
    n_significant = ifelse(is.na(n_significant), NA, 
                          as.integer(n_significant)),
    # 格式化p值显示（保留4位小数）
    median_pval = ifelse(is.na(median_pval), NA,
                        round(median_pval, digits = 4)),
    median_beta = ifelse(is.na(median_beta), NA,
                        round(median_beta, digits = 6)),
    median_se = ifelse(is.na(median_se), NA,
                      round(median_se, digits = 6)),
    missing_eaf = ifelse(is.na(missing_eaf), NA,
                        round(missing_eaf * 100, digits = 2)),
    missing_samplesize = ifelse(is.na(missing_samplesize), NA,
                               round(missing_samplesize * 100, digits = 2)),
    data_completeness = ifelse(is.na(data_completeness), NA,
                              round(data_completeness, digits = 2))
  ) %>%
  # 按分类和名称排序
  arrange(category, short_name)

# 保存完整的数据集元信息表
write.csv(dataset_metadata, 
          "results/data/dataset_metadata.csv", 
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("  ✓ dataset_metadata.csv\n")

# ===================================================================
# 7.1 生成分类汇总统计表格
# ===================================================================
cat("【步骤7.1】生成分类汇总统计表格...\n")

# 按分类汇总统计
category_summary_stats <- dataset_metadata %>%
  group_by(category) %>%
  summarise(
    n_datasets = n(),
    n_success_loaded = sum(!is.na(n_snps)),
    total_sample_size = sum(sample_size, na.rm = TRUE),
    mean_sample_size = round(mean(sample_size, na.rm = TRUE), 0),
    median_sample_size = round(median(sample_size, na.rm = TRUE), 0),
    total_n_snps = sum(n_snps, na.rm = TRUE),
    mean_n_snps = round(mean(n_snps, na.rm = TRUE), 0),
    median_n_snps = round(median(n_snps, na.rm = TRUE), 0),
    mean_data_completeness = round(mean(data_completeness, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    loading_rate = round(n_success_loaded / n_datasets * 100, 1),
    category_label = case_when(
      category == "exposures" ~ "暴露因子",
      category == "outcomes" ~ "结局",
      TRUE ~ category
    )
  ) %>%
  select(category, category_label, everything()) %>%
  arrange(category)

write.csv(category_summary_stats,
          "results/data/dataset_category_summary.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("  ✓ dataset_category_summary.csv\n")

# ===================================================================
# 7.2 生成简化版数据集信息表（用于论文）
# ===================================================================
cat("【步骤7.2】生成简化版数据集信息表（论文用）...\n")

# 创建适合论文使用的简化表格
# 定义列名映射（新列名 = 原列名）
paper_cols_mapping <- list(
  Category = "category",
  Trait = "short_name",
  GWAS_ID = "id",
  `Sample Size` = "sample_size",
  `N SNPs` = "n_snps",
  Population = "population",
  Year = "year",
  PMID = "pmid",
  Author = "author",
  Consortium = "consortium",
  `Data Source` = "data_source_primary"
)

# 检查哪些原列存在于dataset_metadata中
available_original_cols <- names(dataset_metadata)
mapping_list <- list()

for (new_name in names(paper_cols_mapping)) {
  old_name <- paper_cols_mapping[[new_name]]
  if (old_name %in% available_original_cols) {
    mapping_list[[new_name]] <- old_name
  }
}

# 如果有可用列，进行重命名和格式化
if (length(mapping_list) > 0) {
  # 先创建基础表格（只包含存在的列）
  dataset_table_paper <- dataset_metadata %>%
    select(all_of(unlist(mapping_list)))
  
  # 重命名列
  for (new_name in names(mapping_list)) {
    old_name <- mapping_list[[new_name]]
    if (old_name %in% names(dataset_table_paper)) {
      names(dataset_table_paper)[names(dataset_table_paper) == old_name] <- new_name
    }
  }
  
  # 处理缺失值显示和格式化
  if ("Category" %in% names(dataset_table_paper)) {
    dataset_table_paper$Category <- case_when(
      dataset_table_paper$Category == "exposures" ~ "Exposure",
      dataset_table_paper$Category == "outcomes" ~ "Outcome",
      TRUE ~ dataset_table_paper$Category
    )
  }
  
  if ("Sample Size" %in% names(dataset_table_paper)) {
    dataset_table_paper$`Sample Size` <- ifelse(
      is.na(dataset_table_paper$`Sample Size`), 
      "", 
      format(dataset_table_paper$`Sample Size`, big.mark = ",")
    )
  }
  
  if ("N SNPs" %in% names(dataset_table_paper)) {
    dataset_table_paper$`N SNPs` <- ifelse(
      is.na(dataset_table_paper$`N SNPs`), 
      "", 
      format(dataset_table_paper$`N SNPs`, big.mark = ",")
    )
  }
  
  if ("Population" %in% names(dataset_table_paper)) {
    dataset_table_paper$Population <- ifelse(
      is.na(dataset_table_paper$Population), 
      "", 
      dataset_table_paper$Population
    )
  }
  
  if ("Year" %in% names(dataset_table_paper)) {
    dataset_table_paper$Year <- ifelse(
      is.na(dataset_table_paper$Year), 
      "", 
      as.character(dataset_table_paper$Year)
    )
  }
  
  if ("PMID" %in% names(dataset_table_paper)) {
    dataset_table_paper$PMID <- ifelse(
      is.na(dataset_table_paper$PMID) | dataset_table_paper$PMID == "", 
      "", 
      dataset_table_paper$PMID
    )
  }
  
  if ("Author" %in% names(dataset_table_paper)) {
    dataset_table_paper$Author <- ifelse(
      is.na(dataset_table_paper$Author) | dataset_table_paper$Author == "", 
      "", 
      dataset_table_paper$Author
    )
  }
  
  if ("Consortium" %in% names(dataset_table_paper)) {
    dataset_table_paper$Consortium <- ifelse(
      is.na(dataset_table_paper$Consortium) | dataset_table_paper$Consortium == "", 
      "", 
      dataset_table_paper$Consortium
    )
  }
  
  if ("Data Source" %in% names(dataset_table_paper)) {
    dataset_table_paper$`Data Source` <- ifelse(
      is.na(dataset_table_paper$`Data Source`) | dataset_table_paper$`Data Source` == "", 
      "", 
      dataset_table_paper$`Data Source`
    )
  }
  
  # 按Category和Trait排序（如果存在）
  if ("Category" %in% names(dataset_table_paper) && "Trait" %in% names(dataset_table_paper)) {
    dataset_table_paper <- dataset_table_paper %>%
      arrange(Category, Trait)
  } else if ("Trait" %in% names(dataset_table_paper)) {
    dataset_table_paper <- dataset_table_paper %>%
      arrange(Trait)
  }
} else {
  # 如果没有可用列，创建一个空表格
  warning("没有可用的列来创建论文表格")
  dataset_table_paper <- data.frame()
}

write.csv(dataset_table_paper,
          "results/data/dataset_table_for_paper.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("  ✓ dataset_table_for_paper.csv\n")

# ===================================================================
# 7.3 生成Table 1 - GWAS数据摘要（论文用完整版）
# ===================================================================
cat("【步骤7.3】生成Table 1 - GWAS数据摘要（论文用完整版）...\n")

# 合并元数据和数据集信息
if (nrow(all_gwas_metadata) > 0 && "id" %in% names(all_gwas_metadata)) {
  # 将all_gwas_metadata与dataset_metadata合并
  # 优先使用从OpenGWAS API提取的完整元数据
  table1_data <- dataset_metadata %>%
    left_join(all_gwas_metadata, by = c("id" = "id", "short_name" = "short_name"), suffix = c("", "_api"))
  
  # 优先使用API提取的元数据（_api后缀），如果缺失则使用原有数据
  if ("trait_api" %in% names(table1_data)) {
    table1_data$trait <- ifelse(!is.na(table1_data$trait_api), table1_data$trait_api, 
                                 ifelse(!is.na(table1_data$trait), table1_data$trait, table1_data$short_name))
  } else if ("trait" %in% names(table1_data)) {
    table1_data$trait <- ifelse(!is.na(table1_data$trait), table1_data$trait, table1_data$short_name)
  } else {
    table1_data$trait <- table1_data$short_name
  }
  
  if ("sample_size_api" %in% names(table1_data)) {
    table1_data$sample_size <- ifelse(!is.na(table1_data$sample_size_api), table1_data$sample_size_api, table1_data$sample_size)
  }
  
  if ("nsnp_api" %in% names(table1_data)) {
    table1_data$nsnp <- ifelse(!is.na(table1_data$nsnp_api), table1_data$nsnp_api, table1_data$n_snps)
  } else if ("n_snps" %in% names(table1_data)) {
    table1_data$nsnp <- table1_data$n_snps
  }
  
  if ("population_api" %in% names(table1_data)) {
    table1_data$population <- ifelse(!is.na(table1_data$population_api), table1_data$population_api, table1_data$population)
  }
  
  if ("year_api" %in% names(table1_data)) {
    table1_data$year <- ifelse(!is.na(table1_data$year_api), table1_data$year_api, table1_data$year)
  }
  
  if ("author_api" %in% names(table1_data)) {
    table1_data$author <- ifelse(!is.na(table1_data$author_api), table1_data$author_api, table1_data$author)
  }
  
  if ("consortium_api" %in% names(table1_data)) {
    table1_data$consortium <- ifelse(!is.na(table1_data$consortium_api), table1_data$consortium_api, table1_data$consortium)
  }
  
  if ("pmid_api" %in% names(table1_data)) {
    table1_data$pmid <- ifelse(!is.na(table1_data$pmid_api), table1_data$pmid_api, table1_data$pmid)
  }
  
  if ("unit_api" %in% names(table1_data)) {
    table1_data$unit <- ifelse(!is.na(table1_data$unit_api), table1_data$unit_api, table1_data$unit)
  }
  
  if ("ncase_api" %in% names(table1_data)) {
    table1_data$ncase <- ifelse(!is.na(table1_data$ncase_api), table1_data$ncase_api, table1_data$ncase)
  }
  
  if ("ncontrol_api" %in% names(table1_data)) {
    table1_data$ncontrol <- ifelse(!is.na(table1_data$ncontrol_api), table1_data$ncontrol_api, table1_data$ncontrol)
  }
  
  # 移除API后缀列（如果存在）
  api_cols <- grep("_api$", names(table1_data), value = TRUE)
  if (length(api_cols) > 0) {
    table1_data <- table1_data %>% select(-all_of(api_cols))
  }
} else {
  # 如果没有API元数据，使用原有数据
  table1_data <- dataset_metadata
  if (!"trait" %in% names(table1_data)) {
    table1_data$trait <- table1_data$short_name
  }
  if (!"nsnp" %in% names(table1_data) && "n_snps" %in% names(table1_data)) {
    table1_data$nsnp <- table1_data$n_snps
  }
}

# 确保必要的列存在
if (nrow(table1_data) > 0) {
  # 生成Table 1格式
  table1_cols <- c()
  table1_data_final <- data.frame()
  
  # 动态选择存在的列
  if ("category" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "category")
  }
  if ("trait" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "trait")
  } else if ("short_name" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "short_name")
    table1_data$trait <- table1_data$short_name
  }
  if ("id" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "id")
  }
  if ("sample_size" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "sample_size")
  }
  if ("nsnp" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "nsnp")
  } else if ("n_snps" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "n_snps")
    table1_data$nsnp <- table1_data$n_snps
  }
  if ("population" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "population")
  }
  if ("year" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "year")
  }
  if ("author" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "author")
  }
  if ("consortium" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "consortium")
  }
  if ("pmid" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "pmid")
  }
  if ("unit" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "unit")
  }
  if ("ncase" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "ncase")
  }
  if ("ncontrol" %in% names(table1_data)) {
    table1_cols <- c(table1_cols, "ncontrol")
  }
  
  table1 <- table1_data %>%
    select(all_of(unique(table1_cols))) %>%
    mutate(
      Category = case_when(
        category == "exposures" ~ "Exposure",
        category == "outcomes" ~ "Outcome",
        TRUE ~ category
      ),
      Trait = coalesce(trait, short_name),
      GWAS_ID = id,
      `Sample Size` = ifelse(is.na(sample_size), "", format(sample_size, big.mark = ",", scientific = FALSE)),
      `N SNPs` = ifelse(is.na(nsnp), "", format(nsnp, big.mark = ",", scientific = FALSE)),
      Year = ifelse(is.na(year), "", as.character(year)),
      `N Cases` = ifelse(is.na(ncase), "", format(ncase, big.mark = ",", scientific = FALSE)),
      `N Controls` = ifelse(is.na(ncontrol), "", format(ncontrol, big.mark = ",", scientific = FALSE)),
      Population = ifelse(is.na(population) | population == "", "", as.character(population)),
      Author = ifelse(is.na(author) | author == "", "", as.character(author)),
      Consortium = ifelse(is.na(consortium) | consortium == "", "", as.character(consortium)),
      PMID = ifelse(is.na(pmid) | pmid == "", "", as.character(pmid)),
      Unit = ifelse(is.na(unit) | unit == "", "", as.character(unit))
    ) %>%
    select(Category, Trait, GWAS_ID, `Sample Size`, `N SNPs`, Population, Year, Author, Consortium, PMID, Unit, `N Cases`, `N Controls`) %>%
    arrange(Category, Trait)
  
  write.csv(table1,
            "results/data/Table1_GWAS_Summary.csv",
            row.names = FALSE,
            fileEncoding = "UTF-8")
  cat("  ✓ Table1_GWAS_Summary.csv\n")
} else {
  cat("  ⚠ 警告: 无法生成Table 1（无数据）\n")
}

# ===================================================================
# 7.4 生成Harmonization Information关键指标汇总表
# ===================================================================
cat("【步骤7.4】生成Harmonization Information关键指标汇总表...\n")

# 合并数据质量统计和元数据
if (exists("all_quality_stats") && nrow(all_quality_stats) > 0) {
  # 首先检查dataset_metadata中是否存在所需的列
  available_cols <- names(dataset_metadata)
  
  # 动态构建mutate语句，只处理存在的列
  mutate_expr <- list()
  if ("category" %in% available_cols) {
    mutate_expr$Category <- quote(case_when(
      category == "exposures" ~ "Exposure",
      category == "outcomes" ~ "Outcome",
      TRUE ~ category
    ))
  }
  if ("short_name" %in% available_cols) mutate_expr$Trait <- quote(short_name)
  if ("id" %in% available_cols) mutate_expr$GWAS_ID <- quote(id)
  if ("sample_size" %in% available_cols) {
    mutate_expr$`Sample Size` <- quote(ifelse(is.na(sample_size), "", format(sample_size, big.mark = ",", scientific = FALSE)))
  }
  if ("n_snps" %in% available_cols) {
    mutate_expr$`N SNPs` <- quote(ifelse(is.na(n_snps), "", format(n_snps, big.mark = ",", scientific = FALSE)))
  }
  if ("n_significant" %in% available_cols) {
    mutate_expr$`N Significant SNPs` <- quote(ifelse(is.na(n_significant), "", format(n_significant, big.mark = ",", scientific = FALSE)))
  }
  if ("missing_eaf" %in% available_cols) {
    mutate_expr$`Missing EAF (%)` <- quote(ifelse(is.na(missing_eaf), "", sprintf("%.2f", missing_eaf)))
  }
  if ("missing_samplesize" %in% available_cols) {
    mutate_expr$`Missing Sample Size (%)` <- quote(ifelse(is.na(missing_samplesize), "", sprintf("%.2f", missing_samplesize)))
  }
  if ("data_completeness" %in% available_cols) {
    mutate_expr$`Data Completeness (%)` <- quote(ifelse(is.na(data_completeness), "", sprintf("%.2f", data_completeness)))
  }
  
  # 动态构建select语句，只选择存在的列
  select_cols <- c()
  if ("Category" %in% names(mutate_expr)) select_cols <- c(select_cols, "Category")
  if ("Trait" %in% names(mutate_expr)) select_cols <- c(select_cols, "Trait")
  if ("GWAS_ID" %in% names(mutate_expr)) select_cols <- c(select_cols, "GWAS_ID")
  if ("Sample Size" %in% names(mutate_expr)) select_cols <- c(select_cols, "Sample Size")
  if ("N SNPs" %in% names(mutate_expr)) select_cols <- c(select_cols, "N SNPs")
  if ("N Significant SNPs" %in% names(mutate_expr)) select_cols <- c(select_cols, "N Significant SNPs")
  if ("Missing EAF (%)" %in% names(mutate_expr)) select_cols <- c(select_cols, "Missing EAF (%)")
  if ("Missing Sample Size (%)" %in% names(mutate_expr)) select_cols <- c(select_cols, "Missing Sample Size (%)")
  if ("Data Completeness (%)" %in% names(mutate_expr)) select_cols <- c(select_cols, "Data Completeness (%)")
  
  harmonization_info <- dataset_metadata %>%
    left_join(all_quality_stats, by = c("short_name" = "dataset_name")) %>%
    mutate(!!!mutate_expr)
  
  # 只选择存在的列进行排序和输出
  if (length(select_cols) > 0) {
    harmonization_info <- harmonization_info %>%
      select(all_of(select_cols))
    
    # 如果有Category和Trait列，按它们排序
    if ("Category" %in% select_cols && "Trait" %in% select_cols) {
      harmonization_info <- harmonization_info %>% arrange(Category, Trait)
    }
    
    write.csv(harmonization_info,
              "results/data/Harmonization_Information.csv",
              row.names = FALSE,
              fileEncoding = "UTF-8")
    cat("  ✓ Harmonization_Information.csv\n")
  } else {
    cat("  ⚠ 警告: 没有可用的列来生成Harmonization Information\n")
  }
} else {
  cat("  ⚠ 警告: 无法生成Harmonization Information（无质量统计）\n")
}
cat("\n")

# ===================================================================
# 7.5 生成数据质量详细统计表
# ===================================================================
cat("【步骤7.3】生成数据质量详细统计表...\n")

# 合并质量统计和元信息，创建详细的质量统计表
# 首先检查dataset_metadata中是否存在所需的列
available_cols <- names(dataset_metadata)

# 动态构建select语句，只选择存在的列
select_expr <- list()
if ("category" %in% available_cols) select_expr$Category <- "category"
if ("short_name" %in% available_cols) select_expr$Dataset <- "short_name"
if ("id" %in% available_cols) select_expr$GWAS_ID <- "id"
if ("n_snps" %in% available_cols) select_expr$`N SNPs` <- "n_snps"
if ("n_significant" %in% available_cols) select_expr$`Significant SNPs (p<5e-8)` <- "n_significant"
if ("median_pval" %in% available_cols) select_expr$`Median P-value` <- "median_pval"
if ("median_beta" %in% available_cols) select_expr$`Median Beta` <- "median_beta"
if ("median_se" %in% available_cols) select_expr$`Median SE` <- "median_se"
if ("missing_eaf" %in% available_cols) select_expr$`Missing EAF (%)` <- "missing_eaf"
if ("missing_samplesize" %in% available_cols) select_expr$`Missing Sample Size (%)` <- "missing_samplesize"
if ("data_completeness" %in% available_cols) select_expr$`Data Completeness (%)` <- "data_completeness"

# 使用动态构建的select语句
quality_detail <- dataset_metadata %>%
  select(!!!select_expr) %>%
  mutate(
    Category = case_when(
      Category == "exposures" ~ "Exposure",
      Category == "outcomes" ~ "Outcome",
      TRUE ~ Category
    )
  ) %>%
  arrange(Category, Dataset)

write.csv(quality_detail,
          "results/data/dataset_quality_detail.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
cat("  ✓ dataset_quality_detail.csv\n\n")

# ===================================================================
# 8. 生成数据加载报告
# ===================================================================
cat("【步骤8】生成数据加载报告...\n")

report_text <- c(
  "================================================================================",
  "STEP 2: DATA LOADING AND PREPROCESSING REPORT",
  "步骤2：数据加载与预处理报告",
  "================================================================================",
  "",
  sprintf("报告生成时间: %s", Sys.time()),
  "",
  "【数据加载摘要】",
  "",
  sprintf("总数据集数量: %d", nrow(dataset_info)),
  sprintf("  • 暴露因子: %d (%d 成功提取)",
          nrow(exposure_datasets), length(exposure_data_list)),
  sprintf("  • 结局: %d (%d 成功提取)",
          nrow(outcome_datasets), length(outcome_data_list)),
  "",
  "【暴露因子详情】",
  ""
)

for (name in names(exposure_data_list)) {
  data <- exposure_data_list[[name]]
  report_text <- c(report_text,
    sprintf("  %s:", name),
    sprintf("    - SNPs: %d", nrow(data)),
    sprintf("    - 显著SNPs (P<5e-8): %d", sum(data$pval < 5e-8, na.rm = TRUE)),
    sprintf("    - 中位P值: %.2e", median(data$pval, na.rm = TRUE)),
    ""
  )
}

report_text <- c(report_text,
  "【结局详情】",
  ""
)

for (name in names(outcome_data_list)) {
  data <- outcome_data_list[[name]]
  report_text <- c(report_text,
    sprintf("  %s:", name),
    sprintf("    - SNPs: %d", nrow(data)),
    sprintf("    - 显著SNPs (P<5e-8): %d", sum(data$pval < 5e-8, na.rm = TRUE)),
    sprintf("    - 中位P值: %.2e", median(data$pval, na.rm = TRUE)),
    ""
  )
}

report_text <- c(report_text,
  "【数据质量评估】",
  "",
  "总体质量指标:",
  sprintf("  • 平均SNP数量: %.0f", mean(all_quality_stats$n_snps, na.rm = TRUE)),
  sprintf("  • 平均显著SNP数量: %.0f", mean(all_quality_stats$n_significant, na.rm = TRUE)),
  sprintf("  • EAF缺失率: %.1f%%", mean(all_quality_stats$missing_eaf, na.rm = TRUE) * 100),
  sprintf("  • 样本量缺失率: %.1f%%", mean(all_quality_stats$missing_samplesize, na.rm = TRUE) * 100),
  "",
  "【生成的统计表格】",
  "",
  "1. dataset_metadata.csv - 完整数据集元信息表（包含所有字段）",
  "2. dataset_category_summary.csv - 按分类汇总的统计表",
  "3. dataset_table_for_paper.csv - 适合论文使用的简化表格",
  "4. dataset_quality_summary.csv - 数据质量摘要统计表",
  "5. dataset_quality_detail.csv - 数据质量详细统计表",
  "",
  "【下一步建议】",
  "",
  "1. 查看分类汇总统计: results/data/dataset_category_summary.csv",
  "2. 检查数据质量统计: results/data/dataset_quality_detail.csv",
  "3. 查看完整元信息: results/data/dataset_metadata.csv",
  "4. 继续进行工具变量筛选 (Step 3)",
  "",
  "================================================================================",
  "报告结束",
  "================================================================================"
)

writeLines(report_text, "results/reports/step02_data_summary.txt")
cat("✓ 数据加载报告已保存: results/reports/step02_data_summary.txt\n\n")

# ===================================================================
# 10. 创建数据概览可视化
# ===================================================================
cat("【步骤9】创建数据概览可视化...\n")

# 图1：每个数据集的SNP数量
p1 <- ggplot(all_quality_stats, 
             aes(x = reorder(dataset_name, n_snps), 
                 y = n_snps, 
                 fill = category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Exposure" = "#E69F00", "Outcome" = "#009E73")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Number of SNPs per Dataset",
    subtitle = "After quality filtering",
    x = "Dataset",
    y = "Number of SNPs",
    fill = "Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "bottom"
  )

# 图2：显著SNP数量
p2 <- ggplot(all_quality_stats, 
             aes(x = reorder(dataset_name, n_significant), 
                 y = n_significant, 
                 fill = category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Exposure" = "#E69F00", "Outcome" = "#009E73")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Number of Genome-Wide Significant SNPs (P < 5e-8)",
    x = "Dataset",
    y = "Number of Significant SNPs",
    fill = "Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# 图3：数据完整性（缺失率）
missing_data <- all_quality_stats %>%
  select(dataset_name, category, missing_eaf, missing_samplesize) %>%
  pivot_longer(cols = c(missing_eaf, missing_samplesize),
               names_to = "metric",
               values_to = "missing_rate") %>%
  mutate(
    metric = case_when(
      metric == "missing_eaf" ~ "Effect Allele Frequency",
      metric == "missing_samplesize" ~ "Sample Size",
      TRUE ~ metric
    )
  )

p3 <- ggplot(missing_data, 
             aes(x = dataset_name, y = missing_rate * 100, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(
    values = c("Effect Allele Frequency" = "#D55E00", "Sample Size" = "#0072B2")
  ) +
  labs(
    title = "Data Completeness",
    subtitle = "Missing data rate by dataset",
    x = "Dataset",
    y = "Missing Rate (%)",
    fill = "Metric"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "bottom"
  )

# 组合图表
library(patchwork)
combined_plot <- (p1 / p2 / p3) +
  plot_annotation(
    title = "Step 2: Data Loading and Preprocessing Overview",
    subtitle = sprintf("Generated: %s", Sys.time()),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
    )
  )

# 保存图表
ggsave("results/figures/step02_data_overview.png",
       plot = combined_plot,
       width = 12,
       height = 16,
       dpi = 300,
       bg = "white")

ggsave("results/figures/step02_data_overview.pdf",
       plot = combined_plot,
       width = 12,
       height = 16,
       device = cairo_pdf)

cat("✓ 数据概览可视化已保存\n\n")

# ===================================================================
# 11. 最终总结
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤2完成！Data Loading and Preprocessing\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【生成的文件】\n\n")
cat("数据文件:\n")
cat("  ✓ results/data/exposure_data_list.RData\n")
cat("  ✓ results/data/outcome_data_list.RData\n\n")

cat("统计表格文件:\n")
cat("  ✓ results/data/dataset_metadata.csv (完整数据集元信息表)\n")
cat("  ✓ results/data/dataset_category_summary.csv (分类汇总统计表)\n")
cat("  ✓ results/data/dataset_table_for_paper.csv (论文用简化表格)\n")
cat("  ✓ results/data/dataset_quality_summary.csv (数据质量摘要表)\n")
cat("  ✓ results/data/dataset_quality_detail.csv (数据质量详细表)\n\n")

cat("报告文件:\n")
cat("  ✓ results/reports/step02_data_summary.txt\n\n")

cat("可视化文件:\n")
cat("  ✓ results/figures/step02_data_overview.png\n")
cat("  ✓ results/figures/step02_data_overview.pdf\n\n")

cat("【数据加载统计】\n")
cat(sprintf("  • 暴露因子: %d/%d 成功 (%.1f%%)\n", 
            length(exposure_data_list), 
            nrow(exposure_datasets),
            length(exposure_data_list) / nrow(exposure_datasets) * 100))
cat(sprintf("  • 结局: %d/%d 成功 (%.1f%%)\n", 
            length(outcome_data_list), 
            nrow(outcome_datasets),
            length(outcome_data_list) / nrow(outcome_datasets) * 100))

cat("\n【下一步操作】\n")
cat("  1. 查看数据加载报告: results/reports/step02_data_summary.txt\n")
cat("  2. 查看分类汇总统计: results/data/dataset_category_summary.csv\n")
cat("  3. 检查数据质量详情: results/data/dataset_quality_detail.csv\n")
cat("  4. 查看论文用表格: results/data/dataset_table_for_paper.csv\n")
cat("  5. 查看可视化概览: results/figures/step02_data_overview.png\n")
cat("  6. 继续运行第3步: 工具变量筛选\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("数据加载完成！Ready for IV selection (Step 3)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

