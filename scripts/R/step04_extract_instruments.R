############################################################################
# 步骤4：提取工具变量（从Step2数据加载）
# 孟德尔随机化研究 - 代谢性状、炎症标志物与肺癌亚型
# 基于Step2生成的数据文件（exposure_data_list.RData）
############################################################################

cat("步骤4：开始提取工具变量（从Step2数据）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ===================================================================
# 0. 加载必要的包
# ===================================================================
cat("【步骤0】加载必要的包...\n")

required_packages <- c(
  "TwoSampleMR",    # MR分析核心包
  "ieugwasr",       # IEU OpenGWAS数据库接口（用于clump）
  "dplyr",          # 数据处理
  "ggplot2"         # 数据可视化
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  安装缺失的包: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# 确保包已正确加载
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("%s包未安装，请先安装该包", pkg))
  }
}

cat("✓ 所有必要的包已加载\n\n")

# ===================================================================
# 1. 加载Step2生成的数据
# ===================================================================
cat("【步骤1】加载Step2生成的数据...\n")

# 加载暴露数据列表
exposure_data_file <- "results/data/exposure_data_list.RData"
if (!file.exists(exposure_data_file)) {
  stop(sprintf("错误：找不到Step2生成的暴露数据文件: %s\n请先运行Step2脚本", exposure_data_file))
}

load(exposure_data_file)
cat(sprintf("✓ 已加载暴露数据列表：%d 个暴露因子\n", length(exposure_data_list)))

# 加载结局数据列表（用于后续分析）
outcome_data_file <- "results/data/outcome_data_list.RData"
if (file.exists(outcome_data_file)) {
  load(outcome_data_file)
  cat(sprintf("✓ 已加载结局数据列表：%d 个结局\n", length(outcome_data_list)))
} else {
  warning(sprintf("未找到结局数据文件: %s（可选）", outcome_data_file))
  outcome_data_list <- list()
}

# 加载数据集元信息（可选）
metadata_file <- "results/data/dataset_metadata.csv"
if (file.exists(metadata_file)) {
  dataset_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ 已加载数据集元信息：%d 个数据集\n", nrow(dataset_metadata)))
} else {
  warning(sprintf("未找到数据集元信息文件: %s（可选）", metadata_file))
  dataset_metadata <- NULL
}

# 创建输出目录
output_dirs <- c(
  "data",
  "results",
  "results/instruments",
  "results/figures",
  "results/reports"
)

for (dir in output_dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

cat("\n")

# ===================================================================
# 2. 定义工具变量提取和筛选函数
# ===================================================================
cat("【步骤2】定义工具变量提取函数...\n")

extract_instruments_from_data <- function(
  data, 
  trait_name, 
                                        p_threshold = 5e-8,
                                        clump = TRUE,
                                        r2_threshold = 0.001,
  kb_threshold = 10000
) {
  
  tryCatch({
    # 输入验证
    if (is.null(data) || nrow(data) == 0) {
      warning(sprintf("提取工具变量：%s的数据为空", trait_name))
      return(NULL)
    }
    
    cat(sprintf("  → 处理: %s (初始SNP数: %d)\n", trait_name, nrow(data)))
    
    # 检查数据格式，可能需要格式化
    data_formatted <- NULL
    
    # 如果数据已经格式化为TwoSampleMR格式（有exposure相关列）
    if ("id.exposure" %in% names(data) || "exposure" %in% names(data)) {
      data_formatted <- data
      cat("    ✓ 数据已格式化\n")
    } else {
      # 尝试识别列名并格式化
      # 常见的列名映射
      snp_col <- NULL
      pval_col <- NULL
      beta_col <- NULL
      se_col <- NULL
      eaf_col <- NULL
      effect_allele_col <- NULL
      other_allele_col <- NULL
      
      # 查找SNP列
      for (col in c("SNP", "rsid", "rs_id", "variant", "markername")) {
        if (col %in% names(data)) {
          snp_col <- col
          break
        }
      }
      
      # 查找P值列
      for (col in c("pval", "p", "p_value", "P", "pval.exposure")) {
        if (col %in% names(data)) {
          pval_col <- col
          break
        }
      }
      
      # 查找beta列
      for (col in c("beta", "b", "BETA", "beta.exposure")) {
        if (col %in% names(data)) {
          beta_col <- col
          break
        }
      }
      
      # 查找se列
      for (col in c("se", "SE", "std_err", "standard_error", "se.exposure")) {
        if (col %in% names(data)) {
          se_col <- col
          break
        }
      }
      
      # 查找eaf列
      for (col in c("eaf", "EAF", "ea_freq", "effect_allele_frequency", "eaf.exposure")) {
        if (col %in% names(data)) {
          eaf_col <- col
          break
        }
      }
      
      # 查找effect_allele列
      for (col in c("ea", "effect_allele", "A1", "effect_allele.exposure")) {
        if (col %in% names(data)) {
          effect_allele_col <- col
          break
        }
      }
      
      # 查找other_allele列
      for (col in c("nea", "other_allele", "A2", "other_allele.exposure")) {
        if (col %in% names(data)) {
          other_allele_col <- col
          break
        }
      }
      
      # 如果找到了必要的列，进行格式化
      if (!is.null(snp_col) && !is.null(pval_col) && !is.null(beta_col) && !is.null(se_col)) {
        tryCatch({
          data_formatted <- TwoSampleMR::format_data(
            data,
            type = "exposure",
            snp_col = snp_col,
            beta_col = beta_col,
            se_col = se_col,
            effect_allele_col = effect_allele_col,
            other_allele_col = other_allele_col,
            eaf_col = eaf_col,
            pval_col = pval_col
          )
          cat("    ✓ 数据格式化成功\n")
        }, error = function(e) {
          warning(sprintf("数据格式化失败: %s", conditionMessage(e)))
          return(NULL)
        })
      } else {
        warning(sprintf("无法识别数据列名，缺少必要列"))
        return(NULL)
      }
    }
    
    if (is.null(data_formatted)) {
      return(NULL)
    }
    
    # 步骤1：筛选显著性SNPs (p < threshold)
    if ("pval.exposure" %in% names(data_formatted)) {
      data_filtered <- data_formatted[data_formatted$pval.exposure < p_threshold, ]
    } else if ("pval" %in% names(data_formatted)) {
      data_filtered <- data_formatted[data_formatted$pval < p_threshold, ]
    } else {
      warning(sprintf("未找到P值列，跳过筛选"))
      data_filtered <- data_formatted
    }
    
    cat(sprintf("    ✓ 筛选后: %d SNPs (p < %g)\n", nrow(data_filtered), p_threshold))
    
    if (nrow(data_filtered) == 0) {
      warning(sprintf("筛选后无SNPs满足p < %g", p_threshold))
      return(NULL)
    }
    
    # 步骤2：LD clumping（如果需要）
    instruments <- data_filtered
    if (clump && nrow(data_filtered) > 0) {
      cat("    → 进行LD clumping...\n")
      
      # 检查是否需要添加rsid列
      if (!"rsid" %in% names(instruments)) {
        if ("SNP" %in% names(instruments)) {
          instruments$rsid <- instruments$SNP
        } else if ("id" %in% names(instruments)) {
          instruments$rsid <- instruments$id
        } else {
          warning("无法找到rsid列，跳过clumping")
          clump <- FALSE
        }
      }
      
      if (clump && "rsid" %in% names(instruments)) {
        tryCatch({
          # 准备clump数据
          pval_col_name <- if("pval.exposure" %in% names(instruments)) "pval.exposure" else "pval"
          pval_values <- if(pval_col_name %in% names(instruments)) {
            instruments[[pval_col_name]]
          } else if("pval" %in% names(instruments)) {
            instruments$pval
          } else {
            rep(NA, nrow(instruments))
          }
          
          clump_data <- data.frame(
            rsid = instruments$rsid,
            pval = pval_values,
            stringsAsFactors = FALSE
          )
          clump_data <- clump_data[!is.na(clump_data$rsid) & !is.na(clump_data$pval), ]
          
          if (nrow(clump_data) > 0) {
            # 进行clumping
            clumped <- ieugwasr::ld_clump(
              dat = clump_data,
              clump_r2 = r2_threshold,
              clump_kb = kb_threshold,
              pop = "EUR"  # 默认使用欧洲人群
            )
            
            # 保留clumped后的SNPs
            instruments <- instruments[instruments$rsid %in% clumped$rsid, ]
            cat(sprintf("    ✓ Clumping后: %d SNPs (r2 < %.3f, kb > %d)\n", 
                        nrow(instruments), r2_threshold, kb_threshold))
          } else {
            warning("Clump数据准备失败")
          }
        }, error = function(e) {
          warning(sprintf("Clumping失败: %s，使用未clump的数据", conditionMessage(e)))
        })
      }
    }
    
    if (nrow(instruments) == 0) {
      warning(sprintf("Clumping后无SNPs保留"))
      return(NULL)
    }
    
    # 步骤3：计算质量指标
    # 计算F统计量
    if ("beta.exposure" %in% names(instruments) && "se.exposure" %in% names(instruments)) {
    instruments$F_statistic <- (instruments$beta.exposure^2) / (instruments$se.exposure^2)
    } else if ("beta" %in% names(instruments) && "se" %in% names(instruments)) {
      instruments$F_statistic <- (instruments$beta^2) / (instruments$se^2)
    } else {
      instruments$F_statistic <- NA_real_
      warning("无法计算F统计量（缺少beta或se列）")
    }
    
    # 计算R²
    eaf_col <- if("eaf.exposure" %in% names(instruments)) "eaf.exposure" else "eaf"
    beta_col <- if("beta.exposure" %in% names(instruments)) "beta.exposure" else "beta"
    
    if (eaf_col %in% names(instruments) && beta_col %in% names(instruments)) {
      valid_eaf <- !is.na(instruments[[eaf_col]]) & 
                   instruments[[eaf_col]] >= 0 & 
                   instruments[[eaf_col]] <= 1
      instruments$r_squared <- NA_real_
      if (sum(valid_eaf) > 0) {
        instruments$r_squared[valid_eaf] <- 2 * instruments[[eaf_col]][valid_eaf] * 
                                          (1 - instruments[[eaf_col]][valid_eaf]) * 
                                          (instruments[[beta_col]][valid_eaf]^2)
      }
    } else {
      instruments$r_squared <- NA_real_
    }
    
    # 标记弱工具变量 (F < 10)
    instruments$weak_instrument <- ifelse(
      is.na(instruments$F_statistic), 
      NA, 
      instruments$F_statistic < 10
    )
    
    # 计算总解释方差
    total_r2 <- sum(instruments$r_squared, na.rm = TRUE)
    instruments$variance_explained <- total_r2
    
    # 报告结果
    n_weak <- sum(instruments$weak_instrument, na.rm = TRUE)
    mean_f <- mean(instruments$F_statistic, na.rm = TRUE)
    
    cat(sprintf("    ✓ 成功提取 %d 个工具变量\n", nrow(instruments)))
    if (!is.na(mean_f)) {
      cat(sprintf("    ✓ 平均F统计量: %.2f\n", mean_f))
    }
    if (nrow(instruments) > 0) {
      cat(sprintf("    ✓ 弱工具变量: %d (%.1f%%)\n", 
                  n_weak, round(n_weak / nrow(instruments) * 100, 1)))
    }
    
    return(instruments)
    
  }, error = function(e) {
    cat(sprintf("    ✗ 错误: %s\n", conditionMessage(e)))
    return(NULL)
  })
}

cat("✓ 函数定义完成\n\n")

# ===================================================================
# 3. 批量提取工具变量
# ===================================================================
cat("【步骤3】批量提取暴露因子的工具变量...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 初始化结果存储
all_instruments <- list()
instruments_summary <- data.frame()
extraction_log <- list()

total_traits <- length(exposure_data_list)
current_trait <- 0

for (trait_name in names(exposure_data_list)) {
  current_trait <- current_trait + 1
  cat(sprintf("\n[%d/%d] 处理暴露因子: %s\n", current_trait, total_traits, trait_name))
  
  data <- exposure_data_list[[trait_name]]
  
  start_time <- Sys.time()
  
  instruments <- extract_instruments_from_data(
    data = data,
    trait_name = trait_name,
    p_threshold = 5e-8,
    clump = TRUE,
    r2_threshold = 0.001,
    kb_threshold = 10000
  )
  
  end_time <- Sys.time()
  extraction_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (!is.null(instruments) && nrow(instruments) > 0) {
    all_instruments[[trait_name]] <- instruments
    
    # 保存单个性状的工具变量
    save(instruments, 
         file = sprintf("results/instruments/%s_instruments.RData", trait_name))
    
    # 保存CSV
    tryCatch({
    write.csv(instruments, 
                sprintf("results/instruments/%s_instruments.csv", trait_name), 
              row.names = FALSE)
    }, error = function(e) {
      warning(sprintf("保存CSV失败: %s", conditionMessage(e)))
    })
    
    # 计算统计量
    mean_f_stat <- if("F_statistic" %in% names(instruments)) {
      round(mean(instruments$F_statistic, na.rm = TRUE), 2)
    } else {
      NA_real_
    }
    
    median_f_stat <- if("F_statistic" %in% names(instruments)) {
      round(median(instruments$F_statistic, na.rm = TRUE), 2)
    } else {
      NA_real_
    }
    
    min_f_stat <- if("F_statistic" %in% names(instruments)) {
      round(min(instruments$F_statistic, na.rm = TRUE), 2)
  } else {
      NA_real_
    }
    
    max_f_stat <- if("F_statistic" %in% names(instruments)) {
      round(max(instruments$F_statistic, na.rm = TRUE), 2)
    } else {
      NA_real_
    }
    
    weak_n <- if("weak_instrument" %in% names(instruments)) {
      sum(instruments$weak_instrument, na.rm = TRUE)
    } else {
      0
    }
    
    weak_prop <- if(nrow(instruments) > 0) {
      round(weak_n / nrow(instruments), 3)
    } else {
      NA_real_
    }
    
    var_exp <- if("variance_explained" %in% names(instruments)) {
      round(mean(instruments$variance_explained, na.rm = TRUE), 4)
    } else {
      NA_real_
    }
    
    # 获取GWAS ID（如果有）
    gwas_id <- if("id.exposure" %in% names(instruments)) {
      unique(instruments$id.exposure)[1]
    } else {
      NA_character_
    }
    
    # 添加到汇总
    instruments_summary <- rbind(instruments_summary, data.frame(
      trait = trait_name,
      gwas_id = gwas_id,
      n_instruments = nrow(instruments),
      mean_f_stat = mean_f_stat,
      median_f_stat = median_f_stat,
      min_f_stat = min_f_stat,
      max_f_stat = max_f_stat,
      weak_instrument_n = weak_n,
      weak_instrument_prop = weak_prop,
      total_variance_explained = var_exp,
      extraction_time_sec = round(extraction_time, 2),
      status = "Success",
      stringsAsFactors = FALSE
    ))
    
    extraction_log[[trait_name]] <- list(
      status = "Success",
      n_instruments = nrow(instruments),
      time = extraction_time,
      timestamp = Sys.time()
    )
    
    cat(sprintf("  ✓ 完成: %d 个工具变量 (耗时: %.1f秒)\n", 
                nrow(instruments), extraction_time))
    
  } else {
    # 提取失败
    instruments_summary <- rbind(instruments_summary, data.frame(
      trait = trait_name,
      gwas_id = NA_character_,
      n_instruments = 0,
      mean_f_stat = NA_real_,
      median_f_stat = NA_real_,
      min_f_stat = NA_real_,
      max_f_stat = NA_real_,
      weak_instrument_n = NA_integer_,
      weak_instrument_prop = NA_real_,
      total_variance_explained = NA_real_,
      extraction_time_sec = round(extraction_time, 2),
      status = "Failed",
      stringsAsFactors = FALSE
    ))
    
    extraction_log[[trait_name]] <- list(
      status = "Failed",
      n_instruments = 0,
      time = extraction_time,
      timestamp = Sys.time()
    )
    
    cat(sprintf("  ✗ 失败: 未提取到工具变量\n"))
  }
  
  # 定期保存进度
  if (current_trait %% 5 == 0) {
    save(all_instruments, 
         file = "data/step04_all_instruments_temp.RData")
    cat(sprintf("\n>>> 进度已保存：%d/%d 个暴露因子已处理 <<<\n", 
                current_trait, total_traits))
  }
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("✓ 工具变量提取完成: %d/%d 成功\n\n", 
            sum(instruments_summary$status == "Success"), total_traits))

# ===================================================================
# 4. 保存结果
# ===================================================================
cat("【步骤4】保存结果...\n")

# 保存所有工具变量
save(all_instruments, 
     file = "data/step04_all_instruments.RData")
cat("✓ 已保存: data/step04_all_instruments.RData\n")

# 保存汇总表
write.csv(instruments_summary, 
          "results/step04_instruments_summary.csv", 
          row.names = FALSE)
cat("✓ 已保存: results/step04_instruments_summary.csv\n")

# 保存提取日志
save(extraction_log, 
     file = "data/step04_extraction_log.RData")
cat("✓ 已保存: data/step04_extraction_log.RData\n\n")

# ===================================================================
# 5. 生成质量控制报告
# ===================================================================
cat("【步骤5】生成质量控制报告...\n")

success_count <- sum(instruments_summary$status == "Success")
failed_count <- sum(instruments_summary$status == "Failed")
success_rate <- if(total_traits > 0) {
  round(success_count / total_traits * 100, 1)
} else {
  0
}

total_instruments <- sum(instruments_summary$n_instruments, na.rm = TRUE)
mean_instruments <- if(success_count > 0) {
  round(mean(instruments_summary$n_instruments[instruments_summary$status == "Success"], na.rm = TRUE), 1)
} else {
  NA_real_
}

mean_f_stat <- if(success_count > 0) {
  round(mean(instruments_summary$mean_f_stat[instruments_summary$status == "Success"], na.rm = TRUE), 2)
    } else {
  NA_real_
}

mean_weak_prop <- if(success_count > 0) {
  round(mean(instruments_summary$weak_instrument_prop[instruments_summary$status == "Success"], na.rm = TRUE) * 100, 1)
    } else {
  NA_real_
}

total_time <- round(sum(instruments_summary$extraction_time_sec, na.rm = TRUE) / 60, 1)

# 生成报告文本
report_text <- paste(
  paste(rep("=", 80), collapse = ""),
  "        工具变量提取质量控制报告",
  paste(rep("=", 80), collapse = ""),
  "",
  "【总体统计】",
  paste(rep("-", 80), collapse = ""),
  sprintf("%-30s: %d", "总暴露因子数", total_traits),
  sprintf("%-30s: %d", "成功提取", success_count),
  sprintf("%-30s: %d", "提取失败", failed_count),
  sprintf("%-30s: %.1f%%", "成功率", success_rate),
  "",
  "【工具变量数量统计】",
  paste(rep("-", 80), collapse = ""),
  sprintf("%-30s: %d", "总工具变量数", total_instruments),
  if(!is.na(mean_instruments)) {
    sprintf("%-30s: %.1f", "每个暴露因子平均数量", mean_instruments)
  } else {
    sprintf("%-30s: %s", "每个暴露因子平均数量", "N/A")
  },
  "",
  "【工具变量强度（F统计量）】",
  paste(rep("-", 80), collapse = ""),
  if(!is.na(mean_f_stat)) {
    sprintf("%-30s: %.2f", "整体平均F统计量", mean_f_stat)
  } else {
    sprintf("%-30s: %s", "整体平均F统计量", "N/A")
  },
  if(!is.na(mean_weak_prop)) {
    sprintf("%-30s: %.1f%%", "平均弱工具变量比例", mean_weak_prop)
  } else {
    sprintf("%-30s: %s", "平均弱工具变量比例", "N/A")
  },
  "",
  "【提取时间】",
  paste(rep("-", 80), collapse = ""),
  sprintf("%-30s: %.1f 分钟", "总提取时间", total_time),
  "",
  paste(rep("=", 80), collapse = ""),
  "质量评估：",
  paste(rep("-", 80), collapse = ""),
  if(success_rate >= 80) {
    "✓ 提取成功率优秀 (≥80%)"
  } else if(success_rate >= 60) {
    "○ 提取成功率良好 (60-80%)"
  } else {
    "✗ 提取成功率较低 (<60%)，需要检查"
  },
  if(!is.na(mean_f_stat)) {
    if(mean_f_stat >= 30) {
      "✓ 工具变量强度优秀 (平均F>30)"
    } else if(mean_f_stat >= 10) {
      "○ 工具变量强度合格 (平均F>10)"
    } else {
      "✗ 工具变量强度不足 (平均F<10)，存在弱工具变量偏倚风险"
    }
  } else {
    "○ 无法评估工具变量强度（缺少F统计量）"
  },
  "",
  "详细结果请查看: results/step04_instruments_summary.csv",
  "单个性状工具变量保存在: results/instruments/ 目录",
  "",
  paste(rep("=", 80), collapse = ""),
  sep = "\n"
)

# 保存报告
writeLines(report_text, "results/reports/step04_qc_report.txt")
cat("✓ 已保存: results/reports/step04_qc_report.txt\n\n")

# ===================================================================
# 6. 生成可视化图表
# ===================================================================
cat("【步骤6】生成可视化图表...\n")

if (success_count > 0) {
  success_data <- instruments_summary[instruments_summary$status == "Success", ]
  
  # 工具变量数量分布图
  tryCatch({
  png("results/figures/step04_instruments_overview.png", 
      width = 14, height = 10, units = "in", res = 300)
    par(mfrow = c(2, 2), mar = c(5, 4, 3, 2))
  
  # 1. 工具变量数量分布
  hist(success_data$n_instruments, 
         main = "Distribution of Instrument Count",
       xlab = "Number of Instruments (SNPs)",
       ylab = "Frequency",
       col = "#3498db", border = "white",
       breaks = 20)
  abline(v = median(success_data$n_instruments), col = "red", lwd = 2, lty = 2)
  
  # 2. F统计量分布
    if(any(!is.na(success_data$mean_f_stat))) {
  hist(success_data$mean_f_stat,
       main = "Distribution of Mean F-statistics",
       xlab = "Mean F-statistic",
       ylab = "Frequency",
       col = "#2ecc71", border = "white",
       breaks = 20)
  abline(v = 10, col = "red", lwd = 2, lty = 2)
      abline(v = 30, col = "orange", lwd = 2, lty = 2)
    }
    
    # 3. 弱工具变量比例分布
    if(any(!is.na(success_data$weak_instrument_prop))) {
  hist(success_data$weak_instrument_prop * 100,
       main = "Weak Instrument Proportion",
       xlab = "Weak Instruments (%)",
       ylab = "Frequency",
       col = "#9b59b6", border = "white",
       breaks = 20)
  abline(v = 10, col = "red", lwd = 2, lty = 2)
    }
    
    # 4. F统计量vs工具变量数量散点图
    if(any(!is.na(success_data$mean_f_stat))) {
      plot(success_data$n_instruments, success_data$mean_f_stat,
           main = "Instrument Strength vs Count",
           xlab = "Number of Instruments",
           ylab = "Mean F-statistic",
           col = "#e74c3c", pch = 19, cex = 1.2)
      abline(h = 10, col = "red", lwd = 2, lty = 2)
  }
  
  dev.off()
    cat("✓ 已保存: results/figures/step04_instruments_overview.png\n")
  }, error = function(e) {
    warning(sprintf("图表生成失败: %s", conditionMessage(e)))
  })
  
  # 工具变量数量条形图
  tryCatch({
    png("results/figures/step04_instruments_barplot.png",
        width = 16, height = 8, units = "in", res = 300)
    par(mar = c(12, 4, 3, 2))
    
    # 按工具变量数量排序
    success_data_ordered <- success_data[order(success_data$n_instruments, decreasing = TRUE), ]
    
    barplot(success_data_ordered$n_instruments,
            names.arg = success_data_ordered$trait,
            las = 2,
            col = ifelse(success_data_ordered$mean_f_stat >= 30, "#2ecc71",
                        ifelse(success_data_ordered$mean_f_stat >= 10, "#f39c12", "#e74c3c")),
            main = "Number of Instruments per Exposure",
            ylab = "Number of Instruments",
            cex.names = 0.6)
  
  dev.off()
    cat("✓ 已保存: results/figures/step04_instruments_barplot.png\n")
  }, error = function(e) {
    warning(sprintf("条形图生成失败: %s", conditionMessage(e)))
  })
} else {
  warning("无成功提取的数据，跳过图表生成")
}

cat("\n")

# ===================================================================
# 7. 最终总结
# ===================================================================
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("步骤4完成！Instrument Variable Extraction\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("【处理结果汇总】\n")
cat(sprintf("  • 总暴露因子数: %d\n", total_traits))
cat(sprintf("  • 成功提取: %d (%.1f%%)\n", success_count, success_rate))
cat(sprintf("  • 提取失败: %d\n", failed_count))
cat(sprintf("  • 总工具变量数: %d\n", total_instruments))
if(!is.na(mean_instruments)) {
  cat(sprintf("  • 平均每个暴露因子: %.1f 个工具变量\n", mean_instruments))
}
if(!is.na(mean_f_stat)) {
  cat(sprintf("  • 平均F统计量: %.2f\n", mean_f_stat))
}
cat(sprintf("  • 总提取时间: %.1f 分钟\n\n", total_time))

cat("【生成的文件】\n\n")
cat("数据文件:\n")
cat("  ✓ data/step04_all_instruments.RData\n")
cat("  ✓ results/instruments/*_instruments.RData (单个性状工具变量)\n")
cat("  ✓ results/instruments/*_instruments.csv (单个性状工具变量CSV)\n\n")

cat("报告文件:\n")
cat("  ✓ results/step04_instruments_summary.csv\n")
cat("  ✓ results/reports/step04_qc_report.txt\n\n")

cat("可视化文件:\n")
if(file.exists("results/figures/step04_instruments_overview.png")) {
  cat("  ✓ results/figures/step04_instruments_overview.png\n")
}
if(file.exists("results/figures/step04_instruments_barplot.png")) {
  cat("  ✓ results/figures/step04_instruments_barplot.png\n")
}
cat("\n")

cat("【下一步操作】\n")
cat("  1. 查看工具变量汇总: results/step04_instruments_summary.csv\n")
cat("  2. 查看质量控制报告: results/reports/step04_qc_report.txt\n")
cat("  3. 检查数据质量: 查看各暴露因子的F统计量和弱工具变量比例\n")
cat("  4. 继续运行第5步: 单变量孟德尔随机化分析\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("工具变量提取完成！Ready for MR analysis (Step 5)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
