# ==============================================================================
# Script: 03_run_coloc.R (版本 19 - 修复版)
# Purpose: 运行共定位分析 (内存安全 + 10%检查点 + 恢复 + 修复版)
# 修复内容:
# 1. 修复eQTL数据MAF引用错误
# 2. 移除硬编码任务跳过
# 3. 增加数据验证机制
# ==============================================================================

cat("
================================================================================
Running Colocalization Analysis (Single-Core, Resumable, Error-Catching, FIXED)
================================================================================
\n")

# ==============================================================================
# 加载包
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(coloc)
  library(progress)   
  library(logger)
  library(arrow) 
  library(dplyr) 
})

# 设置日志
log_threshold(INFO)
log_appender(appender_tee("logs/03_run_coloc.log"))

log_info("Starting colocalization analysis (SINGLE-CORE, Resumable, Error-Catching, FIXED)")

# ==============================================================================
# 配置参数
# ==============================================================================

PROCESSED_DIR <- "data/processed" 
RESULTS_DIR <- "results"
TABLES_DIR <- file.path(RESULTS_DIR, "tables") 

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

COLOC_PARAMS <- list(
  p1 = 1e-4,    # 先验概率：trait 1有因果变异
  p2 = 1e-4,    # 先验概率：trait 2有因果变异  
  p12 = 1e-5    # 先验概率：共享因果变异
)

WINDOW_KB <- 500

log_info("Running in SINGLE-CORE mode to guarantee memory safety.")

# ==============================================================================
# 函数定义 (修复版)
# ==============================================================================

#' 数据验证函数
validate_data_slice <- function(data, slice_name) {
  required_cols <- c("chr", "pos", "pval", "beta", "se")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if(length(missing_cols) > 0) {
    stop(sprintf("%s missing required columns: %s", slice_name, paste(missing_cols, collapse = ", ")))
  }
  
  if(nrow(data) == 0) {
    stop(sprintf("%s is empty", slice_name))
  }
  
  # 检查MAF/eaf列
  if(!("eaf" %in% names(data)) && !("MAF" %in% names(data))) {
    log_warn(sprintf("%s missing eaf/MAF column, will use default values", slice_name))
  }
}

#' 提取基因区域数据
extract_gene_region <- function(gwas_data, eqtl_data, gene_symbol, window_kb = 500) {
  
  gene_info <- eqtl_data %>%
    filter(gene_symbol == !!gene_symbol) %>%
    distinct(chr = chromosome_name, start = start_position, end = end_position) %>%
    head(1)
  
  if (nrow(gene_info) == 0) {
    return(NULL)
  }
  
  chr <- gene_info$chr
  
  gwas_region <- gwas_data %>% filter(chr == !!chr)
  
  eqtl_region <- eqtl_data %>% filter(gene_symbol == !!gene_symbol, chr == !!chr)
  
  gwas_region <- gwas_region %>%
    mutate(snp_key = paste(chr, pos, sep = ":"))
  
  eqtl_region <- eqtl_region %>%
    mutate(snp_key = paste(chr, pos, sep = ":"))
  
  common_snps <- intersect(gwas_region$snp_key, eqtl_region$snp_key)
  
  if (length(common_snps) < 50) {
    return(NULL)
  }
  
  gwas_filtered <- gwas_region %>%
    filter(snp_key %in% common_snps) %>%
    arrange(chr, pos)
  
  eqtl_filtered <- eqtl_region %>%
    filter(snp_key %in% common_snps) %>%
    arrange(chr, pos)
  
  return(list(
    gene = gene_symbol,
    chr = chr,
    start = gene_info$start, 
    end = gene_info$end,
    n_snps = length(common_snps),
    gwas = gwas_filtered,
    eqtl = eqtl_filtered
  ))
}

#' 运行三路共定位 (修复版)
run_three_way_coloc <- function(exposure_data, eqtl_data, outcome_data, 
                                gene_symbol, exposure_name, outcome_name) {
  
  exposure_region <- extract_gene_region(exposure_data, eqtl_data, gene_symbol, WINDOW_KB)
  
  if (is.null(exposure_region)) {
    return(data.frame(
      gene = gene_symbol,
      exposure = exposure_name,
      outcome = outcome_name,
      status = "failed_extract_exposure",
      n_snps = 0
    ))
  }
  
  outcome_region <- extract_gene_region(outcome_data, eqtl_data, gene_symbol, WINDOW_KB)
  
  if (is.null(outcome_region)) {
    return(data.frame(
      gene = gene_symbol,
      exposure = exposure_name,
      outcome = outcome_name,
      status = "failed_extract_outcome",
      n_snps = exposure_region$n_snps
    ))
  }
  
  common_snps <- Reduce(intersect, list(
    exposure_region$gwas$snp_key,
    exposure_region$eqtl$snp_key,
    outcome_region$gwas$snp_key
  ))
  
  if (length(common_snps) < 50) {
    return(data.frame(
      gene = gene_symbol,
      exposure = exposure_name,
      outcome = outcome_name,
      status = "insufficient_snps",
      n_snps = length(common_snps)
    ))
  }
  
  exp_data <- exposure_region$gwas %>%
    filter(snp_key %in% common_snps) %>%
    arrange(snp_key)
  
  eqtl_data_filt <- exposure_region$eqtl %>%
    filter(snp_key %in% common_snps) %>%
    arrange(snp_key)
  
  out_data <- outcome_region$gwas %>%
    filter(snp_key %in% common_snps) %>%
    arrange(snp_key)
  
  tryCatch({
    
    # 数据验证
    validate_data_slice(exp_data, "Exposure data")
    validate_data_slice(eqtl_data_filt, "eQTL data") 
    validate_data_slice(out_data, "Outcome data")
    
    # 获取MAF值 (修复: 正确处理eQTL数据的MAF)
    get_maf_values <- function(data) {
      if("MAF" %in% names(data) && !is.null(data$MAF)) {
        return(data$MAF)
      } else if("eaf" %in% names(data) && !is.null(data$eaf)) {
        return(data$eaf)
      } else {
        # 使用默认中等频率
        return(rep(0.3, nrow(data)))
      }
    }
    
    dataset_exp <- list(
      pvalues = exp_data$pval,
      beta = exp_data$beta,
      varbeta = exp_data$se^2,
      N = median(exp_data$n, na.rm = TRUE),
      MAF = get_maf_values(exp_data),  # 修复: 正确获取MAF
      type = "quant",
      snp = exp_data$snp_key
    )
    
    dataset_eqtl <- list(
      pvalues = eqtl_data_filt$pval,
      beta = eqtl_data_filt$beta,
      varbeta = eqtl_data_filt$se^2,
      N = 515,
      MAF = get_maf_values(eqtl_data_filt),  # 修复: 使用eQTL自己的MAF
      type = "quant",
      snp = eqtl_data_filt$snp_key
    )
    
    coloc_exp_eqtl <- coloc.abf(
      dataset1 = dataset_exp,
      dataset2 = dataset_eqtl,
      p1 = COLOC_PARAMS$p1,
      p2 = COLOC_PARAMS$p2,
      p12 = COLOC_PARAMS$p12
    )
    
    dataset_out <- list(
      pvalues = out_data$pval,
      beta = out_data$beta,
      varbeta = out_data$se^2,
      N = median(out_data$n, na.rm = TRUE),
      MAF = get_maf_values(out_data),  # 修复: 正确获取MAF
      type = "cc",
      snp = out_data$snp_key
    )
    
    coloc_eqtl_out <- coloc.abf(
      dataset1 = dataset_eqtl,
      dataset2 = dataset_out,
      p1 = COLOC_PARAMS$p1,
      p2 = COLOC_PARAMS$p2,
      p12 = COLOC_PARAMS$p12
    )
    
    coloc_exp_out <- coloc.abf(
      dataset1 = dataset_exp,
      dataset2 = dataset_out,
      p1 = COLOC_PARAMS$p1,
      p2 = COLOC_PARAMS$p2,
      p12 = COLOC_PARAMS$p12
    )
    
    result <- data.frame(
      gene = gene_symbol,
      exposure = exposure_name,
      outcome = outcome_name,
      chr = exposure_region$chr,
      start = exposure_region$start,
      end = exposure_region$end,
      n_snps = length(common_snps),
      exp_eqtl_PP.H0 = coloc_exp_eqtl$summary["PP.H0.abf"],
      exp_eqtl_PP.H1 = coloc_exp_eqtl$summary["PP.H1.abf"],
      exp_eqtl_PP.H2 = coloc_exp_eqtl$summary["PP.H2.abf"],
      exp_eqtl_PP.H3 = coloc_exp_eqtl$summary["PP.H3.abf"],
      exp_eqtl_PP.H4 = coloc_exp_eqtl$summary["PP.H4.abf"],
      eqtl_out_PP.H0 = coloc_eqtl_out$summary["PP.H0.abf"],
      eqtl_out_PP.H1 = coloc_eqtl_out$summary["PP.H1.abf"],
      eqtl_out_PP.H2 = coloc_eqtl_out$summary["PP.H2.abf"],
      eqtl_out_PP.H3 = coloc_eqtl_out$summary["PP.H3.abf"],
      eqtl_out_PP.H4 = coloc_eqtl_out$summary["PP.H4.abf"],
      exp_out_PP.H0 = coloc_exp_out$summary["PP.H0.abf"],
      exp_out_PP.H1 = coloc_exp_out$summary["PP.H1.abf"],
      exp_out_PP.H2 = coloc_exp_out$summary["PP.H2.abf"],
      exp_out_PP.H3 = coloc_exp_out$summary["PP.H3.abf"],
      exp_out_PP.H4 = coloc_exp_out$summary["PP.H4.abf"],
      status = "success"
    )
    
    return(result)
    
  }, error = function(e) {
    # 增强的错误处理，包括内存不足检测
    error_msg <- gsub("[\r\n]", " ", e$message)
    
    # 如果是内存相关错误，标记为内存不足
    if(grepl("memory|out of memory|cannot allocate", error_msg, ignore.case = TRUE)) {
      return(data.frame(
        gene = gene_symbol,
        exposure = exposure_name,
        outcome = outcome_name,
        status = "error: memory_insufficient",
        n_snps = length(common_snps)
      ))
    }
    
    # 其他错误
    return(data.frame(
      gene = gene_symbol,
      exposure = exposure_name,
      outcome = outcome_name,
      status = paste("error:", error_msg),
      n_snps = length(common_snps)
    ))
  })
}

# ==============================================================================
# 主程序
# ==============================================================================

log_info(paste(rep("=", 80), collapse = ""))
log_info("Loading Data Paths")
log_info(paste(rep("=", 80), collapse = ""))

# --- 1. 获取 GWAS 文件路径 (不加载) ---
log_info("Finding GWAS data (Parquet format)...")
gwas_files <- list.files(file.path(PROCESSED_DIR, "gwas"), 
                         pattern = "_gwas.parquet$",  
                         full.names = TRUE)

gwas_data_paths <- gwas_files
names(gwas_data_paths) <- toupper(gsub("_gwas.parquet", "", basename(gwas_files)))

for (trait_name in names(gwas_data_paths)) {
  log_info(sprintf("  Found %s at: %s", trait_name, gwas_data_paths[[trait_name]]))
}

# --- 2. 获取 eQTL 文件路径 (不加载) ---
log_info("Finding eQTL data (Parquet format)...")
eqtl_data_path <- file.path(PROCESSED_DIR, "eqtl", "lung_eqtl_annotated.parquet")
if (!file.exists(eqtl_data_path)) {
  log_fatal(sprintf("eQTL file not found at: %s", eqtl_data_path))
  stop("eQTL file missing, check path and format.")
}
log_info(sprintf("  Found eQTL at: %s", eqtl_data_path))

# --- 3. 加载*小*的分析任务文件 (这必须加载) ---
tasks <- fread(file.path(PROCESSED_DIR, "analysis_tasks.txt"))
n_tasks <- nrow(tasks)
log_info(sprintf("  Loaded %s analysis tasks", 
                 format(n_tasks, big.mark = ",")))

# --- 4. 移除所有并行代码 ---
log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Running Colocalization Analysis (SINGLE-CORE / Memory-Safe Slicing)")
log_info(paste(rep("=", 80), collapse = ""))

# --- 5. 运行单核智能切片循环 ---
log_info("Starting single-core slicing loop...")

# <--- 检查点文件名
checkpoint_file <- file.path(TABLES_DIR, "coloc_results_TEMP.list.rds")

# <--- 计算 10% 的保存间隔
save_interval <- round(n_tasks * 0.10)
save_interval <- max(1, save_interval) 
log_info(sprintf("Checkpointing enabled. Will save every 10%% (approx. every %d tasks).", save_interval))

# <--- 恢复逻辑
start_index <- 1
if (file.exists(checkpoint_file)) {
  log_info("Found previous checkpoint file. Loading...")
  results_list <- readRDS(checkpoint_file)
  
  if (length(results_list) == n_tasks) {
    # 查找最后一个非 NULL (已完成) 的任务
    last_done <- tail(which(!sapply(results_list, is.null)), 1)
    
    if (length(last_done) == 0) {
      log_info("Checkpoint file was empty. Starting from task 1.")
    } else {
      start_index <- last_done + 1
      log_info(sprintf("Resuming from task %d / %d.", start_index, n_tasks))
    }
  } else {
    log_warning("Checkpoint file length mismatches task list. Starting from scratch.")
    results_list <- vector("list", n_tasks)
  }
} else {
  log_info("No checkpoint file found. Starting from scratch.")
  results_list <- vector("list", n_tasks)
}

# a. 打开 (不加载) Parquet 数据集 (在循环外执行一次)
ds_eqtl <- arrow::open_dataset(eqtl_data_path)
ds_exp_list <- lapply(gwas_data_paths, arrow::open_dataset)
names(ds_exp_list) <- names(gwas_data_paths)

# b. 启动进度条 (从已完成的开始)
pb <- progress_bar$new(
  format = "  [:bar] :current/:total (:percent) | ETA: :eta",
  total = n_tasks,
  clear = FALSE
)
pb$tick(start_index - 1) # 快进到我们已完成的进度

# c. 运行标准 for 循环 (从 start_index 开始)
if (start_index <= n_tasks) {
  for (i in start_index:n_tasks) {
    
    task <- tasks[i, ]
    
    # <--- 移除硬编码跳过，改为通用错误处理
    # 不再硬编码跳过特定任务
    
    # 增强的错误处理
    results_list[[i]] <- tryCatch({
      
      # a. 从任务 中获取基因坐标
      gene_chr <- task$chr
      gene_start <- task$start
      gene_end <- task$end
      
      # b. 定义分析窗口
      window_start <- max(1, gene_start - WINDOW_KB * 1000)
      window_end <- gene_end + WINDOW_KB * 1000
      
      # c. 匹配键名 (使用 V12 的修复)
      exp_key <- toupper(task$exposure)
      out_key <- toupper(task$outcome)
      
      # d. 获取已打开的数据集
      ds_exp <- ds_exp_list[[exp_key]]
      ds_out <- ds_exp_list[[out_key]] 
      
      # e. 安全检查
      if (is.null(ds_exp) || is.null(ds_out)) {
        stop(sprintf("GWAS dataset not found (exp: %s, out: %s)", exp_key, out_key))
      }
      
      # f. *只*加载窗口内的 eQTL 数据切片
      eqtl_slice <- ds_eqtl %>%
        filter(chr == gene_chr, 
               pos >= window_start, 
               pos <= window_end) %>%
        collect()
      
      # g. *只*加载窗口内的 Exposure GWAS 数据切片
      exp_slice <- ds_exp %>%
        filter(chr == gene_chr, 
               pos >= window_start, 
               pos <= window_end) %>%
        collect()
      
      # h. *只*加载窗口内的 Outcome GWAS 数据切片
      out_slice <- ds_out %>%
        filter(chr == gene_chr, 
               pos >= window_start, 
               pos <= window_end) %>%
        collect()
      
      # i. 检查切片是否为空
      if (nrow(eqtl_slice) == 0 || nrow(exp_slice) == 0 || nrow(out_slice) == 0) {
        data.frame(
          gene = task$gene_symbol,
          exposure = task$exposure,
          outcome = task$outcome,
          status = "insufficient_snps_in_slice",
          n_snps = 0
        )
      } else {
        # j. 将*小切片*传递给函数
        run_three_way_coloc(
          exposure_data = exp_slice,
          eqtl_data = eqtl_slice,
          outcome_data = out_slice,
          gene_symbol = task$gene_symbol, 
          exposure_name = task$exposure,
          outcome_name = task$outcome
        )
      }
      
    }, error = function(e) {
      # 增强的错误处理
      error_message <- gsub("[\r\n]", " ", e$message)
      log_error(sprintf("Task %d (%s) failed: %s", i, task$gene_symbol, error_message))
      
      data.frame(
        gene = task$gene_symbol,
        exposure = task$exposure,
        outcome = task$outcome,
        status = paste("error:", error_message),
        n_snps = 0
      )
    })
    
    # 更新进度条
    pb$tick()
    
    # 检查点 (CHECKPOINT)
    if (i %% save_interval == 0) {
      log_info(sprintf("\n[CHECKPOINT] Saving list... (Task %d)", i))
      saveRDS(results_list, checkpoint_file)
      log_info("[CHECKPOINT] Save complete. Resuming analysis...")
    }
    
  }
} else {
  log_info("All tasks were already completed in the checkpoint file.")
}

# --- 6. 循环结束 ---

log_info("\nCombining final results...")
results_df <- rbindlist(results_list, fill = TRUE)

log_info("\n✓ Colocalization analysis complete!")

# 保存完整结果
results_file <- file.path(TABLES_DIR, "coloc_results_full.txt")
fwrite(results_df, results_file, sep = "\t")
log_info(sprintf("  Saved final results to %s", results_file))

saveRDS(results_df, file.path(TABLES_DIR, "coloc_results_full.rds"))

log_info(sprintf("  Success! Temporary file '%s' is no longer needed.", basename(checkpoint_file)))

# ==============================================================================
# 结果摘要 (无改动)
# ==============================================================================

log_info(paste0("\n", paste(rep("=", 80), collapse = "")))
log_info("Analysis Summary")
log_info(paste(rep("=", 80), collapse = ""))

# 总体统计
cat(sprintf("\nTotal analyses: %d\n", nrow(results_df)))
cat(sprintf("Successful: %d (%.1f%%)\n",
            sum(results_df$status == "success"),
            100 * mean(results_df$status == "success")))

# 成功的分析
success_results <- results_df %>%
  filter(status == "success")

cat(sprintf("\nSuccessful analyses: %d\n", nrow(success_results)))

# 强共定位 (PP.H4 > 0.8)
strong_coloc <- success_results %>%
  filter(exp_eqtl_PP.H4 > 0.8 | eqtl_out_PP.H4 > 0.8)

cat(sprintf("Strong colocalization (PP.H4 > 0.8): %d\n", nrow(strong_coloc)))

if (nrow(strong_coloc) > 0) {
  cat("\nTop colocalized genes:\n")
  top_genes <- strong_coloc %>%
    select(gene, exposure, outcome, exp_eqtl_PP.H4, eqtl_out_PP.H4) %>%
    arrange(desc(pmax(exp_eqtl_PP.H4, eqtl_out_PP.H4))) %>%
    head(10)
  print(top_genes)
}

# 按暴露因素统计
cat("\nResults by exposure:\n")
exposure_summary <- success_results %>%
  group_by(exposure) %>%
  summarise(
    n_analyses = n(),
    n_strong_exp_sql = sum(exp_eqtl_PP.H4 > 0.8),
    n_strong_eqtl_out = sum(eqtl_out_PP.H4 > 0.8),
    mean_PP.H4_exp_eqtl = mean(exp_eqtl_PP.H4),
    mean_PP.H4_eqtl_out = mean(eqtl_out_PP.H4)
  )
print(exposure_summary)

# 保存摘要
fwrite(exposure_summary, 
       file.path(RESULTS_DIR, "tables", "coloc_summary_by_exposure.txt"),
       sep = "\t")

log_info("\n✓ Analysis complete!")
log_info("Next step: Summarize and visualize results (04_summarize_results.R)")