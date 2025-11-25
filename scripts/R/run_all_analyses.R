#!/usr/bin/env Rscript
# 主表1多组学数据补充 - 主控脚本
# 自动运行所有4个分析步骤

cat("\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║     主表1：50个高质量共定位基因多组学概览                 ║\n")
cat("║              数据补充分析流程                              ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat("\n")

`%R%` <- function(x, times) paste(rep(x, times), collapse = "")

# 设置工作目录
if (!grepl("主表1_多组学整合$", normalizePath(getwd(), winslash = "/"))) setwd("论文图表汇总/主图/主表1_多组学整合")

# 检查必要的包
required_packages <- c("biomaRt", "data.table", "dplyr", "limma", "survival", "survminer", "ggplot2", "cowplot")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  cat("警告：以下R包未安装:\n")
  cat(paste("  -", missing_packages, collapse = "\n"), "\n\n")
  cat("请运行以下命令安装:\n")
  cat("  install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n\n", sep = "")
  cat("对于Bioconductor包:\n")
  cat("  BiocManager::install(c('biomaRt', 'limma'))\n\n")
  stop("缺少必要的R包")
}

cat("✓ 所有必要的R包已安装\n\n")

# 分析步骤
steps <- list(
  list(
    name = "步骤1：基因符号注释",
    script = "step1_annotate_genes.R",
    description = "将ENSG ID转换为HGNC基因符号"
  ),
  list(
    name = "步骤2：甲基化差异分析",
    script = "step2_methylation_analysis.R",
    description = "分析启动子区域甲基化差异"
  ),
  list(
    name = "步骤3：蛋白组差异分析",
    script = "step3_proteome_analysis.R",
    description = "计算蛋白水平差异表达"
  ),
  list(
    name = "步骤4：生存分析",
    script = "step4_survival_analysis.R",
    description = "计算LUSC患者生存风险比"
  )
)

# 记录开始时间
start_time <- Sys.time()
results <- list()

# 运行每个步骤
for (i in seq_along(steps)) {
  step <- steps[[i]]
  
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat(step$name, "\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("描述:", step$description, "\n")
  cat("脚本:", step$script, "\n\n")
  
  if (!file.exists(step$script)) {
    cat("✗ 脚本文件不存在:", step$script, "\n")
    results[[i]] <- list(
      step = step$name,
      status = "失败",
      message = "脚本文件不存在"
    )
    next
  }
  
  # 运行脚本
  step_start <- Sys.time()
  tryCatch({
    source(step$script, echo = FALSE)
    step_end <- Sys.time()
    elapsed <- as.numeric(difftime(step_end, step_start, units = "secs"))
    
    cat("\n✓", step$name, "完成\n")
    cat("  用时:", round(elapsed, 1), "秒\n")
    
    results[[i]] <- list(
      step = step$name,
      status = "成功",
      time = elapsed
    )
    
  }, error = function(e) {
    cat("\n✗", step$name, "失败\n")
    cat("  错误:", e$message, "\n")
    
    results[[i]] <- list(
      step = step$name,
      status = "失败",
      message = e$message
    )
  })
}

# 记录结束时间
end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

# 输出总结
cat("\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║                     分析完成总结                           ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat("\n")

successful <- sum(sapply(results, function(x) x$status == "成功"))
failed <- length(results) - successful

cat("总步骤数:", length(steps), "\n")
cat("成功:", successful, "\n")
cat("失败:", failed, "\n")
cat("总用时:", round(total_time, 2), "分钟\n\n")

# 详细结果
cat("详细结果:\n")
cat("────────────────────────────────────────────────────────────\n")
for (result in results) {
  status_symbol <- ifelse(result$status == "成功", "✓", "✗")
  cat(sprintf(
    "%s %s - %s", 
    status_symbol, 
    result$step,
    result$status
  ))
  if (!is.null(result$time)) {
    cat(sprintf(" (%.1f秒)", result$time))
  }
  cat("\n")
  if (!is.null(result$message)) {
    cat("   ", result$message, "\n")
  }
}
cat("────────────────────────────────────────────────────────────\n\n")

# 检查最终输出文件
final_file <- "主表1_核心版_明星基因证据链_完整版.csv"
if (file.exists(final_file)) {
  cat("✓ 最终主表文件已生成:\n")
  cat("  ", final_file, "\n\n")
  
  # 读取并显示摘要
  final_table <- data.table::fread(final_file)
  cat("主表摘要:\n")
  cat("  基因数:", nrow(final_table), "\n")
  cat("  数据列:", ncol(final_table), "\n")
  
  # 检查数据完整性
  tbd_count <- sum(grepl("TBD|待|无数据", as.matrix(final_table)))
  complete_count <- prod(dim(final_table)) - tbd_count
  completeness <- round(100 * complete_count / prod(dim(final_table)), 1)
  
  cat("  数据完整度:", completeness, "%\n")
  cat("  待补充项:", tbd_count, "\n\n")
  out_dir <- getwd()
  en_file_csv <- file.path(out_dir, "MasterTable_Final_Manuscript.csv")
  en_file_md <- file.path(out_dir, "MasterTable_Final_Manuscript.md")
  ft <- as.data.frame(final_table)
  colnames(ft) <- c("Gene", "MR Exposure", "Colocalization Score", "Promoter Methylation (LUSC vs LUAD)", "RNA Fold Change (LUSC vs LUAD)", "Protein Fold Change (LUSC vs LUAD)", "LUSC OS HR")
  repl <- function(x){
    x <- gsub("待补充", "Pending", x)
    x <- gsub("待完整分析", "Pending full analysis", x)
    x <- gsub("无数据", "No data", x)
    x <- gsub("无差异", "No difference", x)
    x <- gsub("LUSC↑", "LUSC up", x)
    x <- gsub("LUSC↓", "LUSC down", x)
    x <- gsub("LUSC高甲基化", "LUSC hypermethylated", x)
    x <- gsub("LUAD高甲基化", "LUAD hypermethylated", x)
    x <- gsub("炎症", "inflammation", x)
    x <- gsub("代谢", "metabolism", x)
    x
  }
  ft[] <- lapply(ft, function(col){ if (is.character(col)) repl(col) else col })
  ftGene <- gsub("\\.\\d+$", "", ft$Gene)
  suppressWarnings({ ft[["Colocalization Score"]] <- round(as.numeric(ft[["Colocalization Score"]]), 3) })
  if (file.exists("methylation_analysis_results.csv")){
    mr <- data.table::fread("methylation_analysis_results.csv")
    names(mr) <- trimws(names(mr))
    mr$Gene <- gsub("\\.\\d+$", "", mr$Gene)
    for (i in seq_along(ftGene)){
      j <- match(ftGene[i], mr$Gene)
      if (!is.na(j)){
        ft[i, "Promoter Methylation (LUSC vs LUAD)"] <- paste0(
          "Δβ=", sprintf("%.3f", mr$Delta_Beta[j]), " (",
          ifelse(mr$Delta_Beta[j] > 0.1, "LUSC hypermethylated", ifelse(mr$Delta_Beta[j] < -0.1, "LUAD hypermethylated", "No difference")), ")",
          ifelse(!is.na(mr$P_value[j]), paste0("; P=", sprintf("%.3f", mr$P_value[j])), ""),
          ifelse(!is.na(mr$FDR[j]), paste0("; FDR=", sprintf("%.3f", mr$FDR[j])), "")
        )
      }
    }
  }
  if (file.exists("proteome_analysis_results.csv")){
    pr <- data.table::fread("proteome_analysis_results.csv")
    names(pr) <- trimws(names(pr))
    pr$Gene <- gsub("\\.\\d+$", "", pr$Gene)
    for (i in seq_along(ftGene)){
      j <- match(ftGene[i], pr$Gene)
      if (!is.na(j)){
        ft[i, "Protein Fold Change (LUSC vs LUAD)"] <- paste0(
          sprintf("%.3f", pr$Log2FC[j]), " (",
          ifelse(pr$Log2FC[j] > 0, "LUSC up", ifelse(pr$Log2FC[j] < 0, "LUAD up", "No difference")), ")",
          ifelse(!is.na(pr$P_value[j]), paste0("; P=", sprintf("%.3f", pr$P_value[j])), ""),
          ifelse(!is.na(pr$FDR[j]), paste0("; FDR=", sprintf("%.3f", pr$FDR[j])), "")
        )
      }
    }
  }
  if (file.exists("survival_analysis_results.csv")){
    sr <- data.table::fread("survival_analysis_results.csv")
    names(sr) <- trimws(names(sr))
    sr$Gene <- gsub("\\.\\d+$", "", sr$Gene)
    for (i in seq_along(ftGene)){
      j <- match(ftGene[i], sr$Gene)
      if (!is.na(j)){
        ft[i, "LUSC OS HR"] <- paste0(
          sprintf("%.3f", sr$HR[j]), " (", sprintf("%.3f", sr$HR_Lower[j]), "–", sprintf("%.3f", sr$HR_Upper[j]), ")",
          ifelse(!is.na(sr$P_value[j]), paste0("; P=", sprintf("%.3f", sr$P_value[j])), "")
        )
      }
    }
  }

  # RNA差异分析（log2(TPM+1)；t检验；输出P/FDR并同步主表）
  root_dir_local <- normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  tcga_dir_local <- file.path(root_dir_local, "PDC", "TCGA")
  luad_expr_file <- if (file.exists(file.path(tcga_dir_local, "TCGA-LUAD.star_tpm.tsv.gz"))) file.path(tcga_dir_local, "TCGA-LUAD.star_tpm.tsv.gz") else file.path(tcga_dir_local, "TCGA-LUAD.star_tpm.tsv")
  lusc_expr_file <- if (file.exists(file.path(tcga_dir_local, "TCGA-LUSC.star_tpm.tsv.gz"))) file.path(tcga_dir_local, "TCGA-LUSC.star_tpm.tsv.gz") else file.path(tcga_dir_local, "TCGA-LUSC.star_tpm.tsv")
  if (file.exists(luad_expr_file) && file.exists(lusc_expr_file)){
    luad_ex <- data.table::fread(luad_expr_file)
    lusc_ex <- data.table::fread(lusc_expr_file)
    luad_ex[[1]] <- sub("\\.\\d+$", "", luad_ex[[1]])
    lusc_ex[[1]] <- sub("\\.\\d+$", "", lusc_ex[[1]])
    annmap_file <- file.path(getwd(), "gene_annotations_complete.csv")
    ens_ids <- ifelse(grepl("^ENSG", ftGene), ftGene, NA)
    if (file.exists(annmap_file)){
      ann <- data.table::fread(annmap_file)
      ann$ensembl_gene_id <- sub("\\.\\d+$", "", ann$ensembl_gene_id)
      for (i in seq_along(ftGene)){
        if (is.na(ens_ids[i])){
          hit <- ann$ensembl_gene_id[ann$hgnc_symbol == ftGene[i]]
          ens_ids[i] <- if (length(hit) > 0) hit[1] else NA
        }
      }
    }
    rna_res <- data.table::data.table(Gene = ftGene, Log2FC = NA_real_, P_value = NA_real_, FDR = NA_real_, Direction = NA_character_)
    for (i in seq_along(ftGene)){
      eid <- ens_ids[i]
      if (is.na(eid)) next
      r1 <- luad_ex[luad_ex[[1]] == eid, ]
      r2 <- lusc_ex[lusc_ex[[1]] == eid, ]
      if (nrow(r1) == 0 || nrow(r2) == 0) next
      v1 <- log2(as.numeric(r1[, -1]) + 1)
      v2 <- log2(as.numeric(r2[, -1]) + 1)
      if (length(v1) < 3 || length(v2) < 3) next
      luad_mean <- mean(v1, na.rm = TRUE)
      lusc_mean <- mean(v2, na.rm = TRUE)
      log2fc <- lusc_mean - luad_mean
      test_result <- tryCatch(stats::t.test(v2, v1), error = function(e) NULL)
      pval <- if (!is.null(test_result)) test_result$p.value else NA_real_
      rna_res$Log2FC[i] <- log2fc
      rna_res$P_value[i] <- pval
      rna_res$Direction[i] <- ifelse(log2fc > 0, "LUSC up", ifelse(log2fc < 0, "LUAD up", "No difference"))
    }
    if (sum(!is.na(rna_res$P_value)) > 0){
      rna_res$FDR <- p.adjust(rna_res$P_value, method = "BH")
    }
    data.table::fwrite(rna_res, "rna_analysis_results.csv")
    for (i in seq_along(ftGene)){
      j <- match(ftGene[i], rna_res$Gene)
      if (!is.na(j) && !is.na(rna_res$Log2FC[j])){
        ft[i, "RNA Fold Change (LUSC vs LUAD)"] <- paste0(
          sprintf("%.3f", rna_res$Log2FC[j]), " (",
          rna_res$Direction[j], ")",
          ifelse(!is.na(rna_res$P_value[j]), paste0("; P=", sprintf("%.3f", rna_res$P_value[j])), ""),
          ifelse(!is.na(rna_res$FDR[j]), paste0("; FDR=", sprintf("%.3f", rna_res$FDR[j])), "")
        )
      }
    }
  }
  data.table::fwrite(ft, en_file_csv)
  hdr <- paste0("| ", paste(colnames(ft), collapse = " | "), " |\n")
  sep <- paste0("| ", paste(rep("---", ncol(ft)), collapse = " | "), " |\n")
  rows <- apply(ft, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |\n"))
  cat(hdr, sep, paste(rows, collapse = ""), file = en_file_md)
  cat("✓ 英文化主表已生成:\n  ", en_file_csv, "\n  ", en_file_md, "\n\n", sep = "")

  review_csv <- file.path(out_dir, "MasterTable_For_Review.csv")
  fr <- data.table::data.table(
    Gene = ft$Gene,
    `MR Exposure` = ft$`MR Exposure`,
    `Colocalization Score` = suppressWarnings(as.numeric(ft$`Colocalization Score`)),
    `Promoter Methylation DeltaBeta` = NA_real_,
    `Promoter Methylation Direction` = NA_character_,
    `Promoter Methylation P` = NA_real_,
    `Promoter Methylation FDR` = NA_real_,
    `Protein Log2FC` = NA_real_,
    `Protein Direction` = NA_character_,
    `Protein P` = NA_real_,
    `Protein FDR` = NA_real_,
    `RNA Log2FC` = NA_real_,
    `RNA Direction` = NA_character_,
    `RNA P` = NA_real_,
    `RNA FDR` = NA_real_,
    `LUSC OS HR` = NA_real_,
    `LUSC OS P` = NA_real_
  )
  gkey <- gsub("\\.\\d+$", "", final_table$Gene_Symbol)
  if (file.exists("methylation_analysis_results.csv")){
    mr <- data.table::fread("methylation_analysis_results.csv")
    mr$Gene <- gsub("\\.\\d+$", "", mr$Gene)
    mm <- setNames(seq_len(nrow(mr)), mr$Gene)
    for (i in seq_along(gkey)){
      k <- gkey[i]
      j <- mm[k]
      if (!is.na(j)){
        fr$`Promoter Methylation DeltaBeta`[i] <- mr$Delta_Beta[j]
        fr$`Promoter Methylation Direction`[i] <- mr$Direction[j]
        fr$`Promoter Methylation P`[i] <- mr$P_value[j]
        fr$`Promoter Methylation FDR`[i] <- mr$FDR[j]
      }
    }
  }
  if (file.exists("proteome_analysis_results.csv")){
    pr <- data.table::fread("proteome_analysis_results.csv")
    pr$Gene <- gsub("\\.\\d+$", "", pr$Gene)
    pcol <- if ("P" %in% names(pr)) "P" else if ("P_value" %in% names(pr)) "P_value" else NA_character_
    mm <- setNames(seq_len(nrow(pr)), pr$Gene)
    for (i in seq_along(gkey)){
      k <- gkey[i]
      j <- mm[k]
      if (!is.na(j)){
        fr$`Protein Log2FC`[i] <- pr$Log2FC[j]
        fr$`Protein Direction`[i] <- ifelse(pr$Log2FC[j] > 0, "LUSC up", ifelse(pr$Log2FC[j] < 0, "LUAD up", "No difference"))
        if (!is.na(pcol)) fr$`Protein P`[i] <- pr[[pcol]][j]
        if ("FDR" %in% names(pr)) fr$`Protein FDR`[i] <- pr$FDR[j]
      }
    }
  }
  if (file.exists("rna_analysis_results.csv")){
    rr <- data.table::fread("rna_analysis_results.csv")
    rr$Gene <- gsub("\\.\\d+$", "", rr$Gene)
    mm <- setNames(seq_len(nrow(rr)), rr$Gene)
    for (i in seq_along(gkey)){
      k <- gkey[i]
      j <- mm[k]
      if (!is.na(j)){
        fr$`RNA Log2FC`[i] <- rr$Log2FC[j]
        fr$`RNA Direction`[i] <- rr$Direction[j]
        fr$`RNA P`[i] <- rr$P_value[j]
        fr$`RNA FDR`[i] <- rr$FDR[j]
      }
    }
  }
  if (file.exists("survival_analysis_results.csv")){
    sr <- data.table::fread("survival_analysis_results.csv")
    sr$Gene <- gsub("\\.\\d+$", "", sr$Gene)
    mm <- setNames(seq_len(nrow(sr)), sr$Gene)
    for (i in seq_along(gkey)){
      k <- gkey[i]
      j <- mm[k]
      if (!is.na(j)){
        fr$`LUSC OS HR`[i] <- sr$HR[j]
        fr$`LUSC OS P`[i] <- sr$P_value[j]
      }
    }
  }
  fr$`Colocalization Score` <- round(fr$`Colocalization Score`, 3)
  fr$`Promoter Methylation DeltaBeta` <- round(fr$`Promoter Methylation DeltaBeta`, 3)
  fr$`Protein Log2FC` <- round(fr$`Protein Log2FC`, 3)
  fr$`RNA Log2FC` <- round(fr$`RNA Log2FC`, 3)
  fr$`LUSC OS HR` <- round(fr$`LUSC OS HR`, 3)
  fr$`Promoter Methylation P` <- ifelse(is.na(fr$`Promoter Methylation P`), "", sprintf("%.3f", fr$`Promoter Methylation P`))
  fr$`Promoter Methylation FDR` <- ifelse(is.na(fr$`Promoter Methylation FDR`), "", sprintf("%.3f", fr$`Promoter Methylation FDR`))
  fr$`Protein P` <- ifelse(is.na(fr$`Protein P`), "", sprintf("%.3f", fr$`Protein P`))
  fr$`Protein FDR` <- ifelse(is.na(fr$`Protein FDR`), "", sprintf("%.3f", fr$`Protein FDR`))
  fr$`RNA P` <- ifelse(is.na(fr$`RNA P`), "", sprintf("%.3f", fr$`RNA P`))
  fr$`RNA FDR` <- ifelse(is.na(fr$`RNA FDR`), "", sprintf("%.3f", fr$`RNA FDR`))
  fr$`LUSC OS P` <- ifelse(is.na(fr$`LUSC OS P`), "", sprintf("%.3f", fr$`LUSC OS P`))
  data.table::fwrite(fr, review_csv)
  cat("✓ 审稿专用主表已生成:\n  ", review_csv, "\n\n", sep = "")
  
} else {
  cat("✗ 最终主表文件未生成\n")
  cat("  请检查各步骤的执行情况\n\n")
}

# 生成分析报告
report_file <- "analysis_execution_report.txt"
sink(report_file)
cat("主表1多组学数据补充分析报告\n")
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("分析步骤执行情况:\n")
for (result in results) {
  cat("\n", result$step, ":\n", sep = "")
  cat("  状态:", result$status, "\n")
  if (!is.null(result$time)) {
    cat("  用时:", round(result$time, 1), "秒\n")
  }
  if (!is.null(result$message)) {
    cat("  信息:", result$message, "\n")
  }
}
cat("\n总用时:", round(total_time, 2), "分钟\n")
sink()

cat("✓ 分析报告已保存:", report_file, "\n\n")

cat("════════════════════════════════════════════════════════════\n")
cat("所有分析流程执行完毕！\n")
cat("════════════════════════════════════════════════════════════\n\n")
plots_dir <- "plots"
dir.create(file.path(plots_dir, "methylation"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(plots_dir, "protein"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(plots_dir, "survival"), showWarnings = FALSE, recursive = TRUE)
root_dir <- normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
tcga_dir <- file.path(root_dir, "PDC", "TCGA")
ann_file <- file.path(root_dir, "PDC", "output", "V6_methylation", "450k_cpg_to_gene.tsv")
merged_file <- file.path(root_dir, "PDC", "output", "V6_methylation", "merged_beta.tsv")
if (file.exists(ann_file) && file.exists(merged_file)){
  cpg <- data.table::fread(ann_file)
  beta_raw <- data.table::fread(merged_file, header = FALSE, data.table = FALSE)
  sids <- as.character(beta_raw[1, -1])
  beta_raw <- beta_raw[-1, ]
  colnames(beta_raw) <- c("Probe", sids)
  rn <- beta_raw$Probe
  beta_raw$Probe <- NULL
  beta <- as.matrix(beta_raw)
  rownames(beta) <- rn
  luad_hdr <- tryCatch(data.table::fread(file.path(tcga_dir, "TCGA-LUAD.methylation450.tsv"), nrows = 0, data.table = FALSE), error = function(e) NULL)
  lusc_hdr <- tryCatch(data.table::fread(file.path(tcga_dir, "TCGA-LUSC.methylation450.tsv"), nrows = 0, data.table = FALSE), error = function(e) NULL)
  luad_sids <- if (!is.null(luad_hdr)) colnames(luad_hdr)[-1] else character(0)
  lusc_sids <- if (!is.null(lusc_hdr)) colnames(lusc_hdr)[-1] else character(0)
  is_tumor <- function(x) grepl("-01A$", x)
  luad_cols <- intersect(sids[is_tumor(sids)], luad_sids)
  lusc_cols <- intersect(sids[is_tumor(sids)], lusc_sids)
  mf <- data.table::fread(final_file)
  gl <- gsub("\\.\\d+$", "", mf$Gene_Symbol)
  for (g in unique(gl)){
    pro <- cpg$Probe[grepl(paste0("(^|;)", g, "(;|$)"), cpg$GeneSymbol)]
    keep <- intersect(pro, rownames(beta))
    if (length(keep) == 0) next
    sub <- beta[keep, , drop = FALSE]
    x1 <- colMeans(sub[, colnames(sub) %in% luad_cols, drop = FALSE], na.rm = TRUE)
    x2 <- colMeans(sub[, colnames(sub) %in% lusc_cols, drop = FALSE], na.rm = TRUE)
    df <- data.frame(group = c(rep("LUAD", length(x1)), rep("LUSC", length(x2))), value = c(x1, x2))
    pv <- tryCatch(wilcox.test(x1, x2)$p.value, error = function(e) NA_real_)
    lab <- paste0("Wilcoxon P=", signif(pv, 3))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, fill = group)) + ggplot2::geom_boxplot(outlier.size = 0.5) + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15), size = 0.8, alpha = 0.4) + ggplot2::scale_fill_manual(values = c(LUAD = "#1f77b4", LUSC = "#d62728")) + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = paste0(g, " promoter methylation (LUSC vs LUAD)"), x = NULL, y = "Beta value") + ggplot2::annotate("label", x = 2, y = max(df$value, na.rm = TRUE), label = lab, hjust = 1, vjust = 1, size = 3)
    ggplot2::ggsave(file.path(plots_dir, "methylation", paste0(g, "_methylation_boxplot.png")), p, width = 6, height = 4, dpi = 600, bg = "transparent")
    ggplot2::ggsave(file.path(plots_dir, "methylation", paste0(g, "_methylation_boxplot.pdf")), p, width = 6, height = 4)
    ggplot2::ggsave(file.path(plots_dir, "methylation", paste0(g, "_methylation_boxplot.tiff")), p, width = 6, height = 4, dpi = 600, compression = "lzw")
  }
}
luad_pf <- file.path(tcga_dir, "TCGA-LUAD.protein.tsv")
lusc_pf <- file.path(tcga_dir, "TCGA-LUSC.protein.tsv")
if (file.exists(luad_pf) && file.exists(lusc_pf)){
  lu <- data.table::fread(luad_pf)
  ls <- data.table::fread(lusc_pf)
  gc <- colnames(lu)[1]
  mf <- data.table::fread(final_file)
  gl <- gsub("\\.\\d+$", "", mf$Gene_Symbol)
  for (g in unique(gl)){
    v1 <- as.numeric(unlist(lu[lu[[gc]] == g, -1]))
    v2 <- as.numeric(unlist(ls[ls[[gc]] == g, -1]))
    if (length(v1) == 0 && length(v2) == 0) next
    df <- data.frame(group = c(rep("LUAD", length(v1)), rep("LUSC", length(v2))), value = c(v1, v2))
    pv <- tryCatch(stats::t.test(v2, v1)$p.value, error = function(e) NA_real_)
    lab <- paste0("t-test P=", signif(pv, 3))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = value, fill = group)) + ggplot2::geom_boxplot(outlier.size = 0.5) + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15), size = 0.8, alpha = 0.4) + ggplot2::scale_fill_manual(values = c(LUAD = "#1f77b4", LUSC = "#d62728")) + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = paste0(g, " protein abundance (LUSC vs LUAD)"), x = NULL, y = "Abundance") + ggplot2::annotate("label", x = 2, y = max(df$value, na.rm = TRUE), label = lab, hjust = 1, vjust = 1, size = 3)
    ggplot2::ggsave(file.path(plots_dir, "protein", paste0(g, "_protein_boxplot.png")), p, width = 6, height = 4, dpi = 600, bg = "transparent")
    ggplot2::ggsave(file.path(plots_dir, "protein", paste0(g, "_protein_boxplot.pdf")), p, width = 6, height = 4)
    ggplot2::ggsave(file.path(plots_dir, "protein", paste0(g, "_protein_boxplot.tiff")), p, width = 6, height = 4, dpi = 600, compression = "lzw")
  }
}
lusc_expr_file <- file.path(tcga_dir, "TCGA-LUSC.star_tpm.tsv.gz")
lusc_surv_file <- file.path(tcga_dir, "TCGA-LUSC.survival.tsv.gz")
if (file.exists(lusc_expr_file) && file.exists(lusc_surv_file)){
  ex <- data.table::fread(lusc_expr_file)
  su <- data.table::fread(lusc_surv_file)
  gcol <- colnames(ex)[1]
  mf <- data.table::fread(final_file)
  gl <- gsub("\\.\\d+$", "", mf$Gene_Symbol)
  for (g in unique(gl)){
    row <- ex[ex[[gcol]] == g, ]
    if (nrow(row) == 0) next
    vals <- as.numeric(row[, -1])
    sids <- colnames(row)[-1]
    df <- data.frame(sample = sids, expression = vals, stringsAsFactors = FALSE)
    df <- merge(df, su, by = "sample", all.x = FALSE)
    if (!all(c("OS.time", "OS") %in% colnames(df))) next
    cc <- complete.cases(df[, c("expression", "OS.time", "OS")])
    df <- df[cc, ]
    if (nrow(df) < 10) next
    med <- stats::median(df$expression, na.rm = TRUE)
    df$group <- ifelse(df$expression > med, "High", "Low")
    fit <- survival::survfit(survival::Surv(OS.time, OS == 1) ~ group, data = df)
    plt <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, palette = c("Low" = "#1f77b4", "High" = "#d62728"), legend.title = NULL, legend.labs = c(paste0(g, " Low"), paste0(g, " High")), ggtheme = ggplot2::theme_bw(base_size = 12), title = paste0(g, " expression vs OS (LUSC)"), xlab = "Time (months)", ylab = "Overall survival probability")
    cow <- suppressWarnings(cowplot::plot_grid(plt$plot, plt$table, ncol = 1, rel_heights = c(2, 0.7)))
    ggplot2::ggsave(file.path(plots_dir, "survival", paste0(g, "_survival_KM.png")), cow, width = 6, height = 4, dpi = 600, bg = "transparent")
    ggplot2::ggsave(file.path(plots_dir, "survival", paste0(g, "_survival_KM.pdf")), cow, width = 6, height = 4)
    ggplot2::ggsave(file.path(plots_dir, "survival", paste0(g, "_survival_KM.tiff")), cow, width = 6, height = 4, dpi = 600, compression = "lzw")
  }
}

if (successful == length(steps)) {
  cat("恭喜！所有步骤均成功完成。\n")
  cat("请查看最终主表文件获取完整的多组学数据。\n\n")
} else {
  cat("部分步骤执行失败，请查看上述详细结果。\n")
  cat("您可以单独运行失败的步骤进行调试。\n\n")
}

cat("下一步建议:\n")
cat("1. 检查最终主表文件的数据完整性\n")
cat("2. 根据需要进行数据可视化\n")
cat("3. 将结果整合到论文中\n")
cat("4. 进行功能验证实验\n\n")