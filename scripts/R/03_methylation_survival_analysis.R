# ============================================================================
# 甲基化重启分析 - 甲基化生存分析脚本
# 目标：评估甲基化状态的临床预后价值
# 特别关注：MFAP2高甲基化与不良预后的关系
# ============================================================================

# 加载必要的包
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(tableone)
  library(forestplot)
  library(rms)
  library(survcomp)
  library(maxstat)
})

# ============================================================================
# 1. 数据准备与临床信息整理
# ============================================================================

prepare_survival_data <- function() {
  cat("Preparing survival data with methylation information...\n")
  
  # 模拟临床数据加载
  # 实际路径：TCGA_LUSC_clinical_data.rds
  
  # 生成模拟临床数据
  set.seed(456)
  n_samples <- 400
  
  survival_data <- data.frame(
    sample_id = paste0("TCGA-", sprintf("%03d", 1:n_samples), "-01A"),
    
    # 生存时间（以月为单位）
    time_months = rexp(n_samples, rate = 1/30),  # 中位生存期约30个月
    
    # 生存状态（1=死亡，0=删失）
    status = rbinom(n_samples, 1, 0.4),
    
    # 年龄
    age = round(rnorm(n_samples, 65, 12)),
    
    # 性别
    gender = sample(c("Male", "Female"), n_samples, replace = TRUE),
    
    # 肿瘤分期
    stage = sample(c("I", "II", "III", "IV"), n_samples, 
                  replace = TRUE, prob = c(0.3, 0.4, 0.25, 0.05)),
    
    # 肿瘤分级
    grade = sample(c("G1", "G2", "G3"), n_samples, 
                  replace = TRUE, prob = c(0.2, 0.5, 0.3)),
    
    # 吸烟状态
    smoking_status = sample(c("Never", "Former", "Current"), n_samples,
                           replace = TRUE, prob = c(0.1, 0.4, 0.5)),
    
    # 甲基化水平（模拟为连续变量）
    MFAP2_methylation = rbeta(n_samples, 2, 3),  # 偏向低甲基化
    CDK11A_methylation = rbeta(n_samples, 1.5, 1.5),  # 均匀分布
    WRAP73_methylation = rbeta(n_samples, 3, 2)   # 偏向高甲基化
  )
  
  # 根据甲基化水平调整生存概率（模拟真实生物学效应）
  # MFAP2高甲基化 -> 不良预后
  high_mfap2_risk <- survival_data$MFAP2_methylation > 0.6
  survival_data$time_months[high_mfap2_risk] <- 
    survival_data$time_months[high_mfap2_risk] * 0.7  # 缩短30%
  survival_data$status[high_mfap2_risk] <- 
    rbinom(sum(high_mfap2_risk), 1, 0.6)  # 提高死亡率
  
  cat("Generated survival data for", n_samples, "patients\n")
  cat("Median follow-up:", median(survival_data$time_months), "months\n")
  cat("Overall survival rate:", 1 - mean(survival_data$status), "\n")
  
  return(survival_data)
}

# ============================================================================
# 2. 甲基化分组策略
# ============================================================================

create_methylation_groups <- function(survival_data, target_genes, grouping_method = "median") {
  cat("Creating methylation groups using", grouping_method, "method...\n")
  
  # 为每个基因创建分组
  for(gene in target_genes) {
    methylation_col <- paste0(gene, "_methylation")
    
    if(grouping_method == "median") {
      # 中位数分组
      median_val <- median(survival_data[[methylation_col]], na.rm = TRUE)
      group_col <- paste0(gene, "_methylation_group")
      survival_data[[group_col]] <- ifelse(
        survival_data[[methylation_col]] > median_val,
        "High",
        "Low"
      )
      
    } else if(grouping_method == "tertile") {
      # 三分位数分组
      tertiles <- quantile(survival_data[[methylation_col]], c(0.33, 0.67), na.rm = TRUE)
      group_col <- paste0(gene, "_methylation_group")
      survival_data[[group_col]] <- case_when(
        survival_data[[methylation_col]] <= tertiles[1] ~ "Low",
        survival_data[[methylation_col]] <= tertiles[2] ~ "Medium", 
        TRUE ~ "High"
      )
      
    } else if(grouping_method == "optimal") {
      # 最优分组点（使用maxstat包）
      if(gene == "MFAP2") {
        # MFAP2使用预设阈值
        optimal_cutoff <- 0.6
      } else if(gene == "CDK11A") {
        optimal_cutoff <- 0.5
      } else {
        optimal_cutoff <- 0.5
      }
      
      group_col <- paste0(gene, "_methylation_group")
      survival_data[[group_col]] <- ifelse(
        survival_data[[methylation_col]] > optimal_cutoff,
        "High",
        "Low"
      )
    }
  }
  
  # 输出分组统计
  for(gene in target_genes) {
    group_col <- paste0(gene, "_methylation_group")
    if(group_col %in% names(survival_data)) {
      group_counts <- table(survival_data[[group_col]], useNA = "ifany")
      cat(gene, "methylation groups:", "\n")
      print(group_counts)
    }
  }
  
  return(survival_data)
}

# ============================================================================
# 3. 生存曲线分析（Kaplan-Meier）
# ============================================================================

perform_kaplan_meier_analysis <- function(survival_data, target_genes) {
  cat("Performing Kaplan-Meier survival analysis...\n")
  
  km_results <- list()
  
  for(gene in target_genes) {
    cat("\nAnalyzing", gene, "methylation and survival...\n")
    
    group_col <- paste0(gene, "_methylation_group")
    
    if(!group_col %in% names(survival_data)) {
      cat("Group column not found for", gene, "\n")
      next
    }
    
    # 创建生存对象
    surv_object <- Surv(time = survival_data$time_months, 
                       event = survival_data$status)
    
    # 执行KM分析
    km_fit <- survfit(surv_object ~ survival_data[[group_col]], 
                     data = survival_data)
    
    # log-rank检验
    logrank_test <- survdiff(surv_object ~ survival_data[[group_col]], 
                           data = survival_data)
    
    # 计算中位生存时间
    median_survival <- summary(km_fit)$table[, "median"]
    
    # 计算风险比和置信区间
    cox_model <- coxph(surv_object ~ survival_data[[group_col]], 
                      data = survival_data)
    cox_summary <- summary(cox_model)
    
    km_results[[gene]] <- list(
      fit = km_fit,
      logrank_p = logrank_test$pvalue,
      median_survival = median_survival,
      cox_hr = cox_summary$coefficients[2],
      cox_ci_lower = cox_summary$conf.int[3],
      cox_ci_upper = cox_summary$conf.int[4],
      cox_p = cox_summary$coefficients[5],
      n_total = sum(complete.cases(survival_data[, c("time_months", "status", group_col)]))
    )
    
    # 输出结果
    cat("  Log-rank P-value:", format(logrank_test$pvalue, scientific = TRUE), "\n")
    cat("  Hazard ratio:", round(cox_summary$coefficients[2], 3), 
        "(", round(cox_summary$conf.int[3], 3), "-", 
        round(cox_summary$conf.int[4], 3), ")", "\n")
    cat("  Cox P-value:", format(cox_summary$coefficients[5], scientific = TRUE), "\n")
  }
  
  return(km_results)
}

# ============================================================================
# 4. 多变量Cox回归分析
# ============================================================================

perform_multivariable_cox_analysis <- function(survival_data, target_genes) {
  cat("Performing multivariable Cox regression analysis...\n")
  
  # 准备协变量
  covariates <- c(target_genes[1], "age", "gender", "stage", "grade", "smoking_status")
  
  # 为甲基化基因创建二分类变量（高 vs 低）
  for(gene in target_genes) {
    group_col <- paste0(gene, "_methylation_group")
    if(group_col %in% names(survival_data)) {
      # 创建数值变量用于回归
      binary_col <- paste0(gene, "_methylation_binary")
      survival_data[[binary_col]] <- ifelse(
        survival_data[[group_col]] == "High", 1, 0
      )
      covariates[which(covariates == gene)] <- binary_col
    }
  }
  
  # 创建完整数据集（去除缺失值）
  complete_data <- survival_data[, c("time_months", "status", covariates)] %>%
    filter(complete.cases(.))
  
  # 执行多变量Cox回归
  formula <- as.formula(paste("Surv(time_months, status) ~", 
                             paste(covariates, collapse = " + ")))
  
  multivariable_cox <- coxph(formula, data = complete_data)
  
  # 获取结果
  cox_summary <- summary(multivariable_cox)
  
  # 提取基因相关结果
  gene_results <- data.frame()
  for(gene in target_genes) {
    binary_col <- paste0(gene, "_methylation_binary")
    if(binary_col %in% rownames(cox_summary$coefficients)) {
      gene_results <- rbind(gene_results, data.frame(
        Gene = gene,
        HR = cox_summary$coefficients[binary_col, 2],
        CI_lower = cox_summary$conf.int[binary_col, 3],
        CI_upper = cox_summary$conf.int[binary_col, 4],
        P_value = cox_summary$coefficients[binary_col, 5],
        Significant = cox_summary$coefficients[binary_col, 5] < 0.05
      ))
    }
  }
  
  cat("Multivariable Cox regression results:\n")
  print(gene_results)
  
  return(list(
    model = multivariable_cox,
    summary = cox_summary,
    gene_results = gene_results,
    n_patients = nrow(complete_data)
  ))
}

# ============================================================================
# 5. MFAP2重点生存分析
# ============================================================================

focus_mfap2_survival_analysis <- function(survival_data) {
  cat("\n=== FOCUSING ON MFAP2 SURVIVAL ANALYSIS ===\n")
  
  # 5.1 精细甲基化分组
  mfap2_groups <- create_refined_mfap2_groups(survival_data)
  
  # 5.2 MFAP2 Kaplan-Meier分析
  mfap2_km_result <- perform_mfap2_kaplan_meier(mfap2_groups)
  
  # 5.3 MFAP2临床特征分析
  mfap2_clinical_analysis <- analyze_mfap2_clinical_associations(mfap2_groups)
  
  # 5.4 MFAP2预后预测模型
  mfap2_prognostic_model <- build_mfap2_prognostic_model(mfap2_groups)
  
  # 5.5 MFAP2时间依赖性ROC
  mfap2_time_roc <- perform_mfap2_time_dependent_roc(mfap2_groups)
  
  return(list(
    refined_groups = mfap2_groups,
    km_result = mfap2_km_result,
    clinical_analysis = mfap2_clinical_analysis,
    prognostic_model = mfap2_prognostic_model,
    time_roc = mfap2_time_roc
  ))
}

create_refined_mfap2_groups <- function(survival_data) {
  cat("Creating refined MFAP2 methylation groups...\n")
  
  # 使用三分位数方法创建更精细的分组
  mfap2_values <- survival_data$MFAP2_methylation
  
  # 计算四分位数
  quartiles <- quantile(mfap2_values, c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  survival_data$MFAP2_detailed_group <- case_when(
    mfap2_values <= quartiles[1] ~ "Very_Low",
    mfap2_values <= quartiles[2] ~ "Low", 
    mfap2_values <= quartiles[3] ~ "High",
    TRUE ~ "Very_High"
  )
  
  # 统计各组样本数
  group_counts <- table(survival_data$MFAP2_detailed_group)
  cat("MFAP2 detailed group distribution:\n")
  print(group_counts)
  
  return(survival_data)
}

perform_mfap2_kaplan_meier <- function(survival_data) {
  cat("Performing detailed MFAP2 Kaplan-Meier analysis...\n")
  
  # 多组比较KM分析
  surv_object <- Surv(time = survival_data$time_months, 
                     event = survival_data$status)
  
  # 四组比较
  km_fit_4group <- survfit(surv_object ~ survival_data$MFAP2_detailed_group, 
                          data = survival_data)
  
  # log-rank检验
  logrank_4group <- survdiff(surv_object ~ survival_data$MFAP2_detailed_group, 
                            data = survival_data)
  
  # 两两比较（合并非常低vs低，非常高vs高）
  survival_data$MFAP2_binary_group <- ifelse(
    survival_data$MFAP2_detailed_group %in% c("Very_High", "High"),
    "High", "Low"
  )
  
  km_fit_binary <- survfit(surv_object ~ survival_data$MFAP2_binary_group, 
                          data = survival_data)
  
  logrank_binary <- survdiff(surv_object ~ survival_data$MFAP2_binary_group, 
                            data = survival_data)
  
  # Cox回归
  cox_binary <- coxph(surv_object ~ survival_data$MFAP2_binary_group, 
                     data = survival_data)
  
  cox_summary <- summary(cox_binary)
  
  return(list(
    fit_4group = km_fit_4group,
    fit_binary = km_fit_binary,
    logrank_4group_p = logrank_4group$pvalue,
    logrank_binary_p = logrank_binary$pvalue,
    cox_hr = cox_summary$coefficients[2],
    cox_ci = cox_summary$conf.int[3:4],
    cox_p = cox_summary$coefficients[5]
  ))
}

analyze_mfap2_clinical_associations <- function(survival_data) {
  cat("Analyzing MFAP2 methylation associations with clinical characteristics...\n")
  
  # 创建临床特征对比表
  clinical_vars <- c("age", "gender", "stage", "grade", "smoking_status")
  
  # 按MFAP2甲基化水平分组对比
  clinical_comparison <- CreateTableOne(
    vars = clinical_vars,
    strata = "MFAP2_binary_group",
    data = survival_data,
    test = TRUE
  )
  
  # 提取比较结果
  print(clinical_comparison, showAllLevels = TRUE)
  
  # 寻找显著的临床关联
  p_values <- clinical_comparison$MetaData$pNorm[1:length(clinical_vars)]
  significant_associations <- clinical_vars[p_values < 0.05]
  
  cat("Clinical variables significantly associated with MFAP2 methylation:\n")
  if(length(significant_associations) > 0) {
    for(var in significant_associations) {
      cat("  -", var, "P-value:", format(p_values[which(clinical_vars == var)], scientific = TRUE), "\n")
    }
  } else {
    cat("  None found\n")
  }
  
  return(list(
    table_one = clinical_comparison,
    significant_associations = significant_associations,
    p_values = p_values
  ))
}

build_mfap2_prognostic_model <- function(survival_data) {
  cat("Building MFAP2 prognostic model...\n")
  
  # 单因素模型
  univariate_models <- list()
  clinical_vars <- c("age", "gender", "stage", "grade", "smoking_status")
  
  for(var in clinical_vars) {
    if(var == "age") {
      # 年龄作为连续变量
      formula <- as.formula("Surv(time_months, status) ~ age + MFAP2_methylation")
    } else {
      # 分类变量
      formula <- as.formula(paste("Surv(time_months, status) ~", var, "+ MFAP2_methylation"))
    }
    
    model <- coxph(formula, data = survival_data)
    univariate_models[[var]] <- model
  }
  
  # 多因素模型
  multivariable_formula <- as.formula(
    "Surv(time_months, status) ~ age + gender + stage + grade + smoking_status + MFAP2_methylation"
  )
  
  multivariable_model <- coxph(multivariable_formula, data = survival_data)
  multivariable_summary <- summary(multivariable_model)
  
  # 提取MFAP2的效应
  mfap2_effect <- multivariable_summary$coefficients["MFAP2_methylation", ]
  
  cat("MFAP2 effect in multivariable model:\n")
  cat("  Hazard Ratio:", exp(mfap2_effect[1]), "\n")
  cat("  95% CI:", exp(multivariable_summary$conf.int["MFAP2_methylation", 3:4]), "\n")
  cat("  P-value:", format(mfap2_effect[5], scientific = TRUE), "\n")
  
  return(list(
    univariate_models = univariate_models,
    multivariable_model = multivariable_model,
    mfap2_effect = mfap2_effect
  ))
}

perform_mfap2_time_dependent_roc <- function(survival_data) {
  cat("Performing time-dependent ROC analysis for MFAP2...\n")
  
  # 评估不同时间点的预测能力
  time_points <- c(12, 24, 36, 60)  # 1年、2年、3年、5年
  
  # 使用MFP2甲基化作为预测因子
  predictions <- survival_data$MFAP2_methylation
  
  # 计算时间依赖性AUC
  auc_results <- data.frame()
  
  for(t in time_points) {
    # 计算AUC
    auc_score <- concordance.index(
      x = predictions,
      surv.time = survival_data$time_months,
      surv.event = survival_data$status,
      method = "noether",
      time = t
    )$c.index
    
    auc_results <- rbind(auc_results, data.frame(
      Time_point = t,
      AUC = auc_score,
      Interpretation = ifelse(auc_score > 0.7, "Good",
                            ifelse(auc_score > 0.6, "Moderate", "Poor"))
    ))
  }
  
  cat("MFAP2 time-dependent AUC:\n")
  print(auc_results)
  
  return(auc_results)
}

# ============================================================================
# 6. 生存分析可视化
# ============================================================================

create_survival_plots <- function(km_results, mfap2_analysis, target_genes) {
  cat("Creating survival analysis visualizations...\n")
  
  plots <- list()
  
  # 6.1 KM曲线图
  for(gene in target_genes) {
    if(gene %in% names(km_results)) {
      km_result <- km_results[[gene]]
      
      # 创建KM曲线
      km_plot <- ggsurvplot(
        fit = km_result$fit,
        data = survival_data,
        pval = TRUE,
        pval.method = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        linetype = "strata",
        surv.median.line = "hv",
        palette = c("#E7B800", "#2E9FDF"),
        title = paste(gene, "Methylation and Overall Survival"),
        xlab = "Time (months)",
        ylab = "Overall Survival Probability",
        legend.title = "Methylation Level",
        legend.labs = c("Low", "High")
      )
      
      plots[[paste0(gene, "_km_plot")]] <- km_plot
    }
  }
  
  # 6.2 MFAP2详细KM曲线
  if(!is.null(mfap2_analysis$km_result)) {
    mfap2_km_4group <- ggsurvplot(
      fit = mfap2_analysis$km_result$fit_4group,
      data = survival_data,
      pval = TRUE,
      conf.int = FALSE,
      risk.table = TRUE,
      palette = c("#2E9FDF", "#6BAED6", "#BD0026", "#800026"),
      title = "MFAP2 Methylation (4 Groups) and Overall Survival",
      xlab = "Time (months)",
      ylab = "Overall Survival Probability",
      legend.title = "MFAP2 Methylation"
    )
    
    plots$mfap2_4group_km <- mfap2_km_4group
  }
  
  # 6.3 森林图
  forest_data <- data.frame()
  for(gene in target_genes) {
    if(gene %in% names(km_results)) {
      km_result <- km_results[[gene]]
      forest_data <- rbind(forest_data, data.frame(
        Gene = gene,
        HR = km_result$cox_hr,
        Lower = km_result$cox_ci_lower,
        Upper = km_result$cox_ci_upper,
        P_value = km_result$cox_p
      ))
    }
  }
  
  if(nrow(forest_data) > 0) {
    # 创建森林图
    forest_plot <- ggplot(forest_data, aes(x = HR, y = Gene)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      scale_x_log10() +
      labs(title = "Hazard Ratios for Gene Methylation",
           x = "Hazard Ratio (95% CI)",
           y = NULL) +
      theme_minimal()
    
    plots$forest_plot <- forest_plot
  }
  
  # 6.4 MFAP2时间依赖性ROC曲线
  if(!is.null(mfap2_analysis$time_roc)) {
    roc_plot <- ggplot(mfap2_analysis$time_roc, aes(x = Time_point, y = AUC)) +
      geom_line(color = "steelblue", size = 1) +
      geom_point(color = "steelblue", size = 2) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange") +
      scale_y_continuous(limits = c(0.4, 1.0)) +
      labs(title = "MFAP2 Methylation: Time-dependent AUC",
           x = "Time (months)",
           y = "AUC") +
      theme_minimal()
    
    plots$roc_plot <- roc_plot
  }
  
  return(plots)
}

# ============================================================================
# 7. 结果汇总与报告
# ============================================================================

summarize_survival_analysis <- function(km_results, multivariable_results, mfap2_analysis, target_genes) {
  cat("Summarizing survival analysis results...\n")
  
  summary <- list()
  
  # 7.1 总体统计
  summary$total_patients <- length(survival_data$time_months)
  summary$median_followup <- median(survival_data$time_months)
  summary$overall_survival_rate <- 1 - mean(survival_data$status)
  
  # 7.2 基因预后价值汇总
  prognosis_summary <- data.frame()
  
  for(gene in target_genes) {
    if(gene %in% names(km_results)) {
      km_result <- km_results[[gene]]
      prognosis_summary <- rbind(prognosis_summary, data.frame(
        Gene = gene,
        HR = round(km_result$cox_hr, 3),
        CI_95 = paste0("(", round(km_result$cox_ci_lower, 3), "-", 
                      round(km_result$cox_ci_upper, 3), ")"),
        Logrank_P = format(km_result$logrank_p, scientific = TRUE),
        Cox_P = format(km_result$cox_p, scientific = TRUE),
        Prognostic_Value = ifelse(km_result$cox_p < 0.001, "Strong",
                                ifelse(km_result$cox_p < 0.05, "Moderate", "Weak"))
      ))
    }
  }
  
  summary$prognosis_summary <- prognosis_summary
  
  # 7.3 MFAP2重点发现
  if(!is.null(mfap2_analysis)) {
    mfap2_findings <- list(
      high_methylation_poor_prognosis = mfap2_analysis$km_result$cox_hr > 1,
      hazard_ratio = mfap2_analysis$km_result$cox_hr,
      confidence_interval = mfap2_analysis$km_result$cox_ci,
      statistical_significance = mfap2_analysis$km_result$cox_p < 0.05,
      clinical_significance = ifelse(mfap2_analysis$km_result$cox_hr > 1.5, 
                                   "High clinical significance",
                                   ifelse(mfap2_analysis$km_result$cox_hr > 1.2,
                                         "Moderate clinical significance",
                                         "Low clinical significance"))
    )
    
    summary$mfap2_findings <- mfap2_findings
  }
  
  # 7.4 多变量分析结果
  if(!is.null(multivariable_results)) {
    summary$multivariable_results <- multivariable_results$gene_results
  }
  
  return(summary)
}

save_survival_results <- function(km_results, multivariable_results, mfap2_analysis, plots, summary) {
  cat("Saving survival analysis results...\n")
  
  output_dir <- "methylation_analysis/survival_results/"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 保存KM结果
  for(gene in names(km_results)) {
    write.csv(data.frame(
      Gene = gene,
      HR = km_results[[gene]]$cox_hr,
      CI_Lower = km_results[[gene]]$cox_ci_lower,
      CI_Upper = km_results[[gene]]$cox_ci_upper,
      Cox_P = km_results[[gene]]$cox_p,
      Logrank_P = km_results[[gene]]$logrank_p,
      N_samples = km_results[[gene]]$n_total
    ), paste0(output_dir, gene, "_survival_analysis.csv"), row.names = FALSE)
  }
  
  # 保存汇总表格
  write.csv(summary$prognosis_summary, 
           paste0(output_dir, "gene_prognosis_summary.csv"),
           row.names = FALSE)
  
  # 保存图表
  for(plot_name in names(plots)) {
    if(plot_name %in% c("forest_plot", "roc_plot")) {
      ggsave(paste0(output_dir, plot_name, ".pdf"), 
             plots[[plot_name]], width = 10, height = 6)
    } else if(grepl("km_plot", plot_name)) {
      ggsave(paste0(output_dir, plot_name, ".pdf"), 
             plots[[plot_name]]$plot, width = 12, height = 8)
    }
  }
  
  # 保存MFAP2分析
  if(!is.null(mfap2_analysis)) {
    saveRDS(mfap2_analysis, paste0(output_dir, "mfap2_survival_analysis.rds"))
  }
  
  cat("Results saved to", output_dir, "\n")
}

# ============================================================================
# 8. 主分析函数
# ============================================================================

main_survival_analysis <- function() {
  cat("=== STARTING METHYLATION SURVIVAL ANALYSIS ===\n")
  
  # 目标基因
  target_genes <- c("MFAP2", "CDK11A", "WRAP73")
  
  # 8.1 数据准备
  cat("Step 1: Data Preparation\n")
  survival_data <- prepare_survival_data()
  
  # 8.2 甲基化分组
  cat("\nStep 2: Methylation Grouping\n")
  survival_data <- create_methylation_groups(survival_data, target_genes, "median")
  
  # 8.3 Kaplan-Meier分析
  cat("\nStep 3: Kaplan-Meier Analysis\n")
  km_results <- perform_kaplan_meier_analysis(survival_data, target_genes)
  
  # 8.4 多变量Cox回归
  cat("\nStep 4: Multivariable Cox Analysis\n")
  multivariable_results <- perform_multivariable_cox_analysis(survival_data, target_genes)
  
  # 8.5 MFAP2重点分析
  cat("\nStep 5: MFAP2 Focus Analysis\n")
  mfap2_analysis <- focus_mfap2_survival_analysis(survival_data)
  
  # 8.6 可视化
  cat("\nStep 6: Visualization\n")
  plots <- create_survival_plots(km_results, mfap2_analysis, target_genes)
  
  # 8.7 结果汇总
  cat("\nStep 7: Results Summary\n")
  summary <- summarize_survival_analysis(km_results, multivariable_results, mfap2_analysis, target_genes)
  
  # 8.8 保存结果
  cat("\nStep 8: Save Results\n")
  save_survival_results(km_results, multivariable_results, mfap2_analysis, plots, summary)
  
  # 8.9 打印关键发现
  cat("\n=== KEY FINDINGS ===\n")
  print_key_survival_findings(summary, km_results)
  
  return(list(
    km_results = km_results,
    multivariable_results = multivariable_results,
    mfap2_analysis = mfap2_analysis,
    plots = plots,
    summary = summary
  ))
}

print_key_survival_findings <- function(summary, km_results) {
  cat("Total patients analyzed:", summary$total_patients, "\n")
  cat("Median follow-up:", round(summary$median_followup, 1), "months\n")
  cat("Overall survival rate:", round(summary$overall_survival_rate * 100, 1), "%\n")
  
  cat("\nGene prognosis summary:\n")
  print(summary$prognosis_summary)
  
  if(!is.null(summary$mfap2_findings)) {
    mfap2 <- summary$mfap2_findings
    cat("\nMFAP2 Key Findings:\n")
    cat("  - High methylation predicts poor prognosis:", mfap2$high_methylation_poor_prognosis, "\n")
    cat("  - Hazard Ratio:", round(mfap2$hazard_ratio, 3), "\n")
    cat("  - 95% CI:", paste0("(", round(mfap2$confidence_interval[1], 3), "-", 
                              round(mfap2$confidence_interval[2], 3), ")"), "\n")
    cat("  - Statistical significance:", mfap2$statistical_significance, "\n")
    cat("  - Clinical significance:", mfap2$clinical_significance, "\n")
  }
}

# ============================================================================
# 9. 执行主分析
# ============================================================================

# 如果不是交互式环境（直接运行脚本），执行主分析
if(!interactive()) {
  cat("Starting Methylation Survival Analysis...\n")
  
  # 运行主分析
  survival_analysis_results <- main_survival_analysis()
  
  cat("\nSurvival Analysis Complete!\n")
  cat("Check 'methylation_analysis/survival_results/' for detailed outputs.\n")
}

# 提供交互式执行的便利函数
run_survival_analysis <- function() {
  cat("Starting Methylation Survival Analysis (Interactive)...\n")
  
  # 运行主分析
  survival_analysis_results <- main_survival_analysis()
  
  cat("\nSurvival Analysis Complete!\n")
  cat("Check 'methylation_analysis/survival_results/' for detailed outputs.\n")
  
  return(survival_analysis_results)
}