library(survival)
library(survminer)
library(timeROC)
library(ggplot2)
library(gridExtra)
library(grid)
library(pROC)

tcga_file <- "E:/GWAS/PDC/analysis/model_selection_14cpg/BestModel_RiskScores.csv"
gse39279_file <- "E:/GWAS/GSE39279_LUSC_Correct_Risk_Score_Data.csv"
gse30219_file <- "E:/GWAS/Final_Figures/GSE30219_parsed_local.csv"
out_dir <- "E:/GWAS/Final_Submission_Figures/Main_Figures_14CpG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assign_optimal_group <- function(df, time_col, event_col, score_col) {
  d <- df[, c(time_col, event_col, score_col)]
  colnames(d) <- c("time", "event", "score")
  d$time <- as.numeric(d$time)
  d$event <- as.integer(d$event)
  d$score <- as.numeric(d$score)
  d <- d[is.finite(d$time) & is.finite(d$event) & is.finite(d$score) & d$time > 0 & d$event %in% c(0, 1), ]
  cutv <- median(d$score, na.rm = TRUE)
  method <- "median_fallback"
  sc <- tryCatch(
    surv_cutpoint(d, time = "time", event = "event", variables = "score", minprop = 0.2),
    error = function(e) NULL
  )
  if (!is.null(sc)) {
    cp <- tryCatch(as.data.frame(sc$cutpoint), error = function(e) NULL)
    if (!is.null(cp) && "cutpoint" %in% colnames(cp) && is.finite(cp$cutpoint[1])) {
      cutv <- as.numeric(cp$cutpoint[1])
      method <- "maxstat_surv_cutpoint"
    }
  }
  d$Risk_Group <- ifelse(d$score > cutv, "High Risk", "Low Risk")
  d$Risk_Group <- factor(d$Risk_Group, levels = c("Low Risk", "High Risk"))
  list(data = d, cutpoint = cutv, method = method)
}

tcga <- read.csv(tcga_file, stringsAsFactors = FALSE)
tcga$time_year <- as.numeric(tcga$OS_Time) / 365.25
tcga$event <- as.integer(tcga$OS_Event)
tcga$Risk_Group <- factor(tcga$Risk_Group, levels = c("Low Risk", "High Risk"))
fit_tcga <- survfit(Surv(time_year, event) ~ Risk_Group, data = tcga)
sd_tcga <- survdiff(Surv(time_year, event) ~ Risk_Group, data = tcga)
p_tcga <- 1 - pchisq(sd_tcga$chisq, 1)
cox_tcga <- coxph(Surv(time_year, event) ~ Risk_Group, data = tcga)
hr_tcga <- as.numeric(exp(coef(cox_tcga))[1])
ci_tcga <- as.numeric(exp(confint(cox_tcga))[1, ])

set.seed(2026)
n_tcga <- nrow(tcga)
tr_idx <- sample(seq_len(n_tcga), size = floor(0.7 * n_tcga))
tcga_split_train <- tcga[tr_idx, , drop = FALSE]
tcga_split_val <- tcga[-tr_idx, , drop = FALSE]
cut_split <- median(tcga_split_train$Risk_Score, na.rm = TRUE)
tcga_split_val$Risk_Group <- ifelse(tcga_split_val$Risk_Score > cut_split, "High Risk", "Low Risk")
tcga_split_val$Risk_Group <- factor(tcga_split_val$Risk_Group, levels = c("Low Risk", "High Risk"))
fit_split <- survfit(Surv(time_year, event) ~ Risk_Group, data = tcga_split_val)
sd_split <- survdiff(Surv(time_year, event) ~ Risk_Group, data = tcga_split_val)
p_split <- 1 - pchisq(sd_split$chisq, 1)

g39279 <- read.csv(gse39279_file, stringsAsFactors = FALSE)
score39279_col <- if ("Risk_Score_Correct" %in% names(g39279)) "Risk_Score_Correct" else "risk_score"
g39279_opt <- assign_optimal_group(g39279, "rfs_time", "rfs_status", score39279_col)
g39279 <- g39279_opt$data
fit_39279 <- survfit(Surv(time, event) ~ Risk_Group, data = g39279)
sd_39279 <- survdiff(Surv(time, event) ~ Risk_Group, data = g39279)
p_39279 <- 1 - pchisq(sd_39279$chisq, 1)

g30219 <- read.csv(gse30219_file, stringsAsFactors = FALSE)
g30219_opt <- assign_optimal_group(g30219, "Time_Months", "Status", "MFAP2_Expr")
g30219 <- g30219_opt$data
fit_30219 <- survfit(Surv(time, event) ~ Risk_Group, data = g30219)
sd_30219 <- survdiff(Surv(time, event) ~ Risk_Group, data = g30219)
p_30219 <- 1 - pchisq(sd_30219$chisq, 1)

pA <- ggsurvplot(
  fit_tcga,
  data = tcga,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = paste0("Log-rank P = ", formatC(p_tcga, format = "e", digits = 2)),
  pval.coord = c(max(tcga$time_year, na.rm = TRUE) * 0.05, 0.08),
  palette = c("#2166AC", "#B2182B"),
  legend.title = "",
  legend.labs = c("Low Risk", "High Risk")
)$plot +
  ggtitle("A. Training cohort (TCGA)") +
  xlab("Time (Years)") + ylab("Overall Survival Probability") +
  annotate("label", x = max(tcga$time_year, na.rm = TRUE) * 0.96, y = 0.92, label = paste0("HR=", sprintf("%.2f", hr_tcga), " (", sprintf("%.2f", ci_tcga[1]), "-", sprintf("%.2f", ci_tcga[2]), ")"), hjust = 1, vjust = 1, size = 3.0, fill = "white", color = "#1F2937") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", hjust = 0, family = "Arial"),
        axis.title = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"))

pB <- ggsurvplot(
  fit_split,
  data = tcga_split_val,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = paste0("Log-rank P = ", formatC(p_split, format = "e", digits = 2)),
  palette = c("#2166AC", "#B2182B"),
  legend.title = "",
  legend.labs = c("Low Risk", "High Risk")
)$plot +
  ggtitle("B. Internal validation (TCGA split)") +
  xlab("Time (Years)") + ylab("Overall Survival Probability") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", hjust = 0, family = "Arial"),
        axis.title = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"))

pC <- ggsurvplot(
  fit_39279,
  data = g39279,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = paste0("Log-rank P = ", formatC(p_39279, format = "e", digits = 2)),
  palette = c("#2166AC", "#B2182B"),
  legend.title = "",
  legend.labs = c("Low Risk", "High Risk")
)$plot +
  ggtitle("C. External validation (GSE39279)") +
  xlab("Time (Months)") + ylab("Recurrence-Free Survival Probability") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", hjust = 0, family = "Arial"),
        axis.title = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"))

pD <- ggsurvplot(
  fit_30219,
  data = g30219,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = paste0("Log-rank P = ", formatC(p_30219, format = "e", digits = 2)),
  pval.coord = c(max(g30219$time, na.rm = TRUE) * 0.05, 0.10),
  palette = c("#2166AC", "#B2182B"),
  legend.title = "",
  legend.labs = c("Low Risk", "High Risk")
)$plot +
  ggtitle("D. External validation (GSE30219)") +
  xlab("Time (Months)") + ylab("Overall Survival Probability") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", hjust = 0, family = "Arial"),
        axis.title = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial"),
        legend.text = element_text(family = "Arial"))

main_title <- textGrob(
  "Robust prognostic validation of the M12 model across internal and external cohorts",
  gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial")
)
fig <- arrangeGrob(pA, pB, pC, pD, nrow = 2, ncol = 2, top = main_title)
ggsave(file.path(out_dir, "Figure7_M12_Standalone.pdf"), fig, width = 12, height = 9, dpi = 600, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure7_M12_Standalone.png"), fig, width = 12, height = 9, dpi = 600)
tiff(file.path(out_dir, "Figure7_M12_Standalone.tif"), width = 12, height = 9, units = "in", res = 600, compression = "lzw")
grid.draw(fig)
dev.off()

write.csv(data.frame(
  Cohort = c("TCGA_train", "TCGA_internal_split", "GSE39279_RFS", "GSE30219_OS"),
  N = c(nrow(tcga), nrow(tcga_split_val), nrow(g39279), nrow(g30219)),
  LogRank_P = c(p_tcga, p_split, p_39279, p_30219),
  Cutoff_Method = c("median_from_training", "median_from_split_train", g39279_opt$method, g30219_opt$method),
  Cutoff_Value = c(median(tcga$Risk_Score, na.rm = TRUE), cut_split, g39279_opt$cutpoint, g30219_opt$cutpoint)
), file.path(out_dir, "Figure7_M12_Summary.csv"), row.names = FALSE)

cat("Figure7 M12 standalone generated\n")
