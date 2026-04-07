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

roc_tcga <- timeROC(T = tcga$time_year, delta = tcga$event, marker = tcga$Risk_Score, cause = 1, weighting = "marginal", times = c(1, 3, 5), iid = TRUE)

tcga_cal <- tcga[, c("Risk_Score", "time_year", "event")]
tcga_cal <- tcga_cal[is.finite(tcga_cal$Risk_Score) & is.finite(tcga_cal$time_year) & is.finite(tcga_cal$event), ]
tcga_cal$risk01 <- (tcga_cal$Risk_Score - min(tcga_cal$Risk_Score, na.rm = TRUE)) / (max(tcga_cal$Risk_Score, na.rm = TRUE) - min(tcga_cal$Risk_Score, na.rm = TRUE) + 1e-12)
tcga_cal$y3 <- as.integer(tcga_cal$event == 1 & tcga_cal$time_year <= 3)
tcga_cal$bin <- cut(rank(tcga_cal$risk01, ties.method = "random"), breaks = 8, labels = FALSE)
cal_pred <- aggregate(risk01 ~ bin, data = tcga_cal, FUN = mean)
cal_obs <- aggregate(y3 ~ bin, data = tcga_cal, FUN = mean)
cal_n <- aggregate(y3 ~ bin, data = tcga_cal, FUN = length)
cal_df <- merge(merge(cal_pred, cal_obs, by = "bin"), cal_n, by = "bin")
colnames(cal_df) <- c("bin", "pred", "obs", "n")
se_cal <- sqrt(pmax(cal_df$obs * (1 - cal_df$obs) / pmax(cal_df$n, 1), 0))
cal_df$obs_l <- pmax(cal_df$obs - 1.96 * se_cal, 0)
cal_df$obs_u <- pmin(cal_df$obs + 1.96 * se_cal, 1)
cal_df <- cal_df[order(cal_df$pred), ]
cal_lm <- lm(obs ~ pred, data = cal_df)
cal_slope <- unname(coef(cal_lm)[2])
cal_r2 <- summary(cal_lm)$r.squared

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
  ggtitle("A. Survival stratification in the training cohort") +
  xlab("Time (Years)") + ylab("Overall Survival Probability") +
  annotate("label", x = max(tcga$time_year, na.rm = TRUE) * 0.96, y = 0.92, label = paste0("HR=", sprintf("%.2f", hr_tcga), " (", sprintf("%.2f", ci_tcga[1]), "-", sprintf("%.2f", ci_tcga[2]), ")"), hjust = 1, vjust = 1, size = 3.0, fill = "white", color = "#1F2937") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

pB <- ggsurvplot(
  fit_39279,
  data = g39279,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = paste0("Log-rank P = ", formatC(p_39279, format = "e", digits = 2)),
  palette = c("#2166AC", "#B2182B"),
  legend.title = "",
  legend.labs = c("Low Risk", "High Risk")
)$plot +
  ggtitle("B. External validation in GSE39279 (recurrence-free survival)") +
  xlab("Time (Months)") + ylab("Recurrence-Free Survival Probability") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

pC <- ggsurvplot(
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
  ggtitle("C. External validation in GSE30219 (overall survival)") +
  xlab("Time (Months)") + ylab("Overall Survival Probability") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

pD <- ggplot(cal_df, aes(x = pred, y = obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray60") +
  geom_ribbon(aes(ymin = obs_l, ymax = obs_u), fill = "#CBD5E1", alpha = 0.35, inherit.aes = TRUE) +
  geom_point(size = 2.5, color = "#2166AC") +
  geom_line(linewidth = 1.1, color = "#2166AC") +
  annotate("label", x = 0.04, y = 0.96, hjust = 0, vjust = 1, size = 3.0, fill = "white", color = "#1F2937",
           label = paste0("Slope=", sprintf("%.2f", cal_slope), ", R²=", sprintf("%.2f", cal_r2))) +
  labs(title = "D. Calibration (3-year)", x = "Predicted risk", y = "Observed event rate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

main_title <- textGrob(
  "Clinical validation of the M12 CpG risk model across independent cohorts",
  gp = gpar(fontsize = 15, fontface = "bold")
)
fig <- arrangeGrob(pA, pB, pC, pD, nrow = 2, ncol = 2, top = main_title)
ggsave(file.path(out_dir, "Figure7_M12_Standalone.pdf"), fig, width = 12, height = 9, dpi = 600, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure7_M12_Standalone.png"), fig, width = 12, height = 9, dpi = 600)
tiff(file.path(out_dir, "Figure7_M12_Standalone.tif"), width = 12, height = 9, units = "in", res = 600, compression = "lzw")
grid.draw(fig)
dev.off()

write.csv(data.frame(
  Cohort = c("TCGA_train", "GSE39279_RFS", "GSE30219_OS"),
  N = c(nrow(tcga), nrow(g39279), nrow(g30219)),
  LogRank_P = c(p_tcga, p_39279, p_30219),
  Cutoff_Method = c("median_from_training", g39279_opt$method, g30219_opt$method),
  Cutoff_Value = c(median(tcga$Risk_Score, na.rm = TRUE), g39279_opt$cutpoint, g30219_opt$cutpoint)
), file.path(out_dir, "Figure7_M12_Summary.csv"), row.names = FALSE)

cat("Figure7 M12 standalone generated\n")
