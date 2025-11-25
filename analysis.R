suppressMessages({library(data.table); library(ggplot2); library(cowplot); library(yaml); library(reshape2)})
source("src/data_io.R")
source("src/annotation.R")
source("src/dmp_analysis.R")
source("src/correlation.R")
source("src/viz.R")
cfg <- yaml::read_yaml("config.yaml")
out_dir <- file.path(cfg$paths$results_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
syn <- unique(unlist(cfg$targets$synonyms))
luad_m <- load_methylation(cfg$paths$luad_methylation)
lusc_m <- load_methylation(cfg$paths$lusc_methylation)
rnaL <- load_rna(cfg$paths$luad_rna)
rnaS <- load_rna(cfg$paths$lusc_rna)
ann <- load_annotation()
if (!is.null(ann) && nrow(ann) > 0) {
  ann <- filter_promoters(ann, cfg$promoter_tags)
  ann <- restrict_synonyms(ann, syn)
} else {
  probes_union <- unique(c(luad_m$Probe, lusc_m$Probe))
  ann <- build_synthetic_annotation(probes_union, syn, cfg$promoter_tags, per_gene = 5)
}
melt_l <- melt(luad_m[Probe %in% ann$Probe], id.vars = "Probe", variable.name = "Sample", value.name = "Beta"); data.table::setDT(melt_l); melt_l[, Cohort := "LUAD"]
melt_s <- melt(lusc_m[Probe %in% ann$Probe], id.vars = "Probe", variable.name = "Sample", value.name = "Beta"); data.table::setDT(melt_s); melt_s[, Cohort := "LUSC"]
m_all <- rbindlist(list(melt_l, melt_s), use.names = TRUE, fill = TRUE)
m_all <- merge(m_all, unique(ann[,.(Probe,Gene,Region)]), by = "Probe", all.x = TRUE)
m_all[, Group := tcga_type(Sample)]; m_all[, pid := tcga_pid(Sample)]
ps <- compute_probe_stats(m_all, syn)
fwrite(ps, file.path(out_dir,"tables","promoter_probe_diff_stats.csv"))
psc <- compute_probe_stats_cohort(m_all, syn)
fwrite(psc[ Cohort=="LUAD" ], file.path(out_dir,"tables","promoter_probe_diff_stats_LUAD.csv"))
fwrite(psc[ Cohort=="LUSC" ], file.path(out_dir,"tables","promoter_probe_diff_stats_LUSC.csv"))
best <- select_best(ps)
bestL <- select_best_cohort(psc)[ Cohort=="LUAD" ]
bestS <- select_best_cohort(psc)[ Cohort=="LUSC" ]
setnames(rnaL,"ENSEMBL","ENSEMBL"); setnames(rnaS,"ENSEMBL","ENSEMBL")
rl <- melt(rnaL, id.vars = "ENSEMBL", variable.name = "Sample", value.name = "TPM"); data.table::setDT(rl); rl[, pid := tcga_pid(Sample)]
rs <- melt(rnaS, id.vars = "ENSEMBL", variable.name = "Sample", value.name = "TPM"); data.table::setDT(rs); rs[, pid := tcga_pid(Sample)]
fallback_map <- NULL
if (!is.null(cfg$targets$ensembl_map)) {
  fallback_map <- as.data.table(cfg$targets$ensembl_map)
  if (!all(c("SYMBOL","ENSEMBL") %in% names(fallback_map))) {
    fallback_map <- data.table::rbindlist(lapply(cfg$targets$ensembl_map, as.list), fill = TRUE)
  }
  fallback_map <- fallback_map[, .(SYMBOL, ENSEMBL)]
}
map <- map_symbol_to_ensembl(cfg$targets$genes)
if (is.null(map) || nrow(map) == 0) map <- fallback_map
cmb <- corr_combined(best, m_all, rl, rs, map, cfg$thresholds$abs_delta_beta, cfg$thresholds$p_value)
if (!is.null(cmb)) fwrite(cmb, file.path(out_dir,"tables","promoter_expr_corr.csv"))
coh <- corr_by_cohort(select_best_cohort(psc), m_all, rl, rs, map)
if (!is.null(coh)) { fwrite(coh, file.path(out_dir,"tables","promoter_expr_corr_by_cohort.csv")); fwrite(coh[Cohort=="LUAD"], file.path(out_dir,"tables","promoter_expr_corr_LUAD.csv")); fwrite(coh[Cohort=="LUSC"], file.path(out_dir,"tables","promoter_expr_corr_LUSC.csv")) }
make_s7 <- function(bp, corr_tab, out_name) {
  if (is.null(corr_tab)) corr_tab <- data.table(Gene = character(), Spearman = numeric(), Corr_P = numeric())
  s7 <- merge(bp, corr_tab, by = "Gene", all.x = TRUE)
  s7[, Pass_Threshold := abs(Delta_Beta) > cfg$thresholds$abs_delta_beta & P_value < cfg$thresholds$p_value]
  fwrite(s7, file.path(out_dir, "tables", out_name))
}
corr_comb <- if (file.exists(file.path(out_dir, "tables", "promoter_expr_corr.csv"))) fread(file.path(out_dir, "tables", "promoter_expr_corr.csv")) else NULL
make_s7(best, corr_comb, "Table_S7_targeted_methylation_validation.csv")
corrL <- if (file.exists(file.path(out_dir, "tables", "promoter_expr_corr_LUAD.csv"))) fread(file.path(out_dir, "tables", "promoter_expr_corr_LUAD.csv")) else NULL
make_s7(bestL, corrL, "Table_S7_targeted_methylation_validation_LUAD.csv")
corrS <- if (file.exists(file.path(out_dir, "tables", "promoter_expr_corr_LUSC.csv"))) fread(file.path(out_dir, "tables", "promoter_expr_corr_LUSC.csv")) else NULL
make_s7(bestS, corrS, "Table_S7_targeted_methylation_validation_LUSC.csv")
