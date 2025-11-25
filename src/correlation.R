suppressMessages(library(data.table))
map_symbol_to_ensembl <- function(genes) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
  m <- try(AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("ENSEMBL")), silent = TRUE)
  if (inherits(m,"try-error") || is.null(m)) return(NULL)
  as.data.table(m)[!is.na(ENSEMBL), .(SYMBOL, ENSEMBL)]
}
corr_combined <- function(best, m_all, rl, rs, map, abs_thr, p_thr) {
  out <- list()
  for (g in unique(best$Gene)) {
    ens <- map[SYMBOL == g]$ENSEMBL[1]
    bp <- best[Gene == g]$Probe[1]
    gb <- m_all[Probe == bp]
    ra <- rbindlist(list(rl[ENSEMBL == ens], rs[ENSEMBL == ens]), use.names = TRUE, fill = TRUE)
    dt <- merge(gb[, .(pid, Beta)], ra[, .(pid, TPM)], by = "pid")
    if (nrow(dt) >= 10) {
      r <- suppressWarnings(cor(dt$Beta, dt$TPM, method = "spearman", use = "complete.obs"))
      p <- tryCatch(suppressWarnings(cor.test(dt$Beta, dt$TPM, method = "spearman")$p.value), error = function(e) NA_real_)
      out[[length(out) + 1]] <- data.table(Gene = g, Spearman = r, Corr_P = p)
    }
  }
  if (length(out) > 0) rbindlist(out) else NULL
}

corr_by_cohort <- function(bestC, m_all, rl, rs, map) {
  out <- list()
  genes <- unique(bestC$Gene)
  for (g in genes) {
    ens <- map[SYMBOL == g]$ENSEMBL[1]
    for (coh in c("LUAD", "LUSC")) {
      bp <- bestC[Cohort == coh & Gene == g]$Probe[1]
      gb <- m_all[Cohort == coh & Probe == bp]
      rll <- if (coh == "LUAD") rl else rs
      dt <- merge(gb[, .(pid, Beta)], rll[ENSEMBL == ens, .(pid, TPM)], by = "pid")
      if (nrow(dt) >= 10) {
        r <- suppressWarnings(cor(dt$Beta, dt$TPM, method = "spearman", use = "complete.obs"))
        p <- tryCatch(suppressWarnings(cor.test(dt$Beta, dt$TPM, method = "spearman")$p.value), error = function(e) NA_real_)
        out[[length(out) + 1]] <- data.table(Cohort = coh, Gene = g, Spearman = r, Corr_P = p)
      }
    }
  }
  if (length(out) > 0) rbindlist(out) else NULL
}
