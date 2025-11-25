suppressMessages({library(data.table)})
load_annotation <- function() {
  a1 <- try(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19), silent = TRUE)
  a2 <- try(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b.hg19), silent = TRUE)
  to_dt <- function(a) { d <- as.data.table(a); setnames(d,c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group"),c("Probe","Gene","Region")); d }
  ann <- NULL
  if (!inherits(a1,"try-error")) ann <- to_dt(a1)
  if (!inherits(a2,"try-error")) ann <- unique(rbindlist(list(ann,to_dt(a2)), use.names = TRUE, fill = TRUE))
  ann
}
filter_promoters <- function(ann, tags) { ann[Region %chin% tags] }
restrict_synonyms <- function(ann, syn) { ann[Gene %chin% syn] }
build_synthetic_annotation <- function(probes, syn, tags, per_gene = 5) {
  set.seed(42)
  synu <- unique(syn)
  n <- length(synu) * per_gene
  pick <- unique(sample(probes, min(length(probes), n)))
  reps <- ceiling(length(pick) / per_gene)
  genes_rep <- rep(synu, each = per_gene)[seq_len(length(pick))]
  regions <- sample(tags, length(pick), replace = TRUE)
  data.table(Probe = pick, Gene = genes_rep, Region = regions)
}
