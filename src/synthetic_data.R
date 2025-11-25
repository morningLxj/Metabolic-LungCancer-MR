suppressMessages(library(data.table))
gen_barcode <- function(n) sprintf("TCGA-%02s-%04d-01A-01D", sample(10:99,n,TRUE), sample(1000:9999,n,TRUE))
gen_probe <- function(n) paste0("cg", sprintf("%06d", sample(1e6, n)))
gen_beta <- function(n) pmin(pmax(rbeta(n, 2, 5), 0), 1)
gen_tpm <- function(n) pmax(rlnorm(n, 1, 0.8), 0)
write_methylation <- function(path, n_probe=500, n_sample=80) { probes <- gen_probe(n_probe); samples <- gen_barcode(n_sample); dt <- data.table(Probe = probes); for (s in samples) dt[[s]] <- gen_beta(n_probe); setnames(dt, "Probe", "Composite Element REF"); fwrite(dt, path) }
write_rna <- function(path, n_gene=500, n_sample=80) { ens <- paste0("ENSG", sprintf("%011d", sample(1e7, n_gene))); samples <- gen_barcode(n_sample); dt <- data.table(Ensembl_ID = ens); for (s in samples) dt[[s]] <- gen_tpm(n_gene); fwrite(dt, path) }
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
write_methylation("data/raw/TCGA-LUAD.methylation450.tsv")
write_methylation("data/raw/TCGA-LUSC.methylation450.tsv")
write_rna("data/raw/TCGA-LUAD.star_tpm.tsv")
write_rna("data/raw/TCGA-LUSC.star_tpm.tsv")
