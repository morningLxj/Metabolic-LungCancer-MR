suppressMessages(library(data.table))
tcga_pid <- function(x) sub("^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*$","\\1",gsub("\\n","",x))
tcga_type <- function(barcode) { m <- sub("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2}).*$","\\1",gsub("\\s","",barcode)); ifelse(m %in% c("01","02","03","05","06","07","08","09","12","13"),"Tumor",ifelse(m %in% c("10","11"),"Normal",NA_character_)) }
load_methylation <- function(path) { if (!file.exists(path)) stop("missing methylation file: ", path); dt <- fread(path); if (!"Composite Element REF" %in% names(dt)) stop("invalid methylation schema"); setnames(dt,"Composite Element REF","Probe"); dt }
load_rna <- function(path) { if (!file.exists(path)) stop("missing rna file: ", path); dt <- fread(path); if (!"Ensembl_ID" %in% names(dt)) stop("invalid rna schema"); setnames(dt,"Ensembl_ID","ENSEMBL"); dt[, ENSEMBL := sub("\\\\.[0-9]+$","", ENSEMBL)]; dt }
