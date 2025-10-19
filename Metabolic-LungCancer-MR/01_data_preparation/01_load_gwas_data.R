# 加载所有GWAS数据
cat("开始加载GWAS数据...\n")

# 初始化存储列表
exposure_data_list <- list()
outcome_data_list <- list()

# 加载暴露数据
cat("加载暴露数据...\n")
for (trait_name in names(STUDY_VARIABLES$exposures)) {
  gwas_id <- STUDY_VARIABLES$exposures[[trait_name]]
  exposure_data <- load_gwas_data(gwas_id, trait_name, "exposure")
  
  if (!is.null(exposure_data)) {
    exposure_data_list[[trait_name]] <- exposure_data
    cat("✓ 成功加载:", trait_name, "-", nrow(exposure_data), "个SNP\n")
  } else {
    cat("✗ 加载失败:", trait_name, "\n")
  }
}

# 加载结局数据  
cat("加载结局数据...\n")
for (trait_name in names(STUDY_VARIABLES$outcomes)) {
  gwas_id <- STUDY_VARIABLES$outcomes[[trait_name]]
  outcome_data <- load_gwas_data(gwas_id, trait_name, "outcome")
  
  if (!is.null(outcome_data)) {
    outcome_data_list[[trait_name]] <- outcome_data
    cat("✓ 成功加载:", trait_name, "\n")
  } else {
    cat("✗ 加载失败:", trait_name, "\n")
  }
}

# 保存加载的数据
saveRDS(exposure_data_list, "01_data_preparation/exposure_data_list.rds")
saveRDS(outcome_data_list, "01_data_preparation/outcome_data_list.rds")

cat("GWAS数据加载完成!\n")
cat("成功加载暴露:", length(exposure_data_list), "个\n")
cat("成功加载结局:", length(outcome_data_list), "个\n")