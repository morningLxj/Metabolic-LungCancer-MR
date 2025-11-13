# ==============================================================================
# Script: 01_install_packages_Robust.R
# Purpose: 安装所有需要的R包 (优化版，解决SSL连接问题)
# ==============================================================================

cat("
================================================================================
Installing Required R Packages for Colocalization Analysis (Robust Version)
================================================================================
\n")

# --- 核心修改点 1: 使用更稳定的镜像和下载方法 ---
# 设置CRAN镜像为RStudio官方镜像，通常更稳定且兼容性更好
options(repos = c(CRAN = "https://cran.rstudio.com/"))
# 在Windows上，设置下载方法为 'wininet' 可以有效解决很多SSL连接问题
if (.Platform$OS.type == "windows") {
  options(download.file.method = "wininet")
}

# 需要的CRAN包
cran_packages <- c(
  # 数据处理
  "tidyverse",      # 数据处理核心包
  "data.table",     # 快速数据读写
  "arrow",          # 大数据处理
  
  # 统计分析
  "coloc",          # 共定位分析
  "susieR",         # 精细定位
  
  # 遗传学
  "genetics",       # 遗传数据处理
  "snpStats",       # SNP数据处理
  
  # 可视化
  "ggplot2",        # 绘图
  "patchwork",      # 图形组合
  "pheatmap",       # 热图
  "ggrepel",        # 标签优化
  "ggsci",          # 科研配色
  "RColorBrewer",   # 颜色
  
  # 报告生成
  "knitr",          # 文档编译
  "rmarkdown",      # R Markdown
  "kableExtra",     # 表格美化
  "DT",             # 交互表格
  
  # 其他
  "parallel",       # 并行计算
  "foreach",        # 循环
  "doParallel",     # 并行后端
  "progress",       # 进度条
  "logger"          # 日志
)

# 需要的Bioconductor包
bioc_packages <- c(
  "VariantAnnotation",  # VCF处理
  "biomaRt",            # 基因注释
  "GenomicRanges",      # 基因组区间
  "rtracklayer",        # 基因组数据转换
  "org.Hs.eg.db"        # 人类基因注释
)

# ==============================================================================
# 安装CRAN包
# ==============================================================================

cat("\n1. Installing CRAN packages...
")

# --- 核心修改点 2: 使用 tryCatch 进行错误处理，避免因单个包失败而中断 ---
# 首先尝试一次性安装所有包，这通常能更好地处理依赖关系
cat("  Attempting to install all CRAN packages at once for better dependency resolution...
")
install_result <- try(
  install.packages(cran_packages, dependencies = TRUE),
  silent = TRUE
)

if (class(install_result) == "try-error") {
  cat("  Bulk installation failed or some packages are missing. Installing one by one...
")
  # 如果一次性安装失败，再逐个安装
  for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("  Installing %s...
", pkg))
      # 使用 tryCatch 包裹，即使失败也继续下一个
      try(install.packages(pkg, dependencies = TRUE), silent = TRUE)
    } else {
      cat(sprintf("  ✓ %s already installed\n", pkg))
    }
  }
} else {
  cat("  ✓ Bulk installation successful.
")
}


# ==============================================================================
# 安装Bioconductor包
# ==============================================================================

cat("
2. Installing Bioconductor packages...
")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 同样为Bioconductor包添加错误处理
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...
", pkg))
    try(BiocManager::install(pkg, update = FALSE, ask = FALSE), silent = TRUE)
  } else {
    cat(sprintf("  ✓ %s already installed
", pkg))
  }
}

# ==============================================================================
# 验证安装
# ==============================================================================

cat("
3. Verifying installation...
")

all_packages <- c(cran_packages, bioc_packages)
missing_packages <- c()
installed_packages <- c()

for (pkg in all_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    installed_packages <- c(installed_packages, pkg)
  } else {
    missing_packages <- c(missing_packages, pkg)
  }
}

cat(sprintf("  ✓ Successfully installed/verified %d packages.
", length(installed_packages)))

if (length(missing_packages) > 0) {
  cat("\n✗ Failed to install the following packages:
")
  cat(paste("  -", missing_packages, collapse = "\n"))
  cat("
\nThese packages may require manual installation or have system-level dependencies.
Please check the error messages above for more details.
")
} else {
  cat("\n✓ All packages successfully installed and verified!\n\n")
}

# ==============================================================================
# 显示R环境信息
# ==============================================================================

cat("
4. R Environment Information:
")
cat(sprintf("  R version: %s
", R.version.string))
cat(sprintf("  Platform: %s
", R.version$platform))
cat(sprintf("  Working directory: %s
", getwd()))

cat("
✓ Setup complete!\n\n")
