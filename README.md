# 代谢肺癌MR分析仓库

**Metabolic Lung Cancer Mendelian Randomization Analysis Repository**

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 项目概述

本仓库提供了一个完整的Mendelian Randomization (MR) 分析流程，专门用于研究代谢特征与肺癌之间的因果关系。该分析流程包含数据质量控制、工具变量提取、双样本MR分析、敏感性分析和结果可视化等功能。

### 主要特征

- 🔬 **双样本MR分析**: 支持多种MR方法（IVW、MR-Egger、加权中位数等）
- 📊 **数据质量控制**: 全面的GWAS数据QC流程
- 🔍 **敏感性分析**: 包含留一法、MR-PRESSO、水平多效性检验
- 📈 **结果可视化**: 森林图、漏斗图、方法比较图
- 📝 **自动化报告**: 生成详细的分析报告和QC报告
- ⚡ **高性能**: 支持并行处理，适合大规模GWAS数据分析
- 🧪 **单元测试**: 完整的测试套件确保代码质量

## 快速开始

### 环境要求

- R >= 4.0
- Python >= 3.7 (可选，用于某些工具)
- 至少8GB内存（推荐16GB+）

### 安装依赖

首先安装所有必需的R包：

```r
# 安装依赖包
source("R_packages.R")

# 或者手动安装核心包
install.packages(c(
  "MendelianRandomization", "TwoSampleMR", "ieugwasr", 
  "dplyr", "ggplot2", "forestplot", "optparse"
))
```

### 运行示例分析

使用提供的示例脚本快速开始：

```bash
# 基本用法（使用默认参数和演示数据）
Rscript scripts/run_example.R

# 指定自定义参数
Rscript scripts/run_example.R \
  --exposure BMI \
  --outcome Lung_Cancer \
  --threads 4 \
  --method "all"
```

### 参数说明

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `--exposure` | 暴露变量名称 | BMI |
| `--outcome` | 结果变量名称 | Lung_Cancer |
| `--exposure_file` | 暴露数据文件路径 | data/raw/exposure_gwas.txt |
| `--outcome_file` | 结果数据文件路径 | data/raw/outcome_gwas.txt |
| `--output` | 输出目录 | results/ |
| `--threads` | 并行线程数 | 1 |
| `--method` | MR方法选择 | all |
| `--sensitivity` | 是否进行敏感性分析 | TRUE |

## 仓库结构

```
├── README.md                     # 项目文档
├── R_packages.R                  # R包依赖安装脚本
├── scripts/                      # 分析脚本目录
│   ├── main_analysis.R          # 主分析执行脚本
│   ├── run_example.R            # 示例运行脚本
│   ├── config.yaml              # 配置文件
│   └── utils/                   # 工具函数
│       ├── data_processing.R    # 数据处理工具
│       ├── mr_analysis.R        # MR分析工具
│       ├── data_qc.R            # 数据质量控制
│       └── visualization.R      # 可视化工具
├── data/                        # 数据目录
│   ├── raw/                     # 原始数据
│   ├── processed/               # 处理后数据
│   └── reference/               # 参考数据
├── results/                     # 分析结果
│   ├── tables/                  # 结果表格
│   ├── plots/                   # 图表文件
│   └── reports/                 # 分析报告
├── docs/                        # 文档目录
│   ├── methods/                 # 方法学文档
│   └── tutorials/               # 教程文档
├── tests/                       # 测试文件
│   ├── test_main_functions.R    # 单元测试
│   └── integration_tests.R      # 集成测试
├── .gitignore                   # Git忽略文件
└── LICENSE                      # 许可证
```

## 快速开始

### 环境要求

- R >= 4.0
- Python >= 3.8
- R包: TwoSampleMR, MendelianRandomization, ggplot2
- Python包: pandas, numpy, matplotlib, scipy

### 安装依赖

```bash
# 安装Python依赖
pip install -r requirements.txt

# 安装R包 (在R中运行)
install.packages(c("TwoSampleMR", "MendelianRandomization", "ggplot2", "dplyr"))
```

### 运行分析

1. **主MR分析**
   ```bash
   cd scripts/R
   # 运行步骤1-12的分析脚本
   Rscript step01_data_preparation.R
   Rscript step02_exposure_gwas.R
   # ... 继续后续步骤
   ```

2. **富集分析**
   ```bash
   cd scripts/python
   python enrichment_analysis.py
   ```

3. **共定位分析**
   ```bash
   cd scripts/python
   python coloc_analysis.py
   ```

## 主要结果

### MR分析结果
- 主要分析结果表格
- 森林图和漏斗图
- 敏感性分析结果

### 富集分析结果
- FUMA功能注释结果
- 通路富集分析
- 基因集富集分析

### 共定位分析结果
- 共定位概率分析
- 基因位点可视化
- 区域关联图

## 论文和发表

本项目支持完整的学术发表流程：

- 📄 完整的论文草稿
- 📊 60+个专业图表
- 📋 详细的方法学说明
- 📈 补充材料和分析结果

## 贡献

欢迎贡献！请阅读 [CONTRIBUTING.md](CONTRIBUTING.md) 了解详细信息。

## 许可证

本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 引用

如果您在研究中使用了本代码，请引用：

```
# TODO: 添加适当的引用信息
```

## 联系方式

- **项目维护者**: morningLxj
- **GitHub**: https://github.com/morningLxj/Metabolic-LungCancer-MR

---

**最后更新**: 2024年12月1日