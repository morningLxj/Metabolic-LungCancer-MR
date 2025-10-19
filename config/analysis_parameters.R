# 孟德尔随机化分析参数配置
# 基于项目数据集表格

# 分析参数
ANALYSIS_PARAMS <- list(
  
  # 工具变量筛选标准 (对应论文方法部分)
  iv_selection = list(
    pval_threshold = 5e-8,
    clump_r2 = 0.001,
    clump_kb = 10000,
    remove_palindromic = TRUE
  ),
  
  # MR分析方法
  mr_methods = c("mr_ivw", "mr_egger", "mr_weighted_median"),
  
  # 多重检验校正
  multiple_testing = list(
    bonferroni_threshold = 0.05/51,  # 17个暴露 × 3个结局
    fdr_threshold = 0.05
  ),
  
  # F统计量阈值
  f_statistic_threshold = 10
)

# 研究变量定义 (直接从您的表格中提取)
STUDY_VARIABLES <- list(
  exposures = list(
    `body mass index` = "ieu-b-40",
    `LDL cholesterol` = "ieu-b-110", 
    `HDL cholesterol` = "ieu-b-109",
    `Apolipoprotein A1 levels` = "ebi-a-GCST90025955",
    `Apolipoprotein B levels` = "ebi-a-GCST90025952",
    `Ratio of apolipoprotein B to apolipoprotein A1 levels` = "ebi-a-GCST90092810",
    `C-reactive protein levels` = "ebi-a-GCST90029070",
    `circulating leptin levels` = "ebi-a-GCST90007316",
    `Serum 25-Hydroxyvitamin D levels` = "ebi-a-GCST90000618",
    `Glycated haemoglobin HbA1c levels` = "ebi-a-GCST90014006",
    `Insulin-like growth factor 1 levels` = "ebi-a-GCST90025989",
    `Concentration of small LDL particles` = "ebi-a-GCST90092963",
    `Concentration of large HDL particles` = "ebi-a-GCST90092851",
    `Concentration of very large HDL particles` = "ebi-a-GCST90093011",
    `Average diameter for HDL particles` = "ebi-a-GCST90092828",
    `Remnant cholesterol` = "ebi-a-GCST90092943",
    `Total concentration of branched-chain amino acids` = "ebi-a-GCST90092984"
  ),
  
  outcomes = list(
    `Overall Lung Cancer` = "ebi-a-GCST90018875",
    `Lung adenocarcinoma` = "ieu-a-984", 
    `Squamous cell lung cancer` = "ieu-a-989"
  )
)

# 设置随机种子保证可重复性
set.seed(2024)