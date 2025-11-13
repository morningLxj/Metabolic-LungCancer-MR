# GitHub 推送指南

## 仓库准备状态

您的代谢肺癌MR分析仓库已经完全准备好推送到GitHub！

### 当前状态
- ✅ Git仓库已初始化
- ✅ 所有文件已添加并提交
- ✅ 主分支已设置为 'main'
- ✅ 所有配置文件已创建

### 推送到GitHub的步骤

#### 方法1：使用GitHub CLI（推荐）

```bash
# 1. 创建GitHub仓库
gh repo create metabolic-lung-cancer-mr-analysis --public --description "Comprehensive Mendelian Randomization analysis pipeline for metabolic traits and lung cancer risk"

# 2. 推送代码
git push -u origin main
```

#### 方法2：使用GitHub网页界面

1. 访问 [GitHub](https://github.com) 并登录
2. 点击 "New repository"
3. 仓库名称：`metabolic-lung-cancer-mr-analysis`
4. 描述：`Comprehensive Mendelian Randomization analysis pipeline for metabolic traits and lung cancer risk`
5. 设置为 Public（公开）
6. **不要**勾选 "Add a README file"（我们已经有了）
7. 点击 "Create repository"
8. 复制仓库URL
9. 在终端中运行：

```bash
git remote add origin https://github.com/YOUR_USERNAME/metabolic-lung-cancer-mr-analysis.git
git push -u origin main
```

#### 方法3：使用HTTPS

如果需要输入凭据：

```bash
git remote add origin https://github.com/YOUR_USERNAME/metabolic-lung-cancer-mr-analysis.git
git push -u origin main
```

### 仓库特性

此仓库包含：

1. **完整的MR分析流水线**
   - 数据预处理和质量控制
   - 工具变量提取
   - 多种MR方法（IVW、MR-Egger、加权中位数等）
   - 敏感性分析

2. **工具函数库**
   - `scripts/utils/data_processing.R` - 数据处理
   - `scripts/utils/mr_analysis.R` - MR分析
   - `scripts/utils/data_qc.R` - 质量控制
   - `scripts/utils/visualization.R` - 可视化

3. **测试套件**
   - 单元测试：`tests/test_main_functions.R`
   - 集成测试：`tests/integration_tests.R`

4. **CI/CD管道**
   - GitHub Actions工作流
   - 自动化测试和代码检查

5. **完整文档**
   - 详细的README.md
   - 配置说明
   - 使用示例

### 后续步骤

推送成功后，您可以：

1. 设置GitHub Pages（如果需要文档网站）
2. 配置分支保护规则
3. 设置Issue模板
4. 创建Release版本
5. 邀请协作者

### 仓库结构

```
metabolic-lung-cancer-mr-analysis/
├── README.md                 # 项目文档
├── LICENSE                   # 开源协议
├── requirements.txt          # Python依赖
├── R_packages.R             # R包安装脚本
├── .gitignore               # Git忽略规则
├── .github/workflows/ci.yml # CI/CD配置
├── docs/                    # 文档目录
├── scripts/                 # 分析脚本
│   ├── R/                  # R脚本
│   ├── python/             # Python脚本
│   ├── utils/              # 工具函数
│   ├── main_analysis.R     # 主执行脚本
│   ├── config.yaml         # 配置文件
│   └── run_example.R       # 示例脚本
├── tests/                   # 测试文件
├── results/                 # 结果输出目录
└── data/                    # 数据目录（可选）
```

### 支持

如果在推送过程中遇到问题，请检查：

1. GitHub认证是否正确设置
2. 网络连接是否正常
3. 仓库名称是否已存在

祝您推送成功！