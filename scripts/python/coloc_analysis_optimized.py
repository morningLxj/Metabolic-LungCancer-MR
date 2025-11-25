import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.offline as pyo
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置样式
sns.set_style("whitegrid")
plt.style.use('seaborn-v0_8')

class ColocAnalysis:
    """
    共定位分析类，提供数据加载、处理和可视化功能
    """
    
    def __init__(self, data_path=None):
        """
        初始化共定位分析对象
        
        参数:
        data_path: 数据文件路径
        """
        self.data = None
        self.filtered_data = None
        self.stats = {}
        
        if data_path:
            self.load_data(data_path)
    
    def load_data(self, file_path):
        """
        加载共定位分析数据
        
        参数:
        file_path: 数据文件路径
        """
        try:
            self.data = pd.read_csv(file_path, sep='\t')
            print(f"成功加载数据，共 {len(self.data)} 条记录")
            print(f"数据列: {list(self.data.columns)}")
            self.calculate_basic_stats()
            return True
        except Exception as e:
            print(f"加载数据失败: {e}")
            return False
    
    def calculate_basic_stats(self):
        """
        计算基本统计信息
        """
        if self.data is None:
            return
            
        # 检查必要的列是否存在
        if 'gene' not in self.data.columns:
            print("警告: 数据中缺少'gene'列")
            return
            
        if 'PP4' not in self.data.columns:
            # 如果没有PP4列，尝试使用exp_out_PP.H4
            if 'exp_out_PP.H4' in self.data.columns:
                self.data['PP4'] = self.data['exp_out_PP.H4']
            else:
                print("警告: 数据中缺少'PP4'或'exp_out_PP.H4'列")
                return
        
        # 确保有trait列
        if 'trait' not in self.data.columns:
            if 'exposure' in self.data.columns:
                self.data['trait'] = self.data['exposure']
            else:
                print("警告: 数据中缺少'trait'或'exposure'列")
                return
        
        # 确保有pvalue列
        if 'pvalue' not in self.data.columns:
            if 'exp_out_PP.H0' in self.data.columns:
                self.data['pvalue'] = self.data['exp_out_PP.H0']
            else:
                print("警告: 数据中缺少'pvalue'或'exp_out_PP.H0'列")
                return
        
        self.stats = {
            'total_records': len(self.data),
            'unique_genes': self.data['gene'].nunique(),
            'unique_traits': self.data['trait'].nunique(),
            'mean_pp4': self.data['PP4'].mean(),
            'median_pp4': self.data['PP4'].median(),
            'high_pp4_count': (self.data['PP4'] > 0.8).sum(),
            'moderate_pp4_count': ((self.data['PP4'] > 0.5) & (self.data['PP4'] <= 0.8)).sum(),
            'low_pp4_count': (self.data['PP4'] <= 0.5).sum()
        }
    
    def filter_data(self, pp4_threshold=0.5, pvalue_threshold=None, genes=None, traits=None):
        """
        根据条件过滤数据
        
        参数:
        pp4_threshold: PP4阈值，默认0.5
        pvalue_threshold: P值阈值
        genes: 基因列表
        traits: 性状列表
        
        返回:
        过滤后的数据
        """
        if self.data is None:
            print("请先加载数据")
            return None
            
        filtered_data = self.data.copy()
        
        # 根据PP4阈值过滤
        filtered_data = filtered_data[filtered_data['PP4'] >= pp4_threshold]
        
        # 根据P值阈值过滤
        if pvalue_threshold is not None:
            filtered_data = filtered_data[filtered_data['pvalue'] <= pvalue_threshold]
        
        # 根据基因过滤
        if genes is not None:
            filtered_data = filtered_data[filtered_data['gene'].isin(genes)]
        
        # 根据性状过滤
        if traits is not None:
            filtered_data = filtered_data[filtered_data['trait'].isin(traits)]
        
        # 保存过滤后的数据
        self.filtered_data = filtered_data
        
        print(f"过滤后数据: {len(filtered_data)} 条记录")
        return filtered_data
    
    def create_interactive_scatter(self, x='PP4', y='pvalue', color='gene', size=None, 
                                  hover_data=None, title="共定位分析交互式散点图"):
        """
        创建交互式散点图
        
        参数:
        x: X轴变量
        y: Y轴变量
        color: 颜色分组变量
        size: 点大小变量
        hover_data: 悬停显示的数据
        title: 图表标题
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        if hover_data is None:
            hover_data = ['gene', 'trait', 'PP4', 'pvalue']
        
        fig = px.scatter(
            data_to_use, 
            x=x, 
            y=y, 
            color=color,
            size=size,
            hover_data=hover_data,
            title=title,
            color_continuous_scale=px.colors.sequential.Viridis,
            template='plotly_white'
        )
        
        # 更新布局
        fig.update_layout(
            xaxis_title=x,
            yaxis_title=y,
            font=dict(size=12),
            hovermode='closest',
            width=900,
            height=600
        )
        
        # 添加阈值线
        if x == 'PP4':
            fig.add_vline(x=0.5, line_dash="dash", line_color="orange", 
                         annotation_text="PP4=0.5阈值")
            fig.add_vline(x=0.8, line_dash="dash", line_color="red", 
                         annotation_text="PP4=0.8阈值")
        
        if y == 'pvalue':
            fig.add_hline(y=0.05, line_dash="dash", line_color="blue", 
                         annotation_text="P=0.05阈值")
        
        return fig
    
    def create_pp4_distribution(self, bins=30, title="PP4值分布"):
        """
        创建PP4值分布图
        
        参数:
        bins: 直方图箱数
        title: 图表标题
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=("PP4值直方图", "PP4值箱线图"),
            vertical_spacing=0.1
        )
        
        # 直方图
        fig.add_trace(
            go.Histogram(
                x=data_to_use['PP4'],
                nbinsx=bins,
                name="PP4分布",
                marker_color='royalblue',
                opacity=0.7
            ),
            row=1, col=1
        )
        
        # 箱线图
        fig.add_trace(
            go.Box(
                y=data_to_use['PP4'],
                name="PP4箱线图",
                marker_color='royalblue'
            ),
            row=2, col=1
        )
        
        # 添加阈值线
        fig.add_vline(x=0.5, line_dash="dash", line_color="orange", 
                     annotation_text="PP4=0.5阈值", row=1, col=1)
        fig.add_vline(x=0.8, line_dash="dash", line_color="red", 
                     annotation_text="PP4=0.8阈值", row=1, col=1)
        
        fig.update_layout(
            title=title,
            height=600,
            showlegend=False,
            template='plotly_white'
        )
        
        return fig
    
    def create_gene_trait_heatmap(self, title="基因-性状共定位热图"):
        """
        创建基因-性状共定位热图
        
        参数:
        title: 图表标题
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        # 创建透视表
        pivot_table = data_to_use.pivot_table(
            values='PP4', 
            index='gene', 
            columns='trait', 
            aggfunc='mean'
        )
        
        # 填充缺失值
        pivot_table = pivot_table.fillna(0)
        
        fig = px.imshow(
            pivot_table,
            labels=dict(x="性状", y="基因", color="PP4值"),
            x=pivot_table.columns,
            y=pivot_table.index,
            color_continuous_scale='Viridis',
            title=title
        )
        
        fig.update_layout(
            xaxis_title="性状",
            yaxis_title="基因",
            width=900,
            height=600,
            template='plotly_white'
        )
        
        return fig
    
    def create_top_genes_plot(self, top_n=10, metric='PP4', title="Top基因共定位分析"):
        """
        创建Top基因共定位分析图
        
        参数:
        top_n: 显示的基因数量
        metric: 排序指标
        title: 图表标题
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        # 按基因分组并计算平均PP4值
        gene_stats = data_to_use.groupby('gene').agg({
            metric: 'mean',
            'pvalue': 'mean',
            'trait': 'count'
        }).reset_index()
        
        # 按指标排序并取前N个
        top_genes = gene_stats.sort_values(by=metric, ascending=False).head(top_n)
        
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=(f"Top {top_n} 基因平均{metric}值", f"Top {top_n} 基因关联性状数量"),
            specs=[[{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # 平均PP4值条形图
        fig.add_trace(
            go.Bar(
                x=top_genes['gene'],
                y=top_genes[metric],
                name=f"平均{metric}值",
                marker_color='royalblue'
            ),
            row=1, col=1
        )
        
        # 关联性状数量条形图
        fig.add_trace(
            go.Bar(
                x=top_genes['gene'],
                y=top_genes['trait'],
                name="关联性状数量",
                marker_color='firebrick'
            ),
            row=1, col=2
        )
        
        fig.update_layout(
            title=title,
            height=500,
            showlegend=False,
            template='plotly_white'
        )
        
        fig.update_xaxes(tickangle=45)
        
        return fig
    
    def create_trait_network(self, min_pp4=0.5, title="性状关联网络图"):
        """
        创建性状关联网络图
        
        参数:
        min_pp4: 最小PP4阈值
        title: 图表标题
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        # 过滤数据
        filtered = data_to_use[data_to_use['PP4'] >= min_pp4]
        
        # 创建基因-性状关联表
        gene_trait_links = filtered[['gene', 'trait', 'PP4']]
        
        # 创建节点
        genes = gene_trait_links['gene'].unique()
        traits = gene_trait_links['trait'].unique()
        
        # 创建网络图
        fig = go.Figure()
        
        # 添加基因节点
        fig.add_trace(go.Scatter(
            x=[0] * len(genes),
            y=list(range(len(genes))),
            mode='markers+text',
            text=genes,
            textposition="middle center",
            marker=dict(size=20, color='skyblue'),
            name="基因"
        ))
        
        # 添加性状节点
        fig.add_trace(go.Scatter(
            x=[1] * len(traits),
            y=list(range(len(traits))),
            mode='markers+text',
            text=traits,
            textposition="middle center",
            marker=dict(size=20, color='lightcoral'),
            name="性状"
        ))
        
        # 添加连线
        for _, row in gene_trait_links.iterrows():
            gene_idx = list(genes).index(row['gene'])
            trait_idx = list(traits).index(row['trait'])
            
            fig.add_shape(
                type="line",
                x0=0, y0=gene_idx,
                x1=1, y1=trait_idx,
                line=dict(
                    width=row['PP4'] * 5,  # 线宽与PP4值成正比
                    color="rgba(0,0,0,0.5)"
                )
            )
        
        fig.update_layout(
            title=title,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            showlegend=True,
            height=600,
            template='plotly_white'
        )
        
        return fig
    
    def create_correlation_matrix(self, title="共定位分析相关性矩阵"):
        """
        创建相关性矩阵图
        
        参数:
        title: 图表标题
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        # 选择数值列
        numeric_cols = data_to_use.select_dtypes(include=[np.number]).columns
        
        # 计算相关性矩阵
        corr_matrix = data_to_use[numeric_cols].corr()
        
        fig = px.imshow(
            corr_matrix,
            labels=dict(x="变量", y="变量", color="相关系数"),
            x=corr_matrix.columns,
            y=corr_matrix.columns,
            color_continuous_scale='RdBu_r',
            title=title,
            zmin=-1,
            zmax=1
        )
        
        fig.update_layout(
            width=600,
            height=600,
            template='plotly_white'
        )
        
        return fig
    
    def generate_summary_report(self):
        """
        生成数据摘要报告
        """
        data_to_use = self.filtered_data if self.filtered_data is not None else self.data
        
        if data_to_use is None:
            print("请先加载数据")
            return None
        
        # 计算统计信息
        report = {
            '数据概览': {
                '总记录数': len(data_to_use),
                '唯一基因数': data_to_use['gene'].nunique(),
                '唯一性状数': data_to_use['trait'].nunique()
            },
            'PP4值统计': {
                '平均值': data_to_use['PP4'].mean(),
                '中位数': data_to_use['PP4'].median(),
                '最小值': data_to_use['PP4'].min(),
                '最大值': data_to_use['PP4'].max(),
                '标准差': data_to_use['PP4'].std()
            },
            'P值统计': {
                '平均值': data_to_use['pvalue'].mean(),
                '中位数': data_to_use['pvalue'].median(),
                '最小值': data_to_use['pvalue'].min(),
                '最大值': data_to_use['pvalue'].max(),
                '标准差': data_to_use['pvalue'].std()
            },
            'PP4分布': {
                '高共定位(>0.8)': (data_to_use['PP4'] > 0.8).sum(),
                '中等共定位(0.5-0.8)': ((data_to_use['PP4'] > 0.5) & (data_to_use['PP4'] <= 0.8)).sum(),
                '低共定位(≤0.5)': (data_to_use['PP4'] <= 0.5).sum()
            }
        }
        
        # 打印报告
        print("\n===== 共定位分析数据摘要报告 =====")
        for section, stats in report.items():
            print(f"\n{section}:")
            for key, value in stats.items():
                if isinstance(value, float):
                    print(f"  {key}: {value:.4f}")
                else:
                    print(f"  {key}: {value}")
        
        return report
    
    def save_figure(self, fig, filename, output_dir="results/figures", file_format="html"):
        """
        保存单个图表
        
        参数:
        fig: 图表对象
        filename: 文件名
        output_dir: 输出目录
        file_format: 文件格式 (html, png, pdf)
        """
        import os
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存图表
        output_path = os.path.join(output_dir, filename)
        if file_format == "html":
            fig.write_html(output_path)
        elif file_format == "png":
            fig.write_image(output_path)
        elif file_format == "pdf":
            fig.write_image(output_path)
        print(f"已保存图表: {output_path}")
    
    def save_all_plots(self, output_dir="./plots", file_format="html"):
        """
        保存所有图表
        
        参数:
        output_dir: 输出目录
        file_format: 文件格式 (html, png, pdf)
        """
        import os
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 创建并保存图表
        plots = {
            "interactive_scatter": self.create_interactive_scatter(),
            "pp4_distribution": self.create_pp4_distribution(),
            "gene_trait_heatmap": self.create_gene_trait_heatmap(),
            "top_genes": self.create_top_genes_plot(),
            "correlation_matrix": self.create_correlation_matrix()
        }
        
        for name, fig in plots.items():
            if fig is not None:
                output_path = os.path.join(output_dir, f"{name}.{file_format}")
                if file_format == "html":
                    fig.write_html(output_path)
                elif file_format == "png":
                    fig.write_image(output_path)
                elif file_format == "pdf":
                    fig.write_image(output_path)
                print(f"已保存图表: {output_path}")


def main():
    """
    主函数
    """
    # 创建分析实例
    analyzer = ColocAnalysis()
    
    # 加载数据
    if not analyzer.load_data("results/tables/coloc_results_full.txt"):
        return
    
    # 生成摘要报告
    analyzer.generate_summary_report()
    
    # 过滤数据
    filtered_data = analyzer.filter_data(pp4_threshold=0.5)
    print(f"过滤后数据: {len(filtered_data)} 条记录")
    
    # 创建可视化图表
    # 1. 交互式散点图
    scatter_fig = analyzer.create_interactive_scatter()
    analyzer.save_figure(scatter_fig, "interactive_scatter.html")
    
    # 2. PP4分布图
    pp4_hist_fig = analyzer.create_pp4_distribution()
    analyzer.save_figure(pp4_hist_fig, "pp4_distribution.html")
    
    # 3. 基因-性状热图
    heatmap_fig = analyzer.create_gene_trait_heatmap()
    analyzer.save_figure(heatmap_fig, "gene_trait_heatmap.html")
    
    print("所有图表已保存到 results/figures/ 目录")


if __name__ == "__main__":
    main()