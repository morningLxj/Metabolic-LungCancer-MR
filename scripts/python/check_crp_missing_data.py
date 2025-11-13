#!/usr/bin/env python3
"""
检查CRP数据中的缺失值
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def check_missing_data(df, file_name):
    """检查数据框中的缺失值"""
    print(f"\n=== 检查文件: {file_name} ===")
    print(f"数据形状: {df.shape}")
    
    # 检查各列的缺失值
    missing_counts = df.isnull().sum()
    missing_percentages = (df.isnull().sum() / len(df)) * 100
    
    missing_summary = pd.DataFrame({
        '缺失数量': missing_counts,
        '缺失百分比': missing_percentages
    })
    
    # 只显示有缺失值的列
    missing_cols = missing_summary[missing_summary['缺失数量'] > 0]
    
    if len(missing_cols) > 0:
        print("\n存在缺失值的列:")
        print(missing_cols[missing_cols['缺失数量'] > 0].sort_values('缺失数量', ascending=False))
    else:
        print("\n✓ 该文件没有缺失值")
    
    # 显示数据类型信息
    print(f"\n数据类型信息:")
    print(df.dtypes.value_counts())
    
    return missing_summary

def main():
    """主函数：检查所有CRP相关数据文件"""
    
    # 数据文件路径
    data_dir = Path("data")
    processed_dir = data_dir / "processed"
    raw_dir = data_dir / "raw"
    
    print("开始检查CRP数据的缺失值...")
    
    # 1. 检查原始CRP基因文件
    crp_genes_file = raw_dir / "CRP-genes.txt"
    if crp_genes_file.exists():
        try:
            crp_genes_df = pd.read_csv(crp_genes_file, sep='\t')
            check_missing_data(crp_genes_df, "CRP-genes.txt")
        except Exception as e:
            print(f"读取 CRP-genes.txt 时出错: {e}")
    else:
        print("未找到 CRP-genes.txt 文件")
    
    # 2. 检查处理后的GWAS数据
    gwas_crp_file = processed_dir / "gwas" / "crp_gwas.parquet"
    if gwas_crp_file.exists():
        try:
            crp_gwas_df = pd.read_parquet(gwas_crp_file)
            check_missing_data(crp_gwas_df, "crp_gwas.parquet")
        except Exception as e:
            print(f"读取 crp_gwas.parquet 时出错: {e}")
    else:
        print("未找到 crp_gwas.parquet 文件")
    
    # 3. 检查候选基因汇总文件
    candidate_genes_file = processed_dir / "all_candidate_genes.parquet"
    if candidate_genes_file.exists():
        try:
            candidate_df = pd.read_parquet(candidate_genes_file)
            # 如果有CRP相关数据，单独检查
            if 'trait' in candidate_df.columns:
                crp_candidate = candidate_df[candidate_df['trait'].str.contains('CRP', case=False, na=False)]
                if len(crp_candidate) > 0:
                    check_missing_data(crp_candidate, "all_candidate_genes.parquet (CRP相关数据)")
                else:
                    print("\nall_candidate_genes.parquet 中没有找到CRP相关数据")
            else:
                check_missing_data(candidate_df, "all_candidate_genes.parquet")
        except Exception as e:
            print(f"读取 all_candidate_genes.parquet 时出错: {e}")
    else:
        print("未找到 all_candidate_genes.parquet 文件")
    
    # 3.5 检查候选基因汇总文件 (使用不同的编码)
    candidate_genes_file2 = processed_dir / "all_candidate_genes.parquet"
    if candidate_genes_file2.exists():
        try:
            # 尝试重新读取，看看是否是临时问题
            candidate_df2 = pd.read_parquet(candidate_genes_file2)
            check_missing_data(candidate_df2, "all_candidate_genes.parquet (重试)")
        except Exception as e:
            print(f"重试读取 all_candidate_genes.parquet 时出错: {e}")
    
    # 4. 检查分析任务文件
    tasks_file = processed_dir / "analysis_tasks.txt"
    if tasks_file.exists():
        try:
            tasks_df = pd.read_csv(tasks_file, sep='\t')
            check_missing_data(tasks_df, "analysis_tasks.txt")
        except Exception as e:
            print(f"读取 analysis_tasks.txt 时出错: {e}")
    else:
        print("未找到 analysis_tasks.txt 文件")
    
    # 5. 检查其他可能的CRP数据文件
    print("\n=== 检查其他可能的数据文件 ===")
    potential_files = [
        "data/processed/gwas/crp_gwas.parquet",
        "data/raw/gwas/CRP-genes.txt"
    ]
    
    for file_path in potential_files:
        if os.path.exists(file_path):
            try:
                if file_path.endswith('.parquet'):
                    df = pd.read_parquet(file_path)
                elif file_path.endswith('.txt'):
                    df = pd.read_csv(file_path, sep='\t')
                else:
                    continue
                check_missing_data(df, file_path)
            except Exception as e:
                print(f"读取 {file_path} 时出错: {e}")
        else:
            print(f"未找到: {file_path}")

if __name__ == "__main__":
    main()