#!/usr/bin/env python3
"""
简化版CRP数据缺失值检查
"""

import pandas as pd
import sys

def check_file_safely(file_path, file_desc):
    """安全地检查文件缺失值"""
    try:
        print(f"\n=== {file_desc} ===")
        print(f"文件路径: {file_path}")
        
        # 尝试读取文件
        if file_path.endswith('.parquet'):
            df = pd.read_parquet(file_path)
        elif file_path.endswith('.txt'):
            df = pd.read_csv(file_path, sep='\t')
        else:
            print("不支持的文件格式")
            return
        
        print(f"数据形状: {df.shape}")
        print(f"列名: {list(df.columns)}")
        
        # 检查缺失值
        missing_counts = df.isnull().sum()
        total_missing = missing_counts.sum()
        
        print(f"总缺失值数量: {total_missing}")
        
        if total_missing > 0:
            print("\n有缺失值的列:")
            for col, count in missing_counts.items():
                if count > 0:
                    pct = (count / len(df)) * 100
                    print(f"  {col}: {count} 个缺失值 ({pct:.2f}%)")
        else:
            print("✓ 没有缺失值")
        
        # 数据类型
        print(f"\n数据类型分布:")
        print(df.dtypes.value_counts())
        
        return df
        
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return None

def main():
    print("开始简化版CRP数据缺失值检查...")
    
    # 检查列表
    files_to_check = [
        ("data/processed/all_candidate_genes.parquet", "候选基因文件"),
        ("data/processed/analysis_tasks.txt", "分析任务文件"),
        ("data/processed/gwas/crp_gwas.parquet", "CRP GWAS数据"),
    ]
    
    for file_path, desc in files_to_check:
        check_file_safely(file_path, desc)
    
    print("\n=== 汇总 ===")
    print("检查完成")

if __name__ == "__main__":
    main()