#!/usr/bin/env python3
"""
检查原始eQTL数据源，寻找MAF信息
"""

import pandas as pd
import gzip
import os

def check_raw_eqtl_data():
    """检查原始eQTL数据的MAF信息"""
    print("=== 原始eQTL数据MAF检查 ===")
    
    # 检查主要原始文件
    raw_files = [
        "data/raw/eqtl/Lung.v8.egenes.txt.gz",
        "data/raw/eqtl/Lung.signifpairs.txt",
        "data/raw/eqtl/Lung.v8.signif_variant_gene_pairs.txt.gz"
    ]
    
    for file_path in raw_files:
        if not os.path.exists(file_path):
            print(f"[MISSING] {file_path}")
            continue
            
        print(f"\n[CHECK] {file_path}")
        
        try:
            if file_path.endswith('.gz'):
                # 读取压缩文件的头部
                with gzip.open(file_path, 'rt') as f:
                    header = f.readline().strip().split('\t')
                    sample_line = f.readline().strip().split('\t')
                    
                print(f"  列数: {len(header)}")
                print(f"  前10列: {header[:10]}")
                
                # 检查MAF相关列
                maf_cols = []
                for i, col in enumerate(header):
                    if any(keyword in col.lower() for keyword in ['maf', 'eaf', 'freq', 'allele', 'af']):
                        maf_cols.append((i, col))
                
                if maf_cols:
                    print(f"  [FOUND] MAF相关列: {maf_cols}")
                    
                    # 读取一些样本数据
                    with gzip.open(file_path, 'rt') as f:
                        f.readline()  # 跳过header
                        for _ in range(3):  # 读取3行样本
                            line = f.readline().strip().split('\t')
                            for col_idx, col_name in maf_cols:
                                if col_idx < len(line):
                                    print(f"    {col_name}: {line[col_idx]}")
                else:
                    print(f"  [NO] 未找到MAF相关列")
                    
            else:
                # 读取普通文本文件
                df = pd.read_csv(file_path, sep='\t', nrows=10)
                print(f"  形状: {df.shape}")
                print(f"  列名: {list(df.columns)}")
                
                # 检查MAF相关列
                maf_cols = [col for col in df.columns if any(keyword in col.lower() for keyword in ['maf', 'eaf', 'freq', 'allele', 'af'])]
                if maf_cols:
                    print(f"  [FOUND] MAF相关列: {maf_cols}")
                    print(f"  样本值: {df[maf_cols].head(3).to_dict()}")
                else:
                    print(f"  [NO] 未找到MAF相关列")
                    
        except Exception as e:
            print(f"  [ERROR] {e}")
    
    # 检查注释信息
    print(f"\n[NOTE] 如果原始数据没有MAF信息，可以考虑:")
    print(f"1. 使用1000 Genomes项目的MAF数据")
    print(f"2. 使用gnomAD数据库的MAF信息")
    print(f"3. 基于变异ID查询Ensembl数据库")
    print(f"4. 使用默认的合理MAF值 (0.1-0.9范围内)")

if __name__ == "__main__":
    check_raw_eqtl_data()