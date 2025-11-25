#!/usr/bin/env python3
"""
全面检查第3步共定位分析的数据依赖
"""

import pandas as pd
import os
from pathlib import Path
import numpy as np

def check_coloc_data_requirements():
    """检查共定位分析所需的所有数据文件"""
    print("=== 共定位分析数据需求检查 ===")
    print("基于Script 03_run_coloc.R的依赖分析")
    print()
    
    # 必需的文件列表
    required_files = {
        # 1. GWAS数据文件 (必须存在)
        "GWAS - BMI": "data/processed/gwas/bmi_gwas.parquet",
        "GWAS - CRP": "data/processed/gwas/crp_gwas.parquet", 
        "GWAS - WBC": "data/processed/gwas/wbc_gwas.parquet",
        "GWAS - Alcohol": "data/processed/gwas/alcohol_gwas.parquet",
        "GWAS - LUAD": "data/processed/gwas/luad_gwas.parquet",
        "GWAS - LUSC": "data/processed/gwas/lusc_gwas.parquet",
        
        # 2. eQTL数据文件 (关键依赖)
        "eQTL - Lung": "data/processed/eqtl/lung_eqtl_annotated.parquet",
        
        # 3. 分析任务文件
        "Analysis Tasks": "data/processed/analysis_tasks.txt"
    }
    
    missing_files = []
    existing_files = []
    
    print("1. 检查必需文件:")
    for name, path in required_files.items():
        if os.path.exists(path):
            size_mb = os.path.getsize(path) / (1024*1024)
            existing_files.append((name, path, size_mb))
            print(f"   [OK] {name}: {path} ({size_mb:.1f} MB)")
        else:
            missing_files.append((name, path))
            print(f"   [MISSING] {name}: {path} - 文件不存在!")
    
    print(f"\n文件状态: {len(existing_files)}/{len(required_files)} 个文件存在")
    
    if missing_files:
        print(f"\n缺失的关键文件 ({len(missing_files)}个):")
        for name, path in missing_files:
            print(f"   - {name}: {path}")
        return False, missing_files, existing_files
    
    return True, missing_files, existing_files

def check_gwas_data_quality(file_path, trait_name):
    """检查GWAS数据质量"""
    print(f"\n2. 检查GWAS数据质量 - {trait_name}:")
    
    try:
        df = pd.read_parquet(file_path)
        print(f"   数据形状: {df.shape}")
        
        # 必需列检查
        required_cols = ['chr', 'pos', 'ref', 'alt', 'rsid', 'eaf', 'beta', 'se', 'pval', 'n']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"   [ERROR] 缺少必需列: {missing_cols}")
            return False
        
        # 缺失值检查
        missing_summary = df.isnull().sum()
        total_missing = missing_summary.sum()
        
        if total_missing > 0:
            print(f"   [ERROR] 发现缺失值:")
            for col, count in missing_summary.items():
                if count > 0:
                    pct = (count / len(df)) * 100
                    print(f"      {col}: {count:,} ({pct:.1f}%)")
            return False
        else:
            print(f"   [OK] 无缺失值")
        
        # 数据范围检查
        print(f"   [OK] 数据范围检查:")
        print(f"      染色体: {df['chr'].min()}-{df['chr'].max()}")
        print(f"      位置: {df['pos'].min():,}-{df['pos'].max():,}")
        print(f"      beta范围: {df['beta'].min():.4f}到{df['beta'].max():.4f}")
        print(f"      pval范围: {df['pval'].min():.2e}到{df['pval'].max():.2e}")
        print(f"      eaf范围: {df['eaf'].min():.3f}到{df['eaf'].max():.3f}")
        
        return True
        
    except Exception as e:
        print(f"   [ERROR] 读取文件失败: {e}")
        return False

def check_eqtl_data_quality():
    """检查eQTL数据质量"""
    print(f"\n3. 检查eQTL数据质量:")
    
    file_path = "data/processed/eqtl/lung_eqtl_annotated.parquet"
    
    try:
        df = pd.read_parquet(file_path)
        print(f"   数据形状: {df.shape}")
        print(f"   列名: {list(df.columns)}")
        
        # 必需列检查
        # 根据R脚本，eQTL数据需要包含基因信息用于提取区域
        required_base_cols = ['chr', 'pos', 'gene_symbol']
        missing_base_cols = [col for col in required_base_cols if col not in df.columns]
        
        if missing_base_cols:
            print(f"   [WARNING] 缺少基础列: {missing_base_cols}")
            
            # 尝试找到相似的列名
            possible_chr_cols = [col for col in df.columns if 'chr' in col.lower() or 'chrom' in col.lower()]
            possible_gene_cols = [col for col in df.columns if 'gene' in col.lower() or 'symbol' in col.lower()]
            possible_pos_cols = [col for col in df.columns if 'pos' in col.lower() or 'position' in col.lower()]
            
            print(f"   [SEARCH] 可能的替代列:")
            print(f"      染色体: {possible_chr_cols}")
            print(f"      基因: {possible_gene_cols}")
            print(f"      位置: {possible_pos_cols}")
        
        # 检查基因数量
        if 'gene_symbol' in df.columns:
            unique_genes = df['gene_symbol'].nunique()
            print(f"   [OK] 包含 {unique_genes:,} 个独特基因")
        else:
            print(f"   [ERROR] 没有找到gene_symbol列")
            return False
        
        # 缺失值检查
        missing_summary = df.isnull().sum()
        total_missing = missing_summary.sum()
        
        if total_missing > 0:
            print(f"   [WARNING] 发现缺失值:")
            for col, count in missing_summary.items():
                if count > 0 and col in required_base_cols:
                    pct = (count / len(df)) * 100
                    print(f"      {col}: {count:,} ({pct:.1f}%)")
        
        return True
        
    except Exception as e:
        print(f"   [ERROR] 读取文件失败: {e}")
        return False

def check_analysis_tasks():
    """检查分析任务文件"""
    print(f"\n4. 检查分析任务文件:")
    
    file_path = "data/processed/analysis_tasks.txt"
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"   数据形状: {df.shape}")
        print(f"   列名: {list(df.columns)}")
        
        # 必需列检查
        required_cols = ['exposure', 'outcome', 'gene_symbol', 'gene_id', 'chr', 'start', 'end']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"   [ERROR] 缺少必需列: {missing_cols}")
            return False
        
        # 检查exposure和outcome的组合
        exposures = df['exposure'].unique()
        outcomes = df['outcome'].unique()
        
        print(f"   [OK] Exposure traits: {list(exposures)}")
        print(f"   [OK] Outcome traits: {list(outcomes)}")
        print(f"   [OK] 总分析任务: {len(df):,}")
        
        # 检查是否覆盖了所有GWAS数据
        expected_exposures = ['BMI', 'CRP', 'WBC', 'Alcohol']
        expected_outcomes = ['LUAD', 'LUSC']
        
        missing_exposures = set(expected_exposures) - set(exposures)
        missing_outcomes = set(expected_outcomes) - set(outcomes)
        
        if missing_exposures:
            print(f"   [WARNING] 缺少exposure trait: {list(missing_exposures)}")
        if missing_outcomes:
            print(f"   [WARNING] 缺少outcome trait: {list(missing_outcomes)}")
        
        return True
        
    except Exception as e:
        print(f"   [ERROR] 读取文件失败: {e}")
        return False

def main():
    """主检查函数"""
    print("共定位分析数据完整性检查")
    print("="*50)
    
    # 1. 检查必需文件
    files_ok, missing_files, existing_files = check_coloc_data_requirements()
    
    if not files_ok:
        print(f"\n[ERROR] 缺少关键文件，无法进行共定位分析")
        print("请先解决缺失文件问题再继续")
        return False
    
    # 2. 检查GWAS数据质量
    print(f"\n" + "="*50)
    print("数据质量检查")
    print("="*50)
    
    gwas_traits = ['BMI', 'CRP', 'WBC', 'Alcohol', 'LUAD', 'LUSC']
    gwas_ok_count = 0
    
    for trait in gwas_traits:
        file_path = f"data/processed/gwas/{trait.lower()}_gwas.parquet"
        if check_gwas_data_quality(file_path, trait):
            gwas_ok_count += 1
    
    print(f"\nGWAS数据状态: {gwas_ok_count}/{len(gwas_traits)} 个文件通过质量检查")
    
    # 3. 检查eQTL数据
    eqtl_ok = check_eqtl_data_quality()
    
    # 4. 检查分析任务
    tasks_ok = check_analysis_tasks()
    
    # 5. 总结
    print(f"\n" + "="*50)
    print("检查总结")
    print("="*50)
    
    all_ok = files_ok and (gwas_ok_count == len(gwas_traits)) and eqtl_ok and tasks_ok
    
    if all_ok:
        print("[SUCCESS] 所有数据完整且质量良好，共定位分析就绪！")
        return True
    else:
        print("[WARNING] 发现数据问题，需要修复后进行共定位分析")
        
        if gwas_ok_count < len(gwas_traits):
            print(f"   - {len(gwas_traits) - gwas_ok_count}个GWAS文件有问题")
        if not eqtl_ok:
            print("   - eQTL数据有问题")
        if not tasks_ok:
            print("   - 分析任务文件有问题")
        
        return False

if __name__ == "__main__":
    main()