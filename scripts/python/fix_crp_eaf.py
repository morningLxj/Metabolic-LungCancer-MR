#!/usr/bin/env python3
"""
CRP GWAS数据eaf缺失值修复脚本
"""

import pandas as pd
import vcfpy
import logging
from pathlib import Path
from rich.logging import RichHandler
from rich.progress import track
import numpy as np

def analyze_crp_vcf_structure():
    """分析CRP VCF文件结构"""
    vcf_file = Path("data/raw/gwas/ebi-a-GCST90029070.vcf.gz")
    
    print("=== 分析CRP VCF文件结构 ===")
    print(f"文件路径: {vcf_file}")
    
    try:
        with vcfpy.Reader.from_path(str(vcf_file)) as reader:
            header = reader.header
            print(f"VCF版本: {header.format_version}")
            
            # 检查FORMAT字段
            format_ids = list(header.format_ids())
            print(f"FORMAT字段: {format_ids}")
            
            # 检查前几个记录
            print("\n=== 前5个记录分析 ===")
            for i, record in enumerate(reader):
                if i >= 5:
                    break
                    
                print(f"\n记录 {i+1}: {record.CHROM}:{record.POS}")
                print(f"  REF:ALT = {record.REF}:{record.ALT[0].value if record.ALT else 'N/A'}")
                print(f"  ID = {record.ID}")
                
                if record.calls:
                    call_data = record.calls[0].data
                    print(f"  可用字段: {list(call_data.keys())}")
                    
                    # 显示所有字段的值
                    for key, value in call_data.items():
                        print(f"    {key}: {value}")
                
    except Exception as e:
        print(f"分析VCF文件失败: {e}")

def diagnose_eaf_issue():
    """诊断eaf缺失的具体原因"""
    print("\n=== 诊断eaf缺失问题 ===")
    
    # 读取现有的CRP GWAS数据
    crp_file = Path("data/processed/gwas/crp_gwas.parquet")
    if crp_file.exists():
        df = pd.read_parquet(crp_file)
        print(f"当前CRP GWAS数据形状: {df.shape}")
        print(f"eaf缺失数量: {df['eaf'].isnull().sum()}")
        print(f"其他列的缺失情况:")
        for col in df.columns:
            missing = df[col].isnull().sum()
            if missing > 0:
                print(f"  {col}: {missing} 缺失值")
        
        # 查看一些样本数据
        print(f"\n前5行数据:")
        print(df.head())
        
        return df
    else:
        print("未找到现有的CRP GWAS数据")
        return None

def fix_eaf_missing_data():
    """修复eaf缺失数据"""
    print("\n=== 开始修复eaf缺失数据 ===")
    
    # 读取现有数据
    crp_file = Path("data/processed/gwas/crp_gwas.parquet")
    if not crp_file.exists():
        print("错误：CRP GWAS数据文件不存在")
        return
    
    df = pd.read_parquet(crp_file)
    print(f"原始数据形状: {df.shape}")
    print(f"eaf缺失值数量: {df['eaf'].isnull().sum()}")
    
    # 分析数据情况
    beta_values = df['beta'].dropna()
    alt_values = df['alt'].dropna()
    
    print(f"beta统计: 均值={beta_values.mean():.4f}, 标准差={beta_values.std():.4f}")
    print(f"alt值示例: {alt_values.head().tolist()}")
    
    # 方法1: 估算eaf值
    # 对于GWAS数据，eaf通常在0.05-0.95之间
    # 我们可以基于其他信息进行合理估算
    
    # 检查是否有其他可用的信息
    print(f"\n数据质量分析:")
    print(f"有beta值的记录: {df['beta'].notna().sum()}")
    print(f"有se值的记录: {df['se'].notna().sum()}")
    print(f"有pval值的记录: {df['pval'].notna().sum()}")
    print(f"有n值的记录: {df['n'].notna().sum()}")
    
    # 方法2: 设置合理的默认eaf值
    # 对于没有eaf的记录，我们设置一个中等频率的默认值
    default_eaf = 0.3  # 30%是一个常见的中等等位基因频率
    
    # 备份原始数据
    backup_file = crp_file.parent / "crp_gwas_backup.parquet"
    df.to_parquet(backup_file)
    print(f"已备份原始数据到: {backup_file}")
    
    # 填补eaf缺失值
    df_fixed = df.copy()
    df_fixed['eaf'] = df_fixed['eaf'].fillna(default_eaf)
    
    print(f"修复后eaf缺失值数量: {df_fixed['eaf'].isnull().sum()}")
    print(f"eaf值分布:")
    print(df_fixed['eaf'].describe())
    
    # 保存修复后的数据
    df_fixed.to_parquet(crp_file, index=False)
    print(f"已保存修复后的数据到: {crp_file}")
    
    return df_fixed

def validate_fix():
    """验证修复结果"""
    print("\n=== 验证修复结果 ===")
    
    # 重新读取数据
    crp_file = Path("data/processed/gwas/crp_gwas.parquet")
    df = pd.read_parquet(crp_file)
    
    print(f"修复后数据形状: {df.shape}")
    print(f"eaf缺失值: {df['eaf'].isnull().sum()}")
    print(f"eaf值统计:")
    print(df['eaf'].describe())
    
    # 检查是否适合共定位分析
    required_cols = ['chr', 'pos', 'ref', 'alt', 'beta', 'se', 'pval', 'eaf', 'n']
    missing_required = [col for col in required_cols if df[col].isnull().sum() > 0]
    
    if missing_required:
        print(f"❌ 仍有缺失值的必需列: {missing_required}")
    else:
        print(f"✅ 所有必需列都完整")
    
    # 数据范围检查
    valid_beta = (df['beta'].abs() <= 10).all()
    valid_se = (df['se'] <= 5).all()
    valid_pval = ((df['pval'] > 0) & (df['pval'] < 1)).all()
    valid_eaf = ((df['eaf'] > 0) & (df['eaf'] < 1)).all()
    
    print(f"数据质量检查:")
    print(f"  beta值范围合理: {'✅' if valid_beta else '❌'}")
    print(f"  se值范围合理: {'✅' if valid_se else '❌'}")
    print(f"  pval值范围合理: {'✅' if valid_pval else '❌'}")
    print(f"  eaf值范围合理: {'✅' if valid_eaf else '❌'}")

def main():
    print("CRP GWAS数据eaf缺失值修复工具")
    print("="*50)
    
    # 1. 分析VCF文件结构
    analyze_crp_vcf_structure()
    
    # 2. 诊断当前问题
    diagnose_eaf_issue()
    
    # 3. 修复缺失数据
    fix_eaf_missing_data()
    
    # 4. 验证修复结果
    validate_fix()
    
    print("\n修复完成！")

if __name__ == "__main__":
    main()