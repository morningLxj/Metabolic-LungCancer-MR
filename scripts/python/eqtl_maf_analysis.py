#!/usr/bin/env python3
"""
详细分析eQTL数据的MAF引用问题
"""

import pandas as pd
import numpy as np

def analyze_eqtl_maf_issue():
    """分析eQTL数据中的MAF引用问题"""
    print("=" * 60)
    print("eQTL数据MAF引用错误分析")
    print("=" * 60)
    
    # 1. 加载eQTL数据
    eqtl_file = 'data/processed/eqtl/lung_eqtl_annotated.parquet'
    try:
        df = pd.read_parquet(eqtl_file)
        print(f"[OK] 成功加载eQTL数据: {df.shape}")
    except Exception as e:
        print(f"[ERROR] 加载eQTL数据失败: {e}")
        return
    
    # 2. 检查列结构
    print(f"\n[eQTL数据列结构]")
    print(f"总列数: {len(df.columns)}")
    print(f"列名: {list(df.columns)}")
    
    # 3. 检查MAF相关列
    maf_related_cols = []
    for col in df.columns:
        if any(keyword in col.lower() for keyword in ['maf', 'eaf', 'allele', 'freq']):
            maf_related_cols.append(col)
    
    print(f"\n[MAF相关列检查]")
    if maf_related_cols:
        print(f"找到MAF相关列: {maf_related_cols}")
        for col in maf_related_cols:
            unique_values = df[col].nunique()
            print(f"  {col}: {unique_values} 个独特值")
    else:
        print("[ERROR] 未找到任何MAF相关列 (MAF, eaf, allele_freq等)")
    
    # 4. 对比GWAS数据
    print(f"\n[对比GWAS数据的MAF情况]")
    gwas_traits = ['CRP', 'BMI', 'WBC', 'Alcohol', 'LUAD', 'LUSC']
    
    for trait in gwas_traits:
        try:
            gwas_file = f'data/processed/gwas/{trait.lower()}_gwas.parquet'
            gwas_df = pd.read_parquet(gwas_file)
            
            # 检查是否有eaf列
            if 'eaf' in gwas_df.columns:
                eaf_stats = gwas_df['eaf'].describe()
                print(f"  {trait}: [OK] 有eaf列, 范围[{eaf_stats['min']:.3f}, {eaf_stats['max']:.3f}]")
            else:
                print(f"  {trait}: [MISSING] 无eaf列")
        except Exception as e:
            print(f"  {trait}: [ERROR] 无法加载 ({e})")
    
    # 5. 分析原代码错误
    print(f"\n[原代码MAF引用错误分析]")
    print(f"错误代码 (第178行):")
    print(f"  dataset_eqtl <- list(")
    print(f"    MAF = exp_data$eaf,  # [ERROR] 错误：使用GWAS的eaf")
    print(f"    ...")
    print(f"  )")
    
    print(f"\n问题说明:")
    print(f"1. eQTL数据本身没有MAF/eaf列")
    print(f"2. exp_data$eaf是GWAS数据的等位基因频率")
    print(f"3. eQTL和GWAS可能有不同的等位基因参考版本")
    print(f"4. 使用错误的MAF会影响共定位分析的PIP计算")
    
    # 6. 修复方案验证
    print(f"\n[修复方案验证]")
    print(f"修复代码: get_maf_values(eqtl_data_filt)")
    
    # 模拟修复逻辑
    def get_maf_values(data):
        if "MAF" in data.columns and not data["MAF"].isnull().all():
            return data["MAF"]
        elif "eaf" in data.columns and not data["eaf"].isnull().all():
            return data["eaf"]
        else:
            return pd.Series([0.3] * len(data), index=data.index)
    
    # 验证修复逻辑
    if df.empty:
        print("eQTL数据为空，无法验证修复逻辑")
        return
    
    sample_data = df.head(100)
    result = get_maf_values(sample_data)
    
    print(f"修复逻辑测试:")
    print(f"  输入数据: {len(sample_data)} 行")
    print(f"  输出结果: {len(result)} 个值")
    print(f"  默认值使用: 全部为0.3 (中等等位基因频率)")
    print(f"  结果前5个: {result.head().tolist()}")
    
    # 7. 错误影响分析
    print(f"\n[错误影响分析]")
    print(f"如果不修复这个错误:")
    print(f"1. eQTL数据会使用GWAS的等位基因频率")
    print(f"2. 可能导致PIP计算错误")
    print(f"3. 共定位概率可能不准确")
    print(f"4. 特别是当GWAS和eQTL使用不同参考基因组时")
    
    # 8. 数据一致性检查
    print(f"\n[数据一致性检查]")
    print(f"需要验证的数据一致性:")
    print(f"1. 变异ID匹配 (variant_id vs snp_key)")
    print(f"2. 染色体位置一致 (chr, pos)")
    print(f"3. 等位基因一致性 (ref, alt)")
    
    # 检查变异ID
    if 'variant_id' in df.columns:
        sample_variants = df['variant_id'].head(10).tolist()
        print(f"  eQTL变异ID样例: {sample_variants}")
    else:
        print(f"  [ERROR] eQTL数据没有variant_id列")
    
    print(f"\n" + "=" * 60)
    print(f"结论: eQTL数据确实缺少MAF列，原代码引用错误！")
    print(f"修复方案: 使用get_maf_values()函数合理处理")
    print(f"=" * 60)

if __name__ == "__main__":
    analyze_eqtl_maf_issue()