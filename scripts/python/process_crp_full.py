#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
完整处理CRP文件（使用修复后的逻辑）
"""

import logging
import sys
import pandas as pd
import vcfpy
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from rich.logging import RichHandler
from rich.progress import track
from pyliftover import LiftOver

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[RichHandler()]
)
log = logging.getLogger("rich")

def process_crp_full():
    """完整处理CRP文件"""
    
    trait = "CRP"
    vcf_file = Path("data/raw/gwas/ebi-a-GCST90029070.vcf.gz")
    parquet_file = Path("data/processed/gwas/crp_gwas.parquet")
    
    if not vcf_file.exists():
        log.error(f"VCF file not found: {vcf_file}")
        return False
    
    log.info(f"Processing full CRP file: {vcf_file.name}")
    
    # 初始化 Liftover
    log.info("  Initializing liftover from hg19 to hg38...")
    try:
        lo = LiftOver('hg19', 'hg38')
    except Exception as e:
        log.error(f"  Failed to initialize LiftOver: {e}")
        return False
    
    # 打开 VCF 读取器
    try:
        reader = vcfpy.Reader.from_path(str(vcf_file))
    except Exception as e:
        log.error(f"  Failed to read VCF file: {e}")
        return False
    
    schema = None
    writer = None
    batch = []
    batch_size = 500_000
    count = 0
    skipped_info = 0
    skipped_liftover = 0
    
    # 遍历 VCF 记录
    try:
        for record in track(reader, description=f"  Reading {trait}...", total=None):
            count += 1
            
            # Liftover 逻辑
            original_chr = record.CHROM.replace("chr", "")
            original_pos = record.POS
            
            try:
                liftover_chr = f"chr{original_chr}"
                new_coords = lo.convert_coordinate(liftover_chr, original_pos)
                
                if not new_coords or len(new_coords) == 0:
                    skipped_liftover += 1
                    continue
                
                new_chr_str = new_coords[0][0]
                new_pos = new_coords[0][1]
                
                chr_val = new_chr_str.replace("chr", "")
                pos = new_pos
                
                if chr_val not in [str(c) for c in range(1, 23)]:
                    skipped_liftover += 1
                    continue
                
                chr_val = int(chr_val)
                
            except Exception as e:
                skipped_liftover += 1
                continue

            ref = record.REF
            alt = record.ALT[0].value
            rsid = record.ID[0] if record.ID else f"{chr_val}:{pos}"

            # 提取 FORMAT 字段
            try:
                if not hasattr(record, 'calls') or not record.calls:
                    skipped_info += 1
                    continue
                    
                call_data = record.calls[0].data
                
                # 从 CALL DATA 中获取 FORMAT 字段数据
                eaf = call_data.get('AF', [])
                beta = call_data.get('ES', [])
                se = call_data.get('SE', [])
                pval = call_data.get('LP', [])
                n = call_data.get('SS', [])
                
                # 处理列表格式
                if isinstance(eaf, (list, tuple)) and len(eaf) > 0:
                    eaf = eaf[0]
                else:
                    eaf = None
                    
                if isinstance(beta, (list, tuple)) and len(beta) > 0:
                    beta = beta[0]
                else:
                    beta = None
                    
                if isinstance(se, (list, tuple)) and len(se) > 0:
                    se = se[0]
                else:
                    se = None
                    
                if isinstance(pval, (list, tuple)) and len(pval) > 0:
                    pval = pval[0]
                else:
                    pval = None
                    
                if isinstance(n, (list, tuple)) and len(n) > 0:
                    n = n[0]
                else:
                    n = None
                
                # 检查必需字段
                if beta is None or se is None or pval is None:
                    skipped_info += 1
                    continue
                
                # 特殊处理CRP文件的缺失样本量
                if n is None:
                    if trait == "CRP":
                        n = 100000.0  # 10万样本量的默认值
                    else:
                        skipped_info += 1
                        continue
                    
                # 转换为数值类型
                try:
                    if eaf is not None:
                        eaf = float(eaf)
                    beta = float(beta)
                    se = float(se)
                    pval = float(pval)
                    n = float(n)
                except (ValueError, TypeError):
                    skipped_info += 1
                    continue
                    
                # 从 log-pvalue 转换回 p-value
                if pval is not None:
                    pval = 10 ** -pval
                
            except (KeyError, IndexError, TypeError, AttributeError):
                skipped_info += 1
                continue

            # 添加到批处理
            batch.append({
                "chr": chr_val,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "rsid": rsid,
                "eaf": eaf,
                "beta": beta,
                "se": se,
                "pval": pval,
                "n": n
            })

            # 写入批处理
            if len(batch) == batch_size:
                df = pd.DataFrame(batch)
                table = pa.Table.from_pandas(df)
                if writer is None:
                    schema = table.schema
                    writer = pq.ParquetWriter(parquet_file, schema)
                writer.write_table(table)
                batch = []
                log.info(f"  Processed {count:,} records...")

    finally:
        reader.close()
        # 写入最后一个批处理
        if batch:
            df = pd.DataFrame(batch)
            table = pa.Table.from_pandas(df)
            if writer is None:
                schema = table.schema
                writer = pq.ParquetWriter(parquet_file, schema)
            writer.write_table(table)
        
        if writer:
            writer.close()

    total_processed = count - skipped_info - skipped_liftover
    log.info(f"  Finished {trait}:")
    log.info(f"    Total variants read: {count:,}")
    log.info(f"    Skipped (missing info): {skipped_info:,}")
    log.info(f"    Skipped (Liftover fail/non-autosomal): {skipped_liftover:,}")
    log.info(f"    Total variants saved (hg38): {total_processed:,}")
    log.info(f"  Saved to: {parquet_file}")
    
    return total_processed > 0

if __name__ == "__main__":
    success = process_crp_full()
    sys.exit(0 if success else 1)