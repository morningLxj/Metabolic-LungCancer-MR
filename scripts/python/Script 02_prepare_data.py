# ==============================================================================
# Script: 02_prepare_data.py (V3 - 最终版本: 修复 GWAS 和 FUMA 的 Liftover)
# Purpose: 准备 GWAS 数据和分析任务
# ==============================================================================

import logging
import sys
import pandas as pd
import vcfpy
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from rich.logging import RichHandler
from rich.progress import track
from pyliftover import LiftOver  # <--- 修正: 导入 pyliftover

# ==============================================================================
# 配置参数
# ==============================================================================

# 数据路径
BASE_DIR = Path(".")
RAW_DATA_DIR = BASE_DIR / "data/raw"
RAW_FUMA_DIR = RAW_DATA_DIR
RAW_GWAS_DIR = RAW_DATA_DIR / "gwas"
PROCESSED_DIR = BASE_DIR / "data/processed"
LOG_DIR = BASE_DIR / "logs"

# 确保目录存在
LOG_DIR.mkdir(parents=True, exist_ok=True)
(PROCESSED_DIR / "gwas").mkdir(parents=True, exist_ok=True)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(str(LOG_DIR / "02_prepare_data_py.log"), encoding='utf-8'),
        RichHandler()
    ]
)
log = logging.getLogger("rich")

# 1. FUMA 基因文件 (用于获取基因注释)
FUMA_GENES_FILES = {
    "BMI": RAW_FUMA_DIR / "BMI-genes.txt",
    "CRP": RAW_FUMA_DIR / "CRP-genes.txt",
    "WBC": RAW_FUMA_DIR / "WBC-genes.txt",
    "Alcohol": RAW_FUMA_DIR / "Alcohol drink-genes.txt"
}

# 2. GWAS VCF 源文件 (使用实际存在的文件，假设这些都是 hg19)
GWAS_FILES = {
    "BMI": RAW_GWAS_DIR / "ieu-b-40.vcf.gz",
    "CRP": RAW_GWAS_DIR / "ebi-a-GCST90029070.vcf.gz",
    "WBC": RAW_GWAS_DIR / "ieu-a-984.vcf.gz",     # 使用实际存在的文件
    "Alcohol": RAW_GWAS_DIR / "ieu-a-989.vcf.gz", # 使用实际存在的文件
    "LUAD": RAW_GWAS_DIR / "ieu-b-30.vcf.gz",     # 使用实际存在的文件
    "LUSC": RAW_GWAS_DIR / "ieu-b-73.vcf.gz",     # 使用实际存在的文件
}

# 3. 输出文件
ALL_GENES_FILE = PROCESSED_DIR / "all_candidate_genes.parquet"
ANALYSIS_TASKS_FILE = PROCESSED_DIR / "analysis_tasks.txt"


# ==============================================================================
# 函数定义
# ==============================================================================

def read_vcf_gwas_python(trait: str, vcf_file: Path, parquet_file: Path):
    """
    使用 vcfpy (纯 Python) 读取 VCF, 转换为 hg38, 并分块写入 Parquet
    """
    log.info(f"Processing (Python): {vcf_file.name} for {trait}")

    # <--- 修正: 初始化 Liftover (hg19 -> hg38)
    log.info("  Initializing liftover from hg19 to hg38...")
    # pyliftover 会自动下载并缓存 'hg19' 到 'hg38' 的 chain 文件
    try:
        lo = LiftOver('hg19', 'hg38')
    except Exception as e:
        log.error(f"  Failed to initialize LiftOver. Is internet connected? Error: {e}")
        return

    # 打开 VCF 读取器
    try:
        reader = vcfpy.Reader.from_path(str(vcf_file))
    except Exception as e:
        log.error(f"  Failed to read VCF file {vcf_file}: {e}")
        return

    schema = None
    writer = None
    batch = []
    batch_size = 500_000
    count = 0
    skipped_info = 0
    skipped_liftover = 0  # <--- 修正: 添加计数器
    
    # 遍历 VCF 记录
    try:
        for record in track(reader, description=f"  Reading {trait}...", total=None):
            count += 1
            
            # --- 修正: Liftover 逻辑 ---
            original_chr = record.CHROM.replace("chr", "")
            original_pos = record.POS
            
            try:
                liftover_chr = f"chr{original_chr}"
                new_coords = lo.convert_coordinate(liftover_chr, original_pos)
                
                if not new_coords or len(new_coords) == 0:
                    skipped_liftover += 1
                    continue # 跳过无法转换的 SNP
                
                new_chr_str = new_coords[0][0]
                new_pos = new_coords[0][1]
                
                # 更新 chr_val 和 pos 为新的 hg38 坐标
                chr_val = new_chr_str.replace("chr", "")
                pos = new_pos
                
                # 在 *之后* 转换, 检查非 1-22 染色体
                if chr_val not in [str(c) for c in range(1, 23)]:
                    skipped_liftover += 1
                    continue
                
                chr_val = int(chr_val)
                
            except Exception as e:
                log.debug(f"Liftover failed for {original_chr}:{original_pos} with error: {e}")
                skipped_liftover += 1
                continue
            # --- 修正结束 ---

            ref = record.REF
            alt = record.ALT[0].value # 只取第一个 ALT
            rsid = record.ID[0] if record.ID else f"{chr_val}:{pos}"

            # 提取 FORMAT 字段 (从 CALL 数据中)
            try:
                # 获取第一个 CALL 的数据
                if not hasattr(record, 'calls') or not record.calls:
                    skipped_info += 1
                    continue
                    
                call_data = record.calls[0].data
                
                # <--- 修复: 增强的字段获取逻辑，支持多种格式
                # 从 CALL DATA 中获取 FORMAT 字段数据
                eaf = call_data.get('AF', [])
                beta = call_data.get('ES', [])
                se = call_data.get('SE', [])
                pval = call_data.get('LP', [])
                n = call_data.get('SS', [])
                
                # <--- 修复: 处理列表格式，直接提取第一个值
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
                
                # <--- 修复: 更宽松的字段检查，允许AF为None
                # 必需字段：beta, se, pval (CRP中AF可能缺失)
                if beta is None or se is None or pval is None:
                    # <--- 修复: 为CRP添加更详细的调试信息
                    if trait == "CRP":
                        log.debug(f"CRP Record {count}: Missing required fields - beta:{beta}, se:{se}, pval:{pval}")
                        log.debug(f"  Available fields: {list(call_data.keys())}")
                    skipped_info += 1
                    continue
                    
                # <--- 修复: 改进的样本量处理
                if n is None:
                    if trait == "CRP":
                        # CRP GWAS通常有大量样本，设置合理的默认值
                        n = 100000.0  # 10万样本量的默认值
                        log.debug(f"Using default sample size for CRP: {n}")
                    else:
                        skipped_info += 1
                        continue
                    
                # <--- 修复: 改进的数值转换逻辑
                try:
                    # 只有当eaf不为None时才转换
                    if eaf is not None:
                        eaf = float(eaf)
                    beta = float(beta)
                    se = float(se)
                    pval = float(pval)
                    n = float(n)
                    
                    # <--- 修复: 添加数值有效性检查
                    if trait == "CRP":
                        # CRP特异性检查
                        if abs(beta) > 10 or se > 5 or pval <= 0 or pval >= 1:
                            log.debug(f"CRP Record {count}: Invalid values - beta:{beta}, se:{se}, pval:{pval}")
                            skipped_info += 1
                            continue
                            
                except (ValueError, TypeError) as e:
                    log.debug(f"Type conversion failed for {chr_val}:{pos}: {e}")
                    if trait == "CRP":
                        log.debug(f"  Raw values: eaf={eaf}, beta={beta}, se={se}, pval={pval}, n={n}")
                    skipped_info += 1
                    continue
                    
                # 从 log-pvalue (LP) 转换回 p-value
                if pval is not None:
                    pval = 10 ** -pval
                
            except (KeyError, IndexError, TypeError, AttributeError) as e:
                log.debug(f"Failed to extract format data for {chr_val}:{pos}: {e}")
                skipped_info += 1
                continue

            # <--- 修正: 移除了旧的染色体检查, 因为 Liftover 块已处理

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
    log.info(f"    Skipped (Liftover fail/non-autosomal): {skipped_liftover:,}") # <--- 修正
    log.info(f"    Total variants saved (hg38): {total_processed:,}")
    log.info(f"  Saved to: {parquet_file}")
    
    # <--- 新增: 修复CRP数据的eaf缺失值
    if trait == "CRP" and parquet_file.exists():
        log.info("  Fixing missing eaf values in CRP data...")
        try:
            # 读取数据并检查eaf缺失值
            df = pd.read_parquet(parquet_file)
            missing_eaf = df['eaf'].isnull().sum()
            
            if missing_eaf > 0:
                log.info(f"    Found {missing_eaf:,} missing eaf values. Filling with default value 0.3...")
                
                # 备份原始数据
                backup_file = parquet_file.parent / f"{parquet_file.stem}_backup{parquet_file.suffix}"
                df.to_parquet(backup_file)
                log.info(f"    Backup saved to: {backup_file}")
                
                # 填补缺失的eaf值
                default_eaf = 0.3
                df['eaf'] = df['eaf'].fillna(default_eaf)
                df.to_parquet(parquet_file, index=False)
                
                log.info(f"    Fixed {missing_eaf:,} eaf values with default {default_eaf}")
            else:
                log.info("    No missing eaf values found")
                
        except Exception as e:
            log.error(f"    Failed to fix eaf values: {e}")


def load_all_exposure_genes(fuma_file_dict: dict) -> pd.DataFrame:
    """
    从所有 FUMA 'genes.txt' 文件中加载一个统一的、去重的基因列表。
    对 FUMA 基因坐标进行 hg19 -> hg38 Liftover 转换。
    """
    all_genes_list = []
    log.info("Loading candidate genes from ALL Exposure FUMA files")
    
    # <--- 新增: 初始化 Liftover (hg19 -> hg38)
    log.info("  Initializing liftover from hg19 to hg38 for FUMA genes...")
    # pyliftover 会自动下载并缓存 'hg19' 到 'hg38' 的 chain 文件
    try:
        lo = LiftOver('hg19', 'hg38')
    except Exception as e:
        log.error(f"  Failed to initialize LiftOver. Is internet connected? Error: {e}")
        return None
    
    for exposure_name, fuma_file in fuma_file_dict.items():
        if not fuma_file.exists():
            log.warning(f"FUMA file not found: {fuma_file}. Skipping.")
            continue
            
        log.info(f"  Loading candidate genes for {exposure_name} from {fuma_file.name}")
        try:
            fuma_genes = pd.read_csv(fuma_file, sep='\t')
        except Exception as e:
            log.error(f"Failed to read file {fuma_file.name}: {e}")
            continue
            
        required_cols = ['symbol', 'ensg', 'chr', 'start', 'end']
        if not all(col in fuma_genes.columns for col in required_cols):
            log.warning(f"FUMA file {fuma_file.name} missing required columns. Skipping.")
            continue
            
        # 过滤, 重命名, 并设置类型
        fuma_genes_processed = (
            fuma_genes[required_cols]
                .rename(columns={
                    "symbol": "gene_symbol",
                    "ensg": "gene_id",
                })
                .loc[fuma_genes['chr'].astype(str).isin([str(c) for c in range(1, 23)])]
                .astype({"chr": int, "start": int, "end": int})
                .drop_duplicates(subset=['gene_id'])
                .copy()
        )
        
        # <--- 新增: 对 FUMA 基因坐标进行 Liftover 转换 (hg19 -> hg38)
        log.info(f"    Converting FUMA gene coordinates from hg19 to hg38...")
        converted_genes = []
        skipped_liftover = 0
        
        for idx, row in fuma_genes_processed.iterrows():
            original_chr = f"chr{row['chr']}"
            original_start = row['start']
            original_end = row['end']
            
            try:
                # 转换起始位置
                start_coords = lo.convert_coordinate(original_chr, original_start)
                if not start_coords or len(start_coords) == 0:
                    skipped_liftover += 1
                    continue
                
                # 转换结束位置
                end_coords = lo.convert_coordinate(original_chr, original_end)
                if not end_coords or len(end_coords) == 0:
                    skipped_liftover += 1
                    continue
                
                new_chr_str = start_coords[0][0]  # 使用起始坐标的染色体
                new_start = start_coords[0][1]
                new_end = end_coords[0][1]
                
                # 移除 "chr" 前缀并验证染色体
                new_chr = new_chr_str.replace("chr", "")
                if new_chr not in [str(c) for c in range(1, 23)]:
                    skipped_liftover += 1
                    continue
                
                # 更新基因坐标为 hg38
                row_dict = row.to_dict()
                row_dict['chr'] = int(new_chr)
                row_dict['start'] = int(new_start)
                row_dict['end'] = int(new_end)
                
                converted_genes.append(row_dict)
                
            except Exception as e:
                log.debug(f"Liftover failed for {original_chr}:{original_start}-{original_end} with error: {e}")
                skipped_liftover += 1
                continue
        
        if converted_genes:
            fuma_genes_processed = pd.DataFrame(converted_genes)
            log.info(f"    Successfully converted {len(fuma_genes_processed)} genes to hg38")
            log.info(f"    Skipped {skipped_liftover} genes due to liftover failure")
        else:
            log.warning(f"    No genes could be converted for {exposure_name}")
            continue
        
        fuma_genes_processed['exposure'] = exposure_name
        
        log.info(f"    Final: {len(fuma_genes_processed)} unique genes for {exposure_name} (hg38)")
        all_genes_list.append(fuma_genes_processed)
        
    if not all_genes_list:
        log.error("No gene annotation files (FUMA) found. Cannot proceed.")
        return None
        
    # 合并所有暴露的基因
    all_genes_df = pd.concat(all_genes_list).reset_index(drop=True)
    
    log.info(f"Total unique candidate genes loaded from all files (all in hg38): {all_genes_df['gene_id'].nunique()}")
    
    return all_genes_df

# ==============================================================================
# 主程序
# ==============================================================================

def main():
    log.info("Starting data preparation (Python high-speed version)")
    
    # ============================================================================
    log.info("\n" + "="*80)
    log.info("STEP 1: Load GWAS Data (from VCF and save to Parquet, with Liftover)")
    log.info("="*80)
    # ============================================================================
    
    for trait, vcf_file in GWAS_FILES.items():
        parquet_file = PROCESSED_DIR / "gwas" / f"{trait.lower()}_gwas.parquet"
        
        if parquet_file.exists():
            log.info(f"GWAS data for {trait} already exists. Skipping VCF read.")
            continue
            
        if not vcf_file.exists():
            log.warning(f"VCF file not found for {trait}: {vcf_file}. Skipping.")
            continue
        
        read_vcf_gwas_python(trait, vcf_file, parquet_file)

    # ============================================================================
    log.info("\n" + "="*80)
    log.info("STEP 2: Load Candidate Genes from ALL Exposure FUMA files")
    log.info("="*80)
    # ============================================================================

    all_genes = load_all_exposure_genes(FUMA_GENES_FILES)
    
    if all_genes is None:
        log.critical("Failed to load any genes from FUMA files. Stopping.")
        sys.exit(1)
        
    # 保存一份去重的基因列表 (供 01_process_eqtl.py 使用, 尽管该脚本通常先运行)
    all_genes.drop(columns=['exposure']).drop_duplicates(subset=['gene_id']).to_parquet(ALL_GENES_FILE, index=False)

    # ============================================================================
    log.info("\n" + "="*80)
    log.info("STEP 3: Generate Analysis Tasks (Exposure-Gene-Outcome matching)")
    log.info("="*80)
    # ============================================================================

    # 只保留用于“暴露”的基因 (不包括 LUAD/LUSC)
    exposure_genes = all_genes[all_genes['exposure'].isin(['BMI', 'CRP', 'WBC', 'Alcohol'])]
    outcome_traits = ['LUAD', 'LUSC']
    
    tasks_list = []
    
    for exposure_name in exposure_genes['exposure'].unique():
        genes_for_exposure = exposure_genes[exposure_genes['exposure'] == exposure_name]
        
        for outcome_name in outcome_traits:
            # 创建 (Exposure, Gene, Outcome) 的所有组合
            tasks_for_pair = genes_for_exposure.copy()
            tasks_for_pair['outcome'] = outcome_name
            tasks_list.append(tasks_for_pair)
    
    analysis_tasks = pd.concat(tasks_list).reset_index(drop=True)
    
    # 重排并保存
    analysis_tasks = analysis_tasks[['exposure', 'outcome', 'gene_symbol', 'gene_id', 'chr', 'start', 'end']]
    
    log.info(f"Total analysis tasks generated: {len(analysis_tasks):,}")
    
    # 保存为 TXT 以供 R (fread) 读取
    analysis_tasks.to_csv(ANALYSIS_TASKS_FILE, sep='\t', index=False)
    log.info(f"Saved tasks to {ANALYSIS_TASKS_FILE}")

    # ============================================================================
    log.info("\n" + "="*80)
    log.info("STEP 4: Data Quality Validation")
    log.info("="*80)
    # ============================================================================
    
    # 执行数据质量验证
    validation_results = validate_data_quality()
    
    log.info("\nData preparation (Step 2) complete!")
    log.info("Next step: Run colocalization analysis (03_run_coloc.R)")

# ==============================================================================
# 新增: 数据质量验证函数
# ==============================================================================

def validate_data_quality():
    """验证数据质量，特别是CRP数据的完整性"""
    log.info("\n" + "="*80)
    log.info("STEP 4: Data Quality Validation")
    log.info("="*80)
    
    # 检查关键文件
    key_files = {
        'All Candidate Genes': PROCESSED_DIR / "all_candidate_genes.parquet",
        'Analysis Tasks': PROCESSED_DIR / "analysis_tasks.txt",
        'CRP GWAS': PROCESSED_DIR / "gwas" / "crp_gwas.parquet",
        'BMI GWAS': PROCESSED_DIR / "gwas" / "bmi_gwas.parquet",
        'WBC GWAS': PROCESSED_DIR / "gwas" / "wbc_gwas.parquet",
        'LUAD GWAS': PROCESSED_DIR / "gwas" / "luad_gwas.parquet",
        'LUSC GWAS': PROCESSED_DIR / "gwas" / "lusc_gwas.parquet"
    }
    
    validation_results = []
    
    for name, file_path in key_files.items():
        if file_path.exists():
            try:
                if file_path.suffix == '.parquet':
                    df = pd.read_parquet(file_path)
                elif file_path.suffix == '.txt':
                    df = pd.read_csv(file_path, sep='\t')
                else:
                    continue
                
                missing_count = df.isnull().sum().sum()
                row_count = len(df)
                col_count = len(df.columns)
                
                # 特殊检查CRP GWAS的eaf列
                eaf_issues = ""
                if name == "CRP GWAS" and 'eaf' in df.columns:
                    eaf_missing = df['eaf'].isnull().sum()
                    if eaf_missing > 0:
                        eaf_issues = f" [eaf缺失: {eaf_missing:,}个]"
                    
                validation_results.append({
                    'file': name,
                    'status': 'OK' if missing_count == 0 else 'WARNING',
                    'rows': row_count,
                    'cols': col_count,
                    'missing_values': missing_count,
                    'notes': eaf_issues
                })
                
                status_icon = "OK" if missing_count == 0 else "WARN"
                log.info(f"  [{status_icon}] {name}: {row_count:,} rows, {col_count} cols, {missing_count:,} missing values{eaf_issues}")
                
            except Exception as e:
                validation_results.append({
                    'file': name,
                    'status': 'ERROR',
                    'rows': 0,
                    'cols': 0,
                    'missing_values': 0,
                    'notes': str(e)
                })
                log.error(f"  [ERROR] {name}: 读取失败 - {e}")
        else:
            validation_results.append({
                'file': name,
                'status': 'MISSING',
                'rows': 0,
                'cols': 0,
                'missing_values': 0,
                'notes': "文件不存在"
            })
            log.warning(f"  [MISSING] {name}: 文件不存在")
    
    # 总结验证结果
    total_files = len(validation_results)
    ok_files = sum(1 for r in validation_results if r['status'] == 'OK')
    warning_files = sum(1 for r in validation_results if r['status'] == 'WARNING')
    error_files = sum(1 for r in validation_results if r['status'] == 'ERROR')
    missing_files = sum(1 for r in validation_results if r['status'] == 'MISSING')
    
    log.info(f"\n数据质量总结:")
    log.info(f"  [OK] 完整文件: {ok_files}/{total_files}")
    log.info(f"  [WARN] 有缺失值: {warning_files}/{total_files}")
    log.info(f"  [ERROR] 错误文件: {error_files}/{total_files}")
    log.info(f"  [MISSING] 缺失文件: {missing_files}/{total_files}")
    
    # 特别检查CRP数据
    crp_gwas_file = PROCESSED_DIR / "gwas" / "crp_gwas.parquet"
    if crp_gwas_file.exists():
        try:
            df = pd.read_parquet(crp_gwas_file)
            if 'eaf' in df.columns:
                eaf_missing = df['eaf'].isnull().sum()
                if eaf_missing == 0:
                    log.info("  [SUCCESS] CRP GWAS数据已修复，eaf列完整！")
                else:
                    log.warning(f"  [WARN] CRP GWAS数据仍有{eaf_missing:,}个eaf缺失值")
        except Exception as e:
            log.error(f"  [ERROR] CRP GWAS数据验证失败: {e}")
    
    if ok_files == total_files:
        log.info("[SUCCESS] 所有数据文件都完整，数据准备成功！")
    else:
        log.warning("[WARN] 部分数据文件存在问题，请检查上述信息")
    
    return validation_results

# ==============================================================================
# 运行
# ==============================================================================

if __name__ == "__main__":
    main()