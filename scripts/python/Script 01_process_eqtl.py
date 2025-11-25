# ==============================================================================
# Script: 01_process_eqtl.py (V2 - 修复 FUMA 基因坐标 Liftover)
# Purpose: 准备 eQTL 数据 (第1步: 加载 GTEx-v8 肺部 eQTLs 并与 FUMA 基因注释合并)
#
# 依赖:
# 1. Python 包: pandas, rich, pyarrow, pyliftover
#    (请使用 pip install pandas rich pyarrow pyliftover 安装)
#
# ==============================================================================

import logging
import sys
import pandas as pd
from pathlib import Path
from rich.logging import RichHandler
from rich.progress import track
from pyliftover import LiftOver  # <--- 新增: 导入 pyliftover

# ==============================================================================
# 配置参数 (与 02_prepare_data.py 保持一致)
# ==============================================================================

# 数据路径
BASE_DIR = Path(".")
RAW_DATA_DIR = BASE_DIR / "data/raw"
PROCESSED_DIR = BASE_DIR / "data/processed"
LOG_DIR = BASE_DIR / "logs"

# 确保目录存在
LOG_DIR.mkdir(parents=True, exist_ok=True)
(PROCESSED_DIR / "eqtl").mkdir(parents=True, exist_ok=True)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(str(LOG_DIR / "01_process_eqtl_py.log"), encoding='utf-8'),
        RichHandler()
    ]
)
log = logging.getLogger("rich")
log.info("Starting eQTL data preparation (Python high-speed version)")

# 1. FUMA 基因文件 (用于获取基因注释)
#    (我们复用 02_prepare_data.py 脚本中的 FUMA_GENES_FILES 列表)
FUMA_GENES_FILES = {
    "BMI": RAW_DATA_DIR / "BMI-genes.txt",
    "CRP": RAW_DATA_DIR / "CRP-genes.txt",
    "WBC": RAW_DATA_DIR / "WBC-genes.txt",
    "Alcohol": RAW_DATA_DIR / "Alcohol drink-genes.txt"
}

# 2. GTEx eQTL 源文件
EQTL_FILE = RAW_DATA_DIR / "Lung.v8.signif_variant_gene_pairs.txt.gz"

# 3. 输出文件 (供 03_run_coloc.R 使用)
OUTPUT_FILE = PROCESSED_DIR / "eqtl" / "lung_eqtl_annotated.parquet"


# ==============================================================================
# 函数定义
# ==============================================================================

def load_gene_annotation_map(fuma_file_dict: dict) -> pd.DataFrame:
    """
    从所有 FUMA 'genes.txt' 文件中加载一个统一的基因注释表。
    对 FUMA 基因坐标进行 hg19 -> hg38 Liftover 转换。
    """
    all_genes_list = []
    log.info("Loading gene annotations from FUMA files...")
    
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
            log.warning(f"FUMA file not found: {fuma_file}. Skipping for gene map.")
            continue
            
        try:
            fuma_genes = pd.read_csv(fuma_file, sep='\t')
        except Exception as e:
            log.error(f"Failed to read file {fuma_file.name}: {e}")
            continue
            
        required_cols = ['symbol', 'ensg', 'chr', 'start', 'end']
        if not all(col in fuma_genes.columns for col in required_cols):
            log.warning(f"FUMA file {fuma_file.name} missing required columns. Skipping.")
            continue
            
        # 基础处理
        fuma_genes_processed = (
            fuma_genes[required_cols]
                .rename(columns={
                    "symbol": "gene_symbol",
                    "ensg": "gene_id",
                    "chr": "chromosome_name",
                    "start": "start_position",
                    "end": "end_position"
                })
                .loc[fuma_genes['chr'].astype(str).isin([str(c) for c in range(1, 23)])]
                .astype({"chromosome_name": int, "start_position": int, "end_position": int})
                .drop_duplicates(subset=['gene_id'])
                .copy()
        )
        # 移除 gene_id 的版本号 (e.g., from 'ENSG000001.1' to 'ENSG000001')
        fuma_genes_processed['gene_id'] = fuma_genes_processed['gene_id'].str.split('.').str[0]
        
        # <--- 新增: 对 FUMA 基因坐标进行 Liftover 转换 (hg19 -> hg38)
        log.info(f"  Converting FUMA gene coordinates from hg19 to hg38 for {exposure_name}...")
        converted_genes = []
        skipped_liftover = 0
        
        for idx, row in fuma_genes_processed.iterrows():
            original_chr = f"chr{row['chromosome_name']}"
            original_start = row['start_position']
            original_end = row['end_position']
            
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
                row_dict['chromosome_name'] = int(new_chr)
                row_dict['start_position'] = int(new_start)
                row_dict['end_position'] = int(new_end)
                
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
        
        all_genes_list.append(fuma_genes_processed)
        
    if not all_genes_list:
        log.error("No gene annotation files (FUMA) found. Cannot proceed.")
        return None
        
    # 合并并去重
    all_genes_df = pd.concat(all_genes_list).drop_duplicates(subset=['gene_id']).reset_index(drop=True)
    log.info(f"Loaded {all_genes_df['gene_id'].nunique()} unique gene annotations (all in hg38).")
    return all_genes_df

def process_eqtl_data(eqtl_file_path: Path, gene_map: pd.DataFrame) -> pd.DataFrame:
    """
    读取 GTEx eQTL 文件, 解析, 并与基因图谱合并。
    """
    if not eqtl_file_path.exists():
        log.error(f"eQTL file not found: {eqtl_file_path}")
        return None
        
    log.info(f"Reading eQTL data from {eqtl_file_path.name}...")
    try:
        # 读取所有列，包括MAF (修复eQTL MAF问题)
        df = pd.read_csv(
            eqtl_file_path,
            sep='\t',
            compression='gzip'
        )
        log.info(f"  Loaded {len(df)} significant pairs with {len(df.columns)} columns")
        
        # 记录MAF列信息
        if 'maf' in df.columns:
            maf_stats = df['maf'].describe()
            log.info(f"  Found MAF column: range {maf_stats['min']:.4f}-{maf_stats['max']:.4f}, mean {maf_stats['mean']:.4f}")
        else:
            log.warning(f"  MAF column not found in eQTL data")
            
    except Exception as e:
        log.error(f"Failed to read eQTL file: {e}")
        return None
        
    log.info(f"  Loaded {len(df)} significant pairs.")
    
    # 移除 gene_id 的版本号以进行匹配
    df['gene_id'] = df['gene_id'].str.split('.').str[0]
    
    # --- 1. 解析 Variant ID ---
    # GTEx variant_id 格式为: 'chr1_12345_C_T_b38'
    log.info("  Parsing 'variant_id' column to get 'chr' and 'pos'...")
    try:
        var_parts = df['variant_id'].str.split('_', expand=True)
        df['chr'] = var_parts[0].str.replace('chr', '')
        df['pos'] = var_parts[1]
    except Exception as e:
        log.error(f"Failed to parse 'variant_id'. Format must be 'chr_pos_ref_alt_build'. Error: {e}")
        return None

    # --- 2. 清理和重命名 ---
    # 过滤到 1-22 号染色体
    df = df.loc[df['chr'].isin([str(c) for c in range(1, 23)])].copy()
    df['chr'] = df['chr'].astype(int)
    df['pos'] = df['pos'].astype(int)
    
    # 重命名列以匹配 03_run_coloc.R 脚本的期望
    df = df.rename(columns={
        "pval_nominal": "pval",
        "slope": "beta",
        "slope_se": "se"
    })
    
    # <--- 新增: 保留MAF列 (修复eQTL MAF问题)
    if 'maf' in df.columns:
        log.info("  Preserving MAF column for colocalization analysis")
        # MAF列已经存在，不需要重命名
    else:
        log.warning("  MAF column not found, will use default values in R script")
    
    # --- 3. 合并基因注释 ---
    # (注意: GTEx 'gene_id' 和 FUMA 'ensg' 都包含版本号, e.g., ENSG000001.1, 所以它们可以直接合并)
    log.info(f"  Merging {len(df)} eQTLs with {len(gene_map)} gene annotations...")
    
    df_merged = df.merge(gene_map, on='gene_id', how='left')
    
    # 检查合并是否成功
    n_unmapped = df_merged['gene_symbol'].isnull().sum()
    if n_unmapped > 0:
        log.warning(f"  {n_unmapped} eQTL records (out of {len(df_merged)}) could not be mapped to a gene annotation. They will be dropped.")
    
    # 丢弃那些没有在 FUMA 列表中找到注释的 eQTLs
    df_final = df_merged.dropna(subset=['gene_symbol', 'chromosome_name', 'start_position', 'end_position']).copy()
    
    log.info(f"  Created final annotated dataset with {len(df_final)} records.")
    
    return df_final

# ==============================================================================
# 主程序
# ==============================================================================

def main():
    
    # --- STEP 1: 加载基因注释 (来自 FUMA) ---
    gene_annotation_map = load_gene_annotation_map(FUMA_GENES_FILES)
    
    if gene_annotation_map is None:
        log.critical("Could not load any gene annotations from FUMA files. Stopping.")
        sys.exit(1)

    # --- STEP 2: 处理 eQTL 数据 (来自 GTEx) ---
    eqtl_data = process_eqtl_data(EQTL_FILE, gene_annotation_map)
    
    if eqtl_data is None:
        log.critical("Failed to process eQTL data. Stopping.")
        sys.exit(1)
        
    # --- STEP 3: 保存为 Parquet (供 R 脚本使用) ---
    eqtl_data.to_parquet(OUTPUT_FILE, index=False)
    log.info(f"✓ Successfully saved annotated eQTL data to {OUTPUT_FILE}")
    
    log.info("\nNext step: Run 02_prepare_data.py")

# ==============================================================================
# 运行
# ==============================================================================

if __name__ == "__main__":
    main()