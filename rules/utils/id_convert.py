#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
from pathlib import Path
from rich import print as rprint
from typing import List, Dict, Optional

def _validate_df(df: pd.DataFrame, required_cols: List[str], index_col: str) -> None:
    """
    [内部函数] 校验 DataFrame 的完整性和唯一性
    """
    # Import the unified logger
    from snakemake_logger_plugin_rich_loguru import get_analysis_logger
    logger = get_analysis_logger()

    # 1. 校验必填列是否存在
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        error_msg = (
            f"❌ 样本表格式错误！\n"
            f"   缺失列: [bold red]{missing_cols}[/bold red]\n"
            f"   必需列: {required_cols}"
        )
        rprint(error_msg) # 使用 rich 打印高亮错误
        logger.error(f"Sample sheet missing columns: {missing_cols}")
        sys.exit(1)

    # 2. 校验索引列 (Sample ID) 是否有重复
    if df[index_col].duplicated().any():
        duplicated_ids = df[df[index_col].duplicated()][index_col].unique().tolist()
        error_msg = (
            f"❌ 样本ID不唯一！检测到重复样本名 (Sample ID):\n"
            f"   [bold red]{duplicated_ids}[/bold red]"
        )
        rprint(error_msg)
        logger.error(f"Duplicate sample IDs found: {duplicated_ids}")
        sys.exit(1)

    # 3. 校验是否有空值 (NaN)
    # 仅检查必填列中的空值
    if df[required_cols].isnull().any().any():
        nan_rows = df[df[required_cols].isnull().any(axis=1)][index_col].tolist()
        logger.warning(f"⚠️ 警告: 以下样本在必填列中存在空值 (NaN/Empty): {nan_rows}")
        # 这里可以选择是报错退出还是仅警告，目前设为警告


def load_samples(csv_path, required_cols=None, index_col="sample",
                 bam_path : str = "02.mapping/bwa_mem/MarkDup/"):
    """
    读取 CSV，并自动生成 WGS BWA-MEM 的 BAM 路径。
    """
    if required_cols is None:
        required_cols = [index_col, "group"]

    file_path = Path(csv_path)
    if not file_path.exists():
        print(f"❌ Error: 找不到样本表文件: {file_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(file_path, dtype=str, comment='#')
        df.columns = df.columns.str.strip()
        df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        _validate_df(df, required_cols, index_col)

        # 自动构建 BWA-MEM 路径 (与 mapping 规则对齐)
        df['bam'] = df[index_col].apply(
            lambda x: f"{bam_path}{x}.sort.Dup.bam"
        )
        df['bai'] = df[index_col].apply(
            lambda x: f"{bam_path}/{x}.sort.Dup.bam.bai"
        )

        return df.set_index(index_col, drop=False).to_dict(orient="index")

    except Exception as e:
        print(f"❌ Error: load_samples 解析失败: {e}", file=sys.stderr)
        sys.exit(1)


def load_contrasts(csv_path, samples_dict):
    """
    解析对比表 (Somatic Paired)。
    contrasts.csv 格式应为: Tumor, Normal
    """
    file_path = Path(csv_path)
    if not file_path.exists():
        print(f"❌ Error: 找不到对比表文件: {file_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(file_path, dtype=str, comment='#')
        df.columns = df.columns.str.strip()
        df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        if "Normal" not in df.columns or "Tumor" not in df.columns:
            print(f"❌ Error: contrasts.csv 必须包含 'Normal' 和 'Tumor' 列", file=sys.stderr)
            sys.exit(1)

        all_contrasts = []
        contrast_map = {}

        for _, row in df.iterrows():
            tumor_sample = row['Tumor']
            normal_sample = row['Normal']
            
            # 校验样本是否存在于样本表中
            if tumor_sample not in samples_dict or normal_sample not in samples_dict:
                print(f"⚠️ Warning: 样本 {tumor_sample} 或 {normal_sample} 不在样本表中，跳过该配对", file=sys.stderr)
                continue

            c_name = f"{tumor_sample}_vs_{normal_sample}"
            all_contrasts.append(c_name)
            
            # 存储配对的 BAM 路径和样本名 (SM tag)
            contrast_map[c_name] = {
                "tumor_bam": samples_dict[tumor_sample]['bam'],
                "normal_bam": samples_dict[normal_sample]['bam'],
                "tumor_name": tumor_sample,
                "normal_name": normal_sample
            }
            
        return all_contrasts, contrast_map

    except Exception as e:
        print(f"❌ Error: load_contrasts 解析失败: {e}", file=sys.stderr)
        sys.exit(1)

# 测试代码 (只有直接运行此脚本时才会执行)
if __name__ == "__main__":
    # 创建一个伪造的 csv 用于测试
    import io
    csv_content = """sample, sample_name, group, fq1
    s1, sample_1, control, s1_1.fq.gz
    s2, sample_2, treatment, s2_1.fq.gz
    """
    # 模拟读取
    print("--- 开始测试 ---")
    try:
        pass 
    except Exception as e:
        print(e)