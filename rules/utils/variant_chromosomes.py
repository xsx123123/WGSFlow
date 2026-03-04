#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pysam
from loguru import logger

# ----- 使用 pysam 自动获取染色体列表 ----- #
def get_chromosomes(config: dict = None, reference_fasta_parameter=None) -> list:
    """
    从嵌套配置中提取 Fasta 路径并获取染色体列表
    
    Args:
        config: Snakemake 的 config 字典
        reference_fasta_parameter: 字符串（单层）或 列表（多层路径，如 ["bwa_mem", "Lsat_v11", "genome_fa"]）
    """
    reference_fasta = "Unknown Path"
    try:
        if isinstance(reference_fasta_parameter, list):
            temp_config = config
            for key in reference_fasta_parameter:
                temp_config = temp_config[key]
            reference_fasta = temp_config
        else:
            reference_fasta = config[reference_fasta_parameter]
        
        if not reference_fasta or not os.path.exists(reference_fasta):
            raise FileNotFoundError(f"Fasta 路径无效或文件不存在: {reference_fasta}")

        with pysam.FastaFile(reference_fasta) as fa:
            chromosomes = list(fa.references)
            logger.debug(f"Successfully loaded {len(chromosomes)} sequences from {reference_fasta}")
        
        return chromosomes
    
    except KeyError as e:
        logger.error(f"配置文件中找不到指定的键: {e}")
        raise
    except Exception as e:
        logger.error(f"读取参考基因组 '{reference_fasta}' 出错: {e}")
        raise
