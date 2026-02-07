#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
import pysam
from loguru import logger
# ----- 使用 pysam 自动获取染色体列表 ----- #
def get_chromosomes(config):
    """
    从配置文件中读取参考基因组路径，使用 pysam 提取染色体列表
    
    Args:
        config: Snakemake 配置字典
    
    Returns:
        list: 染色体名称列表
    """
    try:
        reference_fasta = config["bcftools"]["reference"]
        with pysam.FastaFile(reference_fasta) as fa:
            chromosomes = list(fa.references)
            logger.debug(f"Found {fa.nreferences} chromosomes: {chromosomes}")
        
        return chromosomes
    
    except KeyError as e:
        logger.error(f"Missing configuration key: {e}")
        raise
    except Exception as e:
        logger.error(f"Error reading reference genome '{reference_fasta}': {e}")
        raise

CHROMOSOMES = get_chromosomes(config)
