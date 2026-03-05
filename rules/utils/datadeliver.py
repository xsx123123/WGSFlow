#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
from pathlib import Path
from typing import Dict, List, Union
from snakemake.io import expand


def qc_clean(samples: Dict = None, data_deliver: List = None) -> List:
    """
    Handle quality control and cleaning steps.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    # fastq-screen
    data_deliver.extend(expand("01.qc/fastq_screen_r1/{sample}_R1_screen.txt", sample=samples.keys()))
    data_deliver.extend(expand("01.qc/fastq_screen_r2/{sample}_R2_screen.txt", sample=samples.keys()))
    data_deliver.append("01.qc/fastq_screen_multiqc_r1/multiqc_r1_fastq_screen_report.html")
    data_deliver.append("01.qc/fastq_screen_multiqc_r2/multiqc_r2_fastq_screen_report.html")
    # fastqc & multiqc
    data_deliver.append("01.qc/short_read_r1_multiqc/multiqc_r1_raw-data_report.html")
    data_deliver.append("01.qc/short_read_r2_multiqc/multiqc_r2_raw-data_report.html")
    # short-read trim & clean result
    data_deliver.append("01.qc/multiqc_short_read_trim/multiqc_short_read_trim_report.html")
    # merge qc
    data_deliver.append("01.qc/multiqc_merge_qc/multiqc_merge_qc_report.html")
    return data_deliver


def mapping(samples: Dict = None, data_deliver: List = None, config: Dict = None) -> List:
    """
    Handle sequence alignment and related outputs.
    
    This function separates 'Core Mapping' (BAM, CRAM, Basic Stats) from 
    'Deep QC' (RSeQC, TIN, BigWig), allowing for flexible control via config.
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []
    if config is None:
        config = {}

    # =========================================================================
    # [1] Core Mapping Outputs (基础核心输出)
    # 这些文件是后续所有分析的基础，无论 config 设置如何，都必须生成。
    # =========================================================================
    
    # BWA Alignment Results
    data_deliver.extend(expand('02.mapping/bwa_mem/{sample}/{sample}.raw.bam',sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/bwa_mem/sort_index/{sample}.sort.bam',sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/bwa_mem/sort_index/{sample}.sort.bam.bai',sample=samples.keys()))

    # MarkDuplicates Output (BWA MEM)
    data_deliver.extend(expand('02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',sample=samples.keys()))

    # BAM2CRAM
    data_deliver.extend(expand('02.mapping/cram/{sample}.Dup.cram',sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/cram/{sample}.Dup.cram.crai',sample=samples.keys()))

    # Mosdepth Coverage
    data_deliver.extend(expand("02.mapping/mosdepth_coverage/{sample}.sort.Dup.global.dist.txt",sample=samples.keys()))
    data_deliver.extend(expand("02.mapping/mosdepth_coverage/{sample}.sort.Dup.mosdepth.summary.txt",sample=samples.keys()))

    # qualimap qc
    data_deliver.extend(expand('02.mapping/qualimap_report/{sample}/qualimapReport.html',sample=samples.keys()))
    data_deliver.extend(expand('02.mapping/qualimap_report/{sample}/genome_results.txt',sample=samples.keys()))

    # samtools flagst
    data_deliver.extend(expand('02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv',sample=samples.keys()))
    
    # samtools stats
    data_deliver.extend(expand('02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv',sample=samples.keys()))

    # mapping_report
    data_deliver.append("02.mapping/mapping_report/mapping_report.html")

    return data_deliver

def bcftools_call_variant(samples: Dict = None, data_deliver: List = None) -> List:
    """
    Handle variant calling analysis.

    Args:
        samples: Dictionary of sample information
        data_deliver: List of deliverable files to extend

    Returns:
        Updated list of deliverable files
    """
    if samples is None:
        samples = {}
    if data_deliver is None:
        data_deliver = []

    # bcftools call variant
    data_deliver.extend(expand('03.call_variant/{sample}.sort.vcf.gz',sample=samples.keys()))
    data_deliver.extend(expand('03.call_variant/{sample}.sort.vcf.csi',sample=samples.keys()))
    data_deliver.extend(expand('03.call_variant/{sample}.sort.vcf.tbi',sample=samples.keys()))

    # merge & filter variant
    data_deliver.append('04.filter_variant/merge_filter.sort.vcf.gz')
    data_deliver.append('04.filter_variant/merge_filter.sort.vcf.gz.csi')
    data_deliver.append('04.filter_variant/merge_filter.sort.vcf.gz.tbi')

    # SnpEff Annotation Variant
    data_deliver.append('04.filter_variant/merge_filter.sort.annotation.csv')
    adata_deliver.append('04.filter_variant/merge_filter.sort.annotation.html')
    data_deliver.append('04.filter_variant/merge_filter.sort.annotation.vcf')

    # bcftools stats multiqc report
    data_deliver.append('03.call_variant/multiqc_variant_stats/multiqc_callvariant_report.html')
    
    return data_deliver


