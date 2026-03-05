#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - GATK Hard Filtering Module

This module implements a rule-based filtering approach (Hard Filtering) for 
joint-called variants. It is the preferred method for non-standard genomes 
(e.g., plants, non-model organisms) or small cohorts where the machine 
learning models of VQSR cannot be effectively trained.

Key Components:
- gatk_hard_filtering: Applies static thresholds to multiple quality metrics 
  to identify and flag low-confidence SNPs and Indels.

The filtering criteria follow GATK's recommended generic thresholds for 
QualityByDepth (QD), FisherStrand (FS), StrandOddsRatio (SOR), and 
MappingQuality (MQ).
"""

import os
from loguru import logger

# ----- rule ----- #

rule gatk_joint_hard_filtering:
    """
    Apply standard GATK Hard Filters to the cohort VCF.

    This rule uses GATK VariantFiltration to tag variants that do not meet 
    static quality thresholds.

    Default SNP Filters:
    - QD < 2.0: Low quality per depth
    - MQ < 40.0: Low mapping quality
    - FS > 60.0: High strand bias
    - SOR > 3.0: High strand odds ratio
    - MQRankSum < -12.5: Mapping quality bias
    - ReadPosRankSum < -8.0: Position bias within reads

    Default INDEL Filters:
    - QD < 2.0, FS > 200.0, SOR > 10.0, ReadPosRankSum < -20.0.
    """
    input:
        vcf = '03.call_variant/join_calling/cohort_raw.vcf.gz',
        tbi = '03.call_variant/join_calling/cohort_raw.vcf.gz.tbi',
    output:
        vcf = '04.filter_variant/cohort_filtered.vcf.gz',
        tbi = '04.filter_variant/cohort_filtered.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/03.call_variant/join_calling/gatk_hard_filtering.log",
    params:
        java_opts = "-Xmx16G",
        snp_filter = "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
        indel_filter = "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0",
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    threads:
        1
    shell:
        """
        gatk --java-options "{params.java_opts}" VariantFiltration \
            -R {params.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-expression "{params.snp_filter}" \
            --filter-name "HardFilteredSNP" \
            --filter-expression "{params.indel_filter}" \
            --filter-name "HardFilteredINDEL" &> {log}
        """
# ----- rule ----- #
