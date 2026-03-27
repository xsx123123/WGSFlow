#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - GATK VQSR Filtering Module

This module implements Variant Quality Score Recalibration (VQSR), a sophisticated 
machine learning-based filtering method. VQSR is the GATK-recommended approach 
for high-quality genomes (e.g., Human, Mouse) where extensive "truth sets" 
and gold-standard variant databases are available.

Key Components:
- VariantRecalibrator: Trains a Gaussian Mixture Model on high-confidence sites 
  to distinguish between true variants and machine artifacts.
- ApplyVQSR: Uses the trained model to score variants and filter the VCF 
  based on target sensitivity tranches.

Workflow:
1. SNP Recalibration & Application
2. INDEL Recalibration & Application (on the SNP-filtered VCF)
"""

import os
from loguru import logger

# ----- rule ----- #

rule VariantRecalibrator_SNP:
    """
    Step 1: Build a high-dimensional recalibration model for SNPs.

    Trains a machine learning model using provided resource datasets (HapMap, 
    Omni, 1000G, dbSNP) to identify the multidimensional profile of true variants.

    Annotations (Metrics):
    - QD (QualByDepth), MQ (MappingQuality), MQRankSum, ReadPosRankSum, FS (FisherStrand), SOR.
    """
    input:
        vcf = '03.call_variant/join_calling/cohort_raw.vcf.gz',
        tbi = '03.call_variant/join_calling/cohort_raw.vcf.gz.tbi',
    output:
        recal = '04.filter_variant/VQSR_Filtering/cohort_snp.recal',
        tranches = '04.filter_variant/VQSR_Filtering/cohort_snp.tranches',
        plots = '04.filter_variant/VQSR_Filtering/cohort_snp_plots.R',
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/06.vqsr/VariantRecalibrator_SNP.log",
    params:
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        hapmap = config['bwa_mem'][config['Genome_Version']]['hapmap'],
        omni = config['bwa_mem'][config['Genome_Version']]['omni'],
        g1000 = config['bwa_mem'][config['Genome_Version']]['1000G_phase1'],
        dbsnp = config['bwa_mem'][config['Genome_Version']]['dbsnp138'],
        gatk_memory = config['parameter']['gatk_memory']["VariantRecalibrator_SNP"],
    threads:
        config['parameter']['threads']['gatk_vqsr'],
    shell:
        """
        gatk --java-options "{params.gatk_memory} -XX:ParallelGCThreads={threads}" VariantRecalibrator \
            -R {params.ref} \
            -V {input.vcf} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.g1000} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.plots} &> {log}
        """

rule ApplyVQSR_SNP:
    """
    Step 2: Apply the SNP recalibration model to the cohort VCF.

    Uses the .recal and .tranches files produced in Step 1 to assign VQSLOD 
    scores to each SNP and filter those that fall below the 99.5% 
    truth sensitivity threshold.
    """
    input:
        recal = '04.filter_variant/VQSR_Filtering/cohort_snp.recal',
        tranches = '04.filter_variant/VQSR_Filtering/cohort_snp.tranches',
        vcf = '03.call_variant/join_calling/cohort_raw.vcf.gz',
        tbi = '03.call_variant/join_calling/cohort_raw.vcf.gz.tbi',
    output:
        vcf = '04.filter_variant/VQSR_Filtering/cohort_snp_recalibrated.vcf.gz',
        tbi = '04.filter_variant/VQSR_Filtering/cohort_snp_recalibrated.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/06.vqsr/ApplyVQSR_SNP.log",
    params:
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        gatk_memory = config['parameter']['gatk_memory']["ApplyVQSR_SNP"],
    threads:
        1
    shell:
        """
        gatk --java-options "{params.gatk_memory}" ApplyVQSR \
            -R {params.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --truth-sensitivity-filter-level 99.5 \
            --tranches-file {input.tranches} \
            --recal-file {input.recal} \
            -mode SNP &> {log}
        """

rule VariantRecalibrator_INDEL:
    """
    Step 3: Build a recalibration model for INDELs.

    Note: The input is the VCF that has already undergone SNP VQSR. 
    Uses Mills and dbSNP as the truth/training sets.
    """
    input:
        vcf = '04.filter_variant/VQSR_Filtering/cohort_snp_recalibrated.vcf.gz',
        tbi = '04.filter_variant/VQSR_Filtering/cohort_snp_recalibrated.vcf.gz.tbi',
    output:
        recal = '04.filter_variant/VQSR_Filtering/cohort_indel.recal',
        tranches = '04.filter_variant/VQSR_Filtering/cohort_indel.tranches',
        plots = '04.filter_variant/VQSR_Filtering/cohort_indel_plots.R',
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/06.vqsr/VariantRecalibrator_INDEL.log",
    params:
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        mills = config['bwa_mem'][config['Genome_Version']]['Mills'],
        dbsnp = config['bwa_mem'][config['Genome_Version']]['dbsnp138'],
        gatk_memory = config['parameter']['gatk_memory']["VariantRecalibrator_INDEL"],
    threads:
        1
    shell:
        """
        gatk --java-options "{params.gatk_memory}" VariantRecalibrator \
            -R {params.ref} \
            -V {input.vcf} \
            -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 \
            -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
            -mode INDEL \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            --rscript-file {output.plots} &> {log}
        """

rule ApplyVQSR_INDEL:
    """
    Step 4: Apply the INDEL recalibration model to produce the final VCF.

    Filters INDELs based on a 99.0% truth sensitivity threshold. The resulting 
    VCF is the final high-quality output of the VQSR pipeline.
    """
    input:
        vcf = '04.filter_variant/VQSR_Filtering/cohort_snp_recalibrated.vcf.gz',
        tbi = '04.filter_variant/VQSR_Filtering/cohort_snp_recalibrated.vcf.gz.tbi',
        recal = '04.filter_variant/VQSR_Filtering/cohort_indel.recal',
        tranches = '04.filter_variant/VQSR_Filtering/cohort_indel.tranches',
    output:
        vcf = '04.filter_variant/cohort_filtered.vcf.gz',
        tbi = '04.filter_variant/cohort_filtered.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/06.vqsr/ApplyVQSR_INDEL.log",
    params:
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        gatk_memory = config['parameter']['gatk_memory']["ApplyVQSR_INDEL"],
    threads:
        1
    shell:
        """
        gatk --java-options {params.gatk_memory} ApplyVQSR \
            -R {params.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file {input.tranches} \
            --recal-file {input.recal} \
            -mode INDEL &> {log}
        """
# ----- rule ----- #