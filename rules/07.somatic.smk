#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - Somatic Variant Calling Module (Mutect2)

This module implements the GATK Best Practices for somatic short variant 
discovery (SNVs and Indels). It handles Tumor-Normal paired analysis, 
including orientation bias modeling, contamination estimation, and 
rigorous filtering.

Key Components:
- mutect2_calling: Somatic variant discovery with HMM acceleration
- learn_read_orientation: FFPE and oxidation artifact modeling
- get_pileup_summaries: Characterization of population-frequency sites
- calculate_contamination: Estimation of sample cross-contamination
- filter_mutect_calls: Final integration of all models for high-confidence calls

The somatic workflow is designed to distinguish true low-frequency somatic 
mutations from sequencing artifacts and germline background.
"""

import os
from loguru import logger

# ----- helper functions ----- #

def get_contrast_inputs(wildcards):
    """
    Retrieve Tumor and Normal BAM paths for a specific contrast.
    """
    if wildcards.contrast not in CONTRAST_MAP:
        raise ValueError(f"Unknown contrast identifier: {wildcards.contrast}")
    
    mapping = CONTRAST_MAP[wildcards.contrast]
    return {
        "tumor_bam": mapping["tumor_bam"],
        "normal_bam": mapping["normal_bam"]
    }

# ----- rules ----- #

rule mutect2_calling:
    """
    Perform somatic variant discovery using GATK Mutect2.

    Operates in Tumor-Normal paired mode to identify candidate somatic 
    mutations while accounting for site-specific noise and germline resource 
    evidence.

    Key Features:
    - -ERC mode: Not used for somatic, direct VCF output
    - Native Pair-HMM: Hardware-accelerated likelihood calculation
    - F1R2: Captures read orientation information for bias modeling
    """
    input:
        unpack(get_contrast_inputs),
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    output:
        vcf = "03.call_variant/somatic/{contrast}.unfiltered.vcf.gz",
        f1r2 = "03.call_variant/somatic/{contrast}.f1r2.tar.gz",
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/03.somatic/{contrast}_mutect2.log",
    message:
        "Running GATK Mutect2 for contrast: {wildcards.contrast}",
    benchmark:
        "benchmarks/somatic/{contrast}_mutect2_benchmark.txt",
    params:
        normal_name = lambda wildcards: CONTRAST_MAP[wildcards.contrast]["normal_name"],
        java_opts = "-Xmx16G -XX:ParallelGCThreads=2",
        pon = config['ploidy_setting'][config['Genome_Version']]['pon'],
        germline = config['ploidy_setting'][config['Genome_Version']]['germline'],
    threads: 
        config['parameter']['threads'].get('mutect2', 8)
    shell:
        """
        gatk --java-options "{params.java_opts}" Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -normal {params.normal_name} \
            --panel-of-normals {params.pon} \
            --germline-resource {params.germline} \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} &> {log}
        """

rule learn_read_orientation:
    """
    Model read orientation bias to filter sequencing artifacts.

    Learns the prior probability of orientation bias (e.g., CtoA/GtoT 
    oxidation) across the genome to improve somatic filtering.
    """
    input:
        f1r2 = "03.call_variant/somatic/{contrast}.f1r2.tar.gz",
    output:
        rom = "04.filter_variant/somatic/{contrast}.read_orientation_model.tar.gz",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/04.somatic/{contrast}_learn_read_orientation.log",
    message:
        "Learning read orientation model for: {wildcards.contrast}",
    benchmark:
        "benchmarks/somatic/{contrast}_learn_read_orientation_benchmark.txt",
    threads: 
        1
    shell:
        """
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.rom} &> {log}
        """

rule get_pileup_summaries:
    """
    Summarize pileups of variants for contamination estimation.

    Counts reads supporting ref and alt alleles at common variant sites 
    to facilitate the calculation of cross-sample contamination.
    """
    input:
        bam = lambda wildcards: CONTRAST_MAP[wildcards.contrast]["tumor_bam"] if wildcards.type == "tumor" else CONTRAST_MAP[wildcards.contrast]["normal_bam"],
    output:
        pileup = "03.call_variant/somatic/{contrast}.{type}_pileup.table",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/03.somatic/{contrast}_{type}_get_pileup_summaries.log",
    message:
        "Summarizing pileups for {wildcards.contrast} ({wildcards.type})",
    benchmark:
        "benchmarks/somatic/{contrast}_{type}_get_pileup_summaries_benchmark.txt",
    wildcard_constraints:
        type = "tumor|normal"
    params:
        common_variants = config['ploidy_setting'][config['Genome_Version']]['common_variants_vcf'],
    threads:
        1
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {params.common_variants} \
            -L {params.common_variants} \
            -O {output.pileup} &> {log}
        """

rule calculate_contamination:
    """
    Estimate sample contamination and tumor segment copy number.

    Uses pileup summaries from tumor and matched normal to quantify 
    foreign DNA and identify segments with variant allele frequency shifts.
    """
    input:
        tumor_pileup = "03.call_variant/somatic/{contrast}.tumor_pileup.table",
        normal_pileup = "03.call_variant/somatic/{contrast}.normal_pileup.table",
    output:
        contamination = "03.call_variant/somatic/{contrast}.contamination.table",
        segmentation = "03.call_variant/somatic/{contrast}.segments.table",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/03.somatic/{contrast}_calculate_contamination.log",
    message:
        "Calculating contamination for contrast: {wildcards.contrast}",
    benchmark:
        "benchmarks/somatic/{contrast}_calculate_contamination_benchmark.txt",
    threads:
        1
    shell:
        """
        gatk CalculateContamination \
            -I {input.tumor_pileup} \
            -matched {input.normal_pileup} \
            -O {output.contamination} \
            --tumor-segmentation {output.segmentation} &> {log}
        """

rule filter_mutect_calls:
    """
    Apply final high-confidence somatic filters.

    Integrates standard filters with orientation bias models and 
    contamination estimates to produce the final somatic VCF.
    """
    input:
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        unfiltered_vcf = "03.call_variant/somatic/{contrast}.unfiltered.vcf.gz",
        contamination = "03.call_variant/somatic/{contrast}.contamination.table",
        segmentation = "03.call_variant/somatic/{contrast}.segments.table",
        rom = "04.filter_variant/somatic/{contrast}.read_orientation_model.tar.gz",
    output:
        filtered_vcf = "04.filter_variant/somatic/{contrast}.somatic_filtered.vcf.gz",
        filtered_tbi = "04.filter_variant/somatic/{contrast}.somatic_filtered.vcf.gz.tbi",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/04.somatic/{contrast}_filter_mutect_calls.log",
    message:
        "Applying final filters to somatic calls for: {wildcards.contrast}",
    benchmark:
        "benchmarks/somatic/{contrast}_filter_mutect_calls_benchmark.txt",
    threads:
        1
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.unfiltered_vcf} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segmentation} \
            --ob-priors {input.rom} \
            -O {output.filtered_vcf} &> {log}
        """
# ----- END ----- #
