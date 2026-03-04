#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - Joint Variant Calling Module (GATK Best Practices)

This module implements the GATK "Joint Calling" workflow, designed to improve 
the sensitivity and specificity of variant discovery across cohorts. Instead 
of calling variants on each sample individually, it produces per-sample GVCFs 
and then performs a joint analysis.

Key Components:
- gatk_haplotypecaller: Generation of per-sample Genomic VCF (GVCF) files.
- create_sample_map: Automated mapping of sample names to GVCF paths for database import.
- gatk_genomicsdbimport: Scalable consolidation of GVCFs into a GenomicsDB workspace.
- gatk_genotypegvcfs: Final joint genotyping to produce a cohort-level VCF.
- concat_joint_vcf: Aggregation of chromosome-level results into a genome-wide VCF.

Joint calling significantly enhances the detection of rare variants and provides 
more accurate genotype calls by leveraging information across all samples in the 
project simultaneously.
"""

import os
from loguru import logger

# ----- rule ----- #
# get chrome list by get_chromosomes function
CHROMOSOMES = get_chromosomes(config=config,
                              reference_fasta=['bwa_mem', config['Genome_Version'], 'genome_fa'])

rule gatk_haplotypecaller:
    """
    Generate per-sample GVCF files using GATK HaplotypeCaller.

    This is the first stage of the Joint Calling workflow. The tool identifies 
    active regions and performs localized assembly to call variants.

    Key Processing Steps:
    - -ERC GVCF: Essential for joint calling; it outputs reference-block 
      information which allows GATK to distinguish between 'no call' and 'ref call'.
    - Localized Assembly: Re-assembles reads in active regions to improve indel accuracy.
    - Native Pair-HMM: Hardware-accelerated genotype likelihood calculation.

    Parameters:
    - ploidy: Sample ploidy (usually 2 for diploid human/mammalian data).
    - java_opts: JVM memory and garbage collection settings.
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    output:
        gvcf = '03.call_variant/join_calling/{sample}.g.vcf.gz',
        tbi = '03.call_variant/join_calling/{sample}.g.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/03.gvcf/{sample}_haplotypecaller.log",
    benchmark:
        "benchmarks/{sample}_haplotypecaller_benchmark.txt",
    message:
        "Running GATK HaplotypeCaller (GVCF mode) for sample: {wildcards.sample}",
    params:
        ploidy = config['ploidy_setting'][config['Genome_Version']]['ploidy'],
        java_opts = "-Xmx16G -XX:ParallelGCThreads=2",
    threads:
        config['parameter']['threads']['gatk_hc'],
    shell:
        """
        gatk --java-options "{params.java_opts}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            --ploidy {params.ploidy} \
            --native-pair-hmm-threads {threads} &> {log}
        """

rule create_sample_map:
    """
    Dynamically generate a sample-name-map file for GenomicsDBImport.

    This map links sample names to their corresponding GVCF absolute paths, 
    serving as a lookup table for the cohort consolidation step.
    """
    input:
        gvcfs = expand('03.call_variant/join_calling/{sample}.g.vcf.gz', sample=samples.keys()),
    output:
        sample_map = '03.call_variant/join_calling/cohort.sample_map',
    run:
        with open(output.sample_map, 'w') as f:
            for gvcf_file in input.gvcfs:
                # Extracts sample name from filename (assumes {sample}.g.vcf.gz)
                sample_name = os.path.basename(gvcf_file).split('.')[0]
                f.write(f"{sample_name}\t{os.path.abspath(gvcf_file)}\n")

rule gatk_genomicsdbimport:
    """
    Consolidate per-sample GVCFs into a GenomicsDB workspace.

    This rule addresses the I/O bottlenecks and scalability issues of 
    multi-sample analysis by importing GVCF data into a specialized database 
    structured for fast access.

    Optimization Features:
    - Interval-based processing: Imported per chromosome to enable parallelization.
    - Batch-size: Limits memory usage by processing a subset of samples at a time.
    - Reader-threads: Parallelizes the reading of GVCF files.
    """
    input:
        sample_map = '03.call_variant/join_calling/cohort.sample_map',
        gvcfs = expand('03.call_variant/join_calling/{sample}.g.vcf.gz', sample=samples.keys()),
        tbis = expand('03.call_variant/join_calling/{sample}.g.vcf.gz.tbi', sample=samples.keys()),
    output:
        db = directory('03.call_variant/join_calling/genomicsdb_{chr}'),
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/04.joint_calling/genomicsdbimport_{chr}.log",
    benchmark:
        "benchmarks/genomicsdbimport_{chr}_benchmark.txt",
    message:
        "Consolidating cohort GVCFs into GenomicsDB for chromosome: {wildcards.chr}",
    params:
        java_opts = "-Xmx16G -Xms16G",
        chr = "{chr}",
        extra_args = "--reader-threads 3 --batch-size 50",
    threads:
        config['parameter']['threads']['gatk_genomicsdb'],
    shell:
        """
        rm -rf {output.db} && \
        gatk --java-options "{params.java_opts}" GenomicsDBImport \
            --sample-name-map {input.sample_map} \
            --genomicsdb-workspace-path {output.db} \
            -L {params.chr} \
            {params.extra_args} &> {log}
        """

rule gatk_genotypegvcfs:
    """
    Perform Joint Genotyping on a consolidated GenomicsDB workspace.

    This is the final step where population-level evidence is used to produce 
    accurate genotype calls and quality scores for all samples simultaneously.

    Key Benefits:
    - Resolves low-confidence sites by looking at other samples in the cohort.
    - Provides a unified multi-sample VCF that identifies where samples match 
      the reference.
    - Enables accurate calculation of population-level metrics (e.g., QualByDepth).
    """
    input:
        db = '03.call_variant/join_calling/genomicsdb_{chr}',
    output:
        vcf = '03.call_variant/join_calling/joint_calling_{chr}.vcf.gz',
        tbi = '03.call_variant/join_calling/joint_calling_{chr}.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'high_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
    log:
        "logs/05.vcf/genotypegvcfs_{chr}.log",
    benchmark:
        "benchmarks/genotypegvcfs_{chr}_benchmark.txt",
    message:
        "Performing joint genotyping on chromosome: {wildcards.chr}",
    params:
        java_opts = "-Xmx32G",
        ref_fa = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    threads:
        config['parameter']['threads']['gatk_genotype'],
    shell:
        """
        gatk --java-options "{params.java_opts}" GenotypeGVCFs \
            -R {params.ref_fa} \
            -V gendb://{input.db} \
            -O {output.vcf} &> {log}
        """

rule concat_joint_vcf:
    """
    Merge chromosome-level joint-called VCFs into a single cohort-level file.

    Aggregates the parallelized outputs into a continuous genomic VCF, 
    maintaining coordinate sorting and generating a global index.
    """
    input:
        vcfs = expand('03.call_variant/join_calling/joint_calling_{chr}.vcf.gz', chr=CHROMOSOMES),
        tbis = expand('03.call_variant/join_calling/joint_calling_{chr}.vcf.gz.tbi', chr=CHROMOSOMES),
    output:
        vcf = '03.call_variant/join_calling/cohort_raw.vcf.gz',
        tbi = '03.call_variant/join_calling/cohort_raw.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    log:
        "logs/05.vcf/concat_joint_vcf.log",
    threads:
        config['parameter']['threads']['bcftools'],
    shell:
        """
        ( bcftools concat --threads {threads} -O z -o {output.vcf} {input.vcfs} && \
          bcftools index --threads {threads} -t {output.vcf} ) &> {log}
        """

# Inclusion of filtering sub-rules based on configuration
if config['join_calling']['filter_type'] == 'hard_filter':
    include "./sub_rules/gatk_hard_filtering.smk"
elif config['join_calling']['filter_type'] == 'VQSR':
    include "./sub_rules/gatk_VQSR_filtering.smk"
else:
    info("No filtering method specified. Skipping filtering step.")

# ----- rule ----- #