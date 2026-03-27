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
                              reference_fasta_parameter = ['bwa_mem', config['Genome_Version'], 'genome_fa'])

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
    output:
        gvcf = '03.call_variant/join_calling/{sample}.g.vcf.gz',
        tbi = '03.call_variant/join_calling/{sample}.g.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/03.gvcf/{sample}_haplotypecaller.log",
    benchmark:
        "benchmarks/{sample}_haplotypecaller_benchmark.txt",
    message:
        "Running GATK HaplotypeCaller (GVCF mode) for sample: {wildcards.sample}",
    params:
        ploidy = config['ploidy_setting'][config['Genome_Version']]['ploidy'],
        java_opts = "-Xmx16G -XX:ParallelGCThreads=2",
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    threads:
        config['parameter']['threads']['gatk_hc'],
    shell:
        """
        gatk --java-options "{params.java_opts}" HaplotypeCaller \
            -R {params.ref} \
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
    threads:
        1
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
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
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
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
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
        **rule_resource(config,'medium_resource',skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "bcftools.yaml"),
    log:
        "logs/05.vcf/concat_joint_vcf.log",
    threads:
        config['parameter']['threads']['bcftools'],
    shell:
        """
        ( bcftools concat --threads {threads} -O z -o {output.vcf} {input.vcfs} && \
          bcftools index --threads {threads} -t {output.vcf} ) &> {log}
        """
# ---------------------------- VQSR ---------------------------- #
# wheat if we want to use VQSR
# get the list of supported genomes
vqsr_supported_genomes = config.get('join_calling_VQSR_filter', ['hg38.p14'])

# loading VQSR or hard filtering rules
if config['Genome_Version'] in vqsr_supported_genomes and len(samples) > 30:
    logger.info(f"🧬 Genome {config['Genome_Version']} supports VQSR and cohort size ({len(samples)}) is sufficient. Including VQSR rules.")
    workflow.include(os.path.join(workflow.basedir, "rules", "sub_rules", "gatk_VQSR_filtering.smk"))
else:
    logger.info(f"🌱 Genome lacks truth sets or cohort size ({len(samples)}) is too small. Fallback to Hard Filtering.")
    workflow.include(os.path.join(workflow.basedir, "rules", "sub_rules", "gatk_hard_filtering.smk"))

# ---------------------------- VQSR ---------------------------- #

rule gatk_joint_SnpEff_annotation:
    """
    Annotate variants and predict functional impacts using SnpEff.

    This rule assigns biological meaning to the identified variants by 
    comparing them against a gene annotation database.

    Outputs:
    - VCF: Variants with added ANN fields (Effect, Impact, Gene Name, etc.)
    - CSV/HTML Reports: Summary statistics of variant types and predicted impacts
    - Upstream/Downstream (ud) 500bp: Considers flanking regions for annotation
    """
    input:
        vcf = '04.filter_variant/cohort_filtered.vcf.gz',
        tbi = '04.filter_variant/cohort_filtered.vcf.gz.tbi',
    output:
        annotated_csv = '04.filter_variant/cohort_filtered.sort.annotation.csv',
        annotated_html = '04.filter_variant/cohort_filtered.sort.annotation.html',
        annotated_vcf = '04.filter_variant/cohort_filtered.sort.annotation.vcf',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        os.path.join(workflow.basedir, "envs", "snpEff.yaml"),
    message:
        "Running call variant annotation by snpEff",
    log:
        "logs/03.call_variant/merge_variant_filter_annotation.log",
    benchmark:
        "benchmarks/merge_variant_filter_annotation_benchmark.txt",
    threads:    
        config['parameter']["threads"]["snpEff"],
    params:
        config = os.path.join(config['snpEff_db_path'],'snpEff.config'),
        genome_name = config['snpEff_name'][config['Genome_Version']]['name'],
    shell:
        """
        snpEff  -csvStats {output.annotated_csv} \
                -s {output.annotated_html}  \
                -c {params.config} -v \
                -ud 500 {params.genome_name} \
                {input.vcf} > {output.annotated_vcf} 2>{log}
        """

rule gatk_joint_variant_stats:
    """
    Generate quantitative statistics for VCF files using bcftools stats.

    Produces a detailed summary of the VCF content, including transition/transversion 
    ratios (Ti/Tv), allele frequency distributions, and per-sample counts.
    """
    input:
        vcf = '04.filter_variant/cohort_filtered.vcf.gz',
        tbi = '04.filter_variant/cohort_filtered.vcf.gz.tbi',
    output:
        stats = '04.filter_variant/cohort_filtered.vcf.stats.txt',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        os.path.join(workflow.basedir, "envs", "bcftools.yaml"),
    message:
        "Running bcftools stats for {input.vcf} file ",
    log:
        "logs/03.call_variant/cohort_filtered_stats.log",
    benchmark:
        "benchmarks/cohort_filtered_variant_stats_benchmark.txt",
    params:
        reference = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    threads: 
        config['parameter']["threads"]["bcftools"],
    shell:
        """
        bcftools stats --threads {threads} \
                            -f {params.reference} \
                            {input.vcf}  > {output.stats} 2>{log}
        """

rule gatk_joint_variant_stats_multiqc:
    """
    Aggregate variant statistics across all samples using MultiQC.

    This rule collects the output from bcftools stats and generates an 
    aggregated report, enabling cross-sample comparisons of variant quality 
    and distribution.

    Key Visualizations:
    - Number of SNPs/Indels per sample
    - Ti/Tv ratios for quality monitoring
    - Missing data distribution
    """
    input:
        stats = '04.filter_variant/cohort_filtered.vcf.stats.txt',
    output:
        report = '04.filter_variant/multiqc_callvariant_report.html',
    resources:
        **rule_resource(config,'low_resource',skip_queue_on_local=True,logger = logger),
    conda:
        os.path.join(workflow.basedir, "envs", "multiqc.yaml"),
    message:
        "Running MultiQC to aggregate callvariant reports",
    benchmark:
        "benchmarks/multiqc_callvariant_benchmark.txt",
    params:
        fastqc_reports = "04.filter_variant",
        report_dir = "04.filter_variant",
        report = "multiqc_callvariant_report.html",
        title = "multiqc_callvariant_report",
    log:
        "logs/03.call_variant/multiqc_call_varinat.log",
    threads:
        config['parameter']['threads']['multiqc'],
    shell:
        """
        multiqc {params.fastqc_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- END ----- #