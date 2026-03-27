#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from loguru import logger

rule gatk_haplotypecaller_gvcf_single:
    """
    Generate per-sample GVCF files using GATK HaplotypeCaller.

    Even for single-sample analysis, producing a GVCF first is recommended 
    as it allows for future joint-calling and better handling of reference 
    confidence evidence.
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        gvcf = temp('03.call_variant/single_calling/{sample}.g.vcf.gz'),
        tbi = temp('03.call_variant/single_calling/{sample}.g.vcf.gz.tbi'),
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/03.call_variant/single_calling/{sample}_haplotypecaller_gvcf.log",
    benchmark:
        "benchmarks/{sample}_single_haplotypecaller_gvcf_benchmark.txt",
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

rule gatk_genotypegvcfs_single:
    """
    Perform genotyping on a single-sample GVCF to produce a standard VCF.
    """
    input:
        gvcf = '03.call_variant/single_calling/{sample}.g.vcf.gz',
        tbi = '03.call_variant/single_calling/{sample}.g.vcf.gz.tbi',
    output:
        vcf = '03.call_variant/single_calling/{sample}.vcf.gz',
        tbi = '03.call_variant/single_calling/{sample}.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/03.call_variant/single_calling/{sample}_genotypegvcfs.log",
    benchmark:
        "benchmarks/{sample}_single_genotypegvcfs_benchmark.txt",
    params:
        java_opts = "-Xmx16G",
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    threads:
        1
    shell:
        """
        gatk --java-options "{params.java_opts}" GenotypeGVCFs \
            -R {params.ref} \
            -V {input.gvcf} \
            -O {output.vcf} &> {log}
        """

rule gatk_single_hard_filtering:
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
        vcf = '03.call_variant/single_calling/{sample}.vcf.gz',
        tbi = '03.call_variant/single_calling/{sample}.vcf.gz.tbi',
    output:
        vcf = '04.filter_variant/single_sample/{sample}_filtered.vcf.gz',
        tbi = '04.filter_variant/single_sample/{sample}_filtered.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml"),
    log:
        "logs/03.call_variant/single_sample/{sample}_gatk_hard_filtering.log",
    benchmark:
        "benchmarks/{sample}_single_gatk_hard_filtering_benchmark.txt",
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

rule gatk_SnpEff_annotation:
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
        vcf = '04.filter_variant/single_sample/{sample}_filtered.vcf.gz',
        tbi = '04.filter_variant/single_sample/{sample}_filtered.vcf.gz.tbi',
    output:
        annotated_csv = '04.filter_variant/single_sample/{sample}_filtered.sort.annotation.csv',
        annotated_html = '04.filter_variant/single_sample/{sample}_filtered.sort.annotation.html',
        annotated_vcf = '04.filter_variant/single_sample/{sample}_filtered.sort.annotation.vcf',
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

rule gatk_variant_stats:
    """
    Generate quantitative statistics for VCF files using bcftools stats.

    Produces a detailed summary of the VCF content, including transition/transversion 
    ratios (Ti/Tv), allele frequency distributions, and per-sample counts.
    """
    input:
        vcf = '04.filter_variant/single_sample/{sample}_filtered.vcf.gz',
        tbi = '04.filter_variant/single_sample/{sample}_filtered.vcf.gz.tbi',
    output:
        stats = '04.filter_variant/single_sample/{sample}_filtered.vcf.stats.txt',
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

rule gatk_variant_stats_multiqc:
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
        stats = expand('04.filter_variant/single_sample/{sample}_filtered.vcf.stats.txt',sample=samples.keys()),
    output:
        report = '04.filter_variant/single_sample/multiqc_callvariant_report.html',
    resources:
        **rule_resource(config,'low_resource',skip_queue_on_local=True,logger = logger),
    conda:
        os.path.join(workflow.basedir, "envs", "multiqc.yaml"),
    message:
        "Running MultiQC to aggregate callvariant reports",
    benchmark:
        "benchmarks/multiqc_callvariant_benchmark.txt",
    params:
        fastqc_reports = "04.filter_variant/single_sample",
        report_dir = "04.filter_variant/single_sample",
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
# ----- rule ----- #

