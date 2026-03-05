#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - Germline Variant Calling and Annotation Module

This module implements a robust variant discovery workflow using the bcftools 
suite and functional annotation using SnpEff. It is designed for high-performance 
parallel processing across chromosomes.

Key Components:
- call_variant_by_chromosome: Parallelized per-chromosome variant discovery
- concat_sample_chromosomes: Reconstruction of sample-level VCF files
- merg_vcf: Joint analysis preparation via multi-sample merging
- merg_vcf_filter: Stringent quality filtering (MAF, missingness)
- SnpEff_annotation: Functional impact prediction and genomic annotation
- variant_stats_multiqc: Aggregated variant quality control reporting

The module transforms alignment data (BAM) into actionable biological insights 
by identifying SNPs and Indels, filtering out low-quality sites, and predicting 
the functional consequences of each variant.
"""

import os
from loguru import logger

# ----- rule ----- #
# get chrome list by get_chromosomes function
CHROMOSOMES = get_chromosomes(config = config,
                              reference_fasta = ['bwa_mem',config['Genome_Version'],'genome_fa'])

rule call_variant_by_chromosome:
    """
    Perform per-chromosome variant calling using bcftools mpileup and call.

    This rule parallelizes the computationally intensive calling process by 
    processing individual chromosomes independently.

    Key Processing Steps:
    - mpileup: Generates genotype likelihoods from alignment data
    - call: Performs actual genotype calling using a multiallelic caller model
    - Region-limited processing (-r) for efficient resource utilization
    
    Parameters:
    - ploidy: Genomic ploidy setting (e.g., diploid)
    - reference: Path to the indexed reference genome FASTA
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        vcf = temp('03.call_variant/by_chr/{sample}.{chr}.raw.vcf.gz'),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    log:
        "logs/03.call_variant/{sample}_{chr}_variant.log",
    benchmark:
        "benchmarks/{sample}_{chr}_variant_benchmark.txt",
    message:
        "Running bcftools mpileup & call on {wildcards.sample} chromosome {wildcards.chr}",
    params:
        reference = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        ploidy = config['ploidy_setting'][config['Genome_Version']]['ploidy'],
        verbosity = config['parameter']["bcftools"]['verbosity'],
        chr = "{chr}",
    threads:
        config['parameter']["threads"]["bcftools"],
    shell:
        """
        bcftools mpileup \
                 --threads {threads} \
                 --verbosity {params.verbosity} \
                 -r {params.chr} \
                 -f {params.reference} \
                 -Ou \
                 {input.bam} | bcftools call \
                 --threads {threads} \
                 --ploidy {params.ploidy} \
                 --verbosity {params.verbosity} \
                 -mv -Oz -o {output.vcf} &>{log}
        """

rule concat_sample_chromosomes:
    """
    Concatenate per-chromosome VCF files into a single sample-level VCF.

    Reassembles the fragmented chromosome VCFs produced in the parallel 
    calling step into a complete genomic record for each sample.
    """
    input:
        vcf_list = lambda wildcards: expand('03.call_variant/by_chr/{{sample}}.{chr}.raw.vcf.gz',
                                             chr=CHROMOSOMES),
    output:
        vcf = '03.call_variant/{sample}.raw.vcf.gz',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    log:
        "logs/03.call_variant/{sample}_concat_chr.log",
    benchmark:
        "benchmarks/{sample}_concat_chr_benchmark.txt",
    message:
        "Concatenating chromosome VCFs for sample: {wildcards.sample}",
    threads: 
        config['parameter']["threads"]["bcftools"],
    shell:
        """
        bcftools concat \
          --threads {threads} \
          -Oz -o {output.vcf} \
          {input.vcf_list} &>{log}
        """

rule sort_index_bcftools:
    """
    Coordinate-sort and index sample-level VCF files.

    Ensures the VCF file is properly ordered by chromosome and position, 
    and generates both TBI and CSI indices for broad tool compatibility.
    """
    input:
        vcf = '03.call_variant/{sample}.raw.vcf.gz',
    output:
        sort_vcf = '03.call_variant/{sample}.sort.vcf.gz',
        sort_vcf_index = '03.call_variant/{sample}.sort.vcf.csi',
        sort_vcf_tbi_index = '03.call_variant/{sample}.sort.vcf.tbi',
    resources:
        **rule_resource(config,'medium_resource',skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    log:
        "logs/03.call_variant/{sample}_variant_sort_index.log",
    message:
        "Running bcftools sort & index for sample : {wildcards.sample}",
    benchmark:
        "benchmarks/{sample}_variant_sort_index_benchmark.txt",
    threads: 
        config['parameter']["threads"]["bcftools"],
    shell:
        """
        ( bcftools sort \
          {input.vcf} -O z -o {output.sort_vcf}  && \
          bcftools index --threads {threads} \
          -t {output.sort_vcf} -o {output.sort_vcf_tbi_index} && \
          bcftools index --threads {threads} \
          -c {output.sort_vcf} -o {output.sort_vcf_index} ) &>{log}
        """

rule merg_vcf:
    """
    Merge multiple sample VCFs into a single multi-sample VCF.

    Combines individual sample results into a unified matrix, allowing 
    for population-level analysis and joint filtering.

    Process:
    - Merges files while maintaining allele consistency
    - Re-sorts the merged output
    - Generates required indices for the unified file
    """
    input: 
        vcf_list = expand('03.call_variant/{sample}.sort.vcf.gz', sample=samples.keys()),
    output:
        merge_vcf = '03.call_variant/merge.vcf.gz',
        merge_sort_vcf = '03.call_variant/merge.sort.vcf.gz',
        merge_index_vcf = '03.call_variant/merge.sort.vcf.gz.csi',
        merge_tbi_index_vcf = '03.call_variant/merge.sort.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    message:
        "Running bcftools merge & sort & index on vcf : {input.vcf_list}",
    log:
        "logs/03.call_variant/merge_variant_sort_index.log",
    benchmark:
        "benchmarks/merge_variant_sort_index_benchmark.txt",
    threads: 
        config['parameter']["threads"]["bcftools"],
    shell:
        """
        ( bcftools merge --threads {threads} \
          {input.vcf_list} -O z -o {output.merge_vcf}  && \
          bcftools sort \
          {output.merge_vcf} -O z -o {output.merge_sort_vcf} && \
          bcftools index --threads {threads} \
          -t {output.merge_sort_vcf} -o {output.merge_tbi_index_vcf} && \
          bcftools index --threads {threads} \
          -c {output.merge_sort_vcf} -o {output.merge_index_vcf} ) &>{log}
        """

def get_dynamic_bcftools_filter(wildcards):
    """
    根据样本数量动态生成 bcftools 的过滤参数。
    - 小样本 (<=15): 仅进行基础质量和相对于样本量的总深度过滤。
    - 大群体 (>15): 增加 MAF 和缺失率 (F_MISSING) 过滤。
    """
    num_samples = len(samples)
    
    # 基础硬过滤条件：
    # QUAL > 30: 质量门槛
    # DP > {num_samples * 5}: 确保全队列平均深度至少在 5x 以上
    min_total_dp = num_samples * 5
    base_filter = f"QUAL > 30 && DP > {min_total_dp}"
    
    if num_samples > 15:
        # 大样本队列：应用群体遗传学过滤标准
        logger.info(f"Cohort size is {num_samples} (>15). Applying population-level filters (MAF 0.05, Missingness 10%).")
        return f"-i '{base_filter} && F_MISSING <= 0.1 && MAF >= 0.05'"
    else:
        # 小样本队列：放宽标准以保留稀有/私有变异
        logger.info(f"Cohort size is {num_samples} (<=15). Applying relaxed filters to preserve private variants.")
        return f"-i '{base_filter}'"

rule merg_vcf_filter:
    """
    Apply quality filters to the multi-sample VCF.

    Dynamically refines the variant set:
    - If samples > 15: Applies base quality (QUAL, DP) + Population metrics (F_MISSING, MAF)
    - If samples <= 15: Applies base quality only (preserves private variants)
    """
    input: 
        merge_vcf = '03.call_variant/merge.sort.vcf.gz',
    output:
        filter_vcf = '04.filter_variant/merge_filter.sort.vcf.gz',
        filter_index_vcf = '04.filter_variant/merge_filter.sort.vcf.gz.csi',
        filter_tbi_index_vcf = '04.filter_variant/merge_filter.sort.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    params:
        filter_expr = get_dynamic_bcftools_filter
    message:
        "Running dynamic bcftools filter for {num_samples} samples: {params.filter_expr}"
    log:
        "logs/03.call_variant/merge_variant_filter.log",
    benchmark:
        "benchmarks/merge_variant_filter_benchmark.txt",
    threads: 
        config['parameter']["threads"]["bcftools"],
    shell:
        """
        ( bcftools view --threads {threads} \
          -v snps,indels \
          {params.filter_expr} \
          {input.merge_vcf} -O z -o {output.filter_vcf} && \
          bcftools index --threads {threads} \
          -t {output.filter_vcf} -o {output.filter_tbi_index_vcf} && \
          bcftools index --threads {threads} \
          -c {output.filter_vcf} -o {output.filter_index_vcf}
          ) &>{log}
        """

rule bcf_SnpEff_annotation:
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
        filter_vcf = '04.filter_variant/merge_filter.sort.vcf.gz',
        filter_index_vcf = '04.filter_variant/merge_filter.sort.vcf.gz.csi',
        filter_tbi_index_vcf = '04.filter_variant/merge_filter.sort.vcf.gz.tbi',
    output:
        annotated_csv = '04.filter_variant/merge_filter.sort.annotation.csv',
        annotated_html = '04.filter_variant/merge_filter.sort.annotation.html',
        annotated_vcf = '04.filter_variant/merge_filter.sort.annotation.vcf',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/snpEff.yaml"),
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
                {input.filter_vcf} > {output.annotated_vcf} 2>{log}
        """

rule bcf_variant_stats:
    """
    Generate quantitative statistics for VCF files using bcftools stats.

    Produces a detailed summary of the VCF content, including transition/transversion 
    ratios (Ti/Tv), allele frequency distributions, and per-sample counts.
    """
    input:
        vcf = '03.call_variant/{sample}.sort.vcf.gz',
    output:
        stats = '03.call_variant/variant_stats/{sample}.vcf.stats.txt',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bcftools.yaml"),
    message:
        "Running bcftools stats for {input.vcf} file ",
    log:
        "logs/03.call_variant/{sample}_variant_stats.log",
    benchmark:
        "benchmarks/{sample}_variant_stats_benchmark.txt",
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

rule bcf_variant_stats_multiqc:
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
        vcf_stats = expand("03.call_variant/variant_stats/{sample}.vcf.stats.txt", sample=samples.keys()),
    output:
        report = '03.call_variant/multiqc_variant_stats/multiqc_callvariant_report.html',
    resources:
        **rule_resource(config,'low_resource',skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate callvariant reports",
    benchmark:
        "benchmarks/multiqc_callvariant_benchmark.txt",
    params:
        fastqc_reports = "03.call_variant/variant_stats/",
        report_dir = "03.call_variant/multiqc_variant_stats/",
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