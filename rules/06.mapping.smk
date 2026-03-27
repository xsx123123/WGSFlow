#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - Short Read Mapping and Alignment Module

This module handles the core alignment workflow, transforming cleaned FASTQ reads 
into coordinate-sorted, deduplicated BAM/CRAM files. It employs industry-standard 
tools for alignment and comprehensive quality assessment.

Key Components:
- bwa_mapping: Primary alignment using the BWA-MEM algorithm
- sort_index: Coordinate sorting and indexing for random access
- sambamba_MarkDuplicates: Identification and marking of PCR/optical duplicates
- bam2cram: Reference-based compression for efficient storage
- mosdepth/qualimap/samtools: Multi-dimensional alignment quality control
- mapping_report: Aggregated QC visualization using MultiQC

The alignment process is the foundation of the WGS/WES pipeline, determining 
the genomic origin of each read. Proper handling of duplicates and rigorous 
QC at this stage are critical for the accuracy of downstream variant calling.
"""

import os
from loguru import logger
from utils.variant_chromosomes import get_chromosomes

# ----- rule ----- #

rule record_ref_metadata:
    """
    Record reference genome metadata for reproducibility.

    This rule captures the exact genomic context used for alignment by 
    documenting the version, file paths, and sequence dictionary.

    Key Features:
    - Extracts @SQ lines (sequence names and lengths) from the FASTA file
    - Stores metadata in a centralized text file
    - Ensures traceability of reference versions across analysis runs
    """
    input:
        ref = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    output:
        ref_info = '02.mapping/cram/reference_version.txt',
    resources:
        **rule_resource(config, 'low_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    log:
        "logs/02.mapping/cram_reference.log",
    benchmark:
        "benchmarks/cram_reference_benchmark.txt",
    threads:
        1
    shell:
        """
        echo "Genome_Version: {config[Genome_Version]}" > {output.ref_info}
        echo "Fasta_Path: {input.ref}" >> {output.ref_info}
        samtools dict {input.ref} | grep '^@SQ' >> {output.ref_info}
        """

rule bwa_mapping:
    """
    Perform primary sequence alignment using BWA-MEM.

    Maps trimmed reads to the reference genome using a Burrows-Wheeler Transform 
    based algorithm. It is optimized for high-quality queries and long reads.

    Key Processing Steps:
    - Addition of Read Group (RG) information (ID, SM, PL, LB)
    - Automatic mapping of paired-end reads
    - Conversion of SAM output to compressed BAM format on-the-fly
    
    Parameters:
    - index: Path to BWA index files
    - pl: Sequencing platform (e.g., ILLUMINA)
    - threads: Parallel processing limit for BWA and Samtools view
    """
    input:
        r1_trimmed = "01.qc/short_read_trim/{sample}.R1.trimed.fq.gz",
        r2_trimmed = "01.qc/short_read_trim/{sample}.R2.trimed.fq.gz",
    output:
        bam = temp('02.mapping/bwa_mem/{sample}/{sample}.raw.bam'),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa.yaml"),
    log:
        "logs/02.mapping/bwa_mem_{sample}.log",
    message:
        "Running bwa-mem mapping on {wildcards.sample} R1 and R2",
    benchmark:
        "benchmarks/{sample}_bwa_mem_benchmark.txt",
    params:
        index = config['bwa_mem'][config['Genome_Version']]['index'],
        pl = config['parameter']["bwa_mem"]["library_platform"],
    threads: 
        config['parameter']["threads"]["bwa_mem"],
    shell:
        """
        ( bwa-mem mem -t {threads} -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:{params.pl}\tLB:{wildcards.sample}' \
                      {params.index} \
                      {input.r1_trimmed} {input.r2_trimmed}  | samtools view -@ {threads} -bh  -o {output.bam} - ) 2>{log}
        """

rule sort_index:
    """
    Sort BAM files by genomic coordinates and generate indices.

    Rearranges reads in the BAM file to follow the order of the reference 
    chromosomes and positions.

    Importance:
    - Required for almost all downstream tools (MarkDuplicates, GATK, IGV)
    - Enables efficient random access via the companion .bai index file
    - Reduces the memory footprint of multi-sample analysis
    """
    input:
        bam = '02.mapping/bwa_mem/{sample}/{sample}.raw.bam',
    output:
        sort_bam =  temp('02.mapping/bwa_mem/sort_index/{sample}.sort.bam'),
        sort_bam_bai = temp('02.mapping/bwa_mem/sort_index/{sample}.sort.bam.bai'),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    log:
        "logs/02.mapping/bwa_sort_index_{sample}.log",
    message:
        "Running samtools sort & index for {wildcards.sample}",
    benchmark:
        "benchmarks/{sample}_bam_sort_index_benchmark.txt",
    threads: 
        config['parameter']["threads"]["samtools"],
    shell:
        """
        (samtools sort -@ {threads} -o {output.sort_bam} {input.bam} &&
        samtools index -@ {threads} {output.sort_bam}) 2>{log}
        """

rule sambamba_MarkDuplicates:
    """
    Identify and mark PCR or optical duplicate reads.

    Duplicate reads can lead to false confidence in variant calls. This rule 
    uses Sambamba, a high-performance alternative to Picard MarkDuplicates.

    Impact on Data:
    - Flags reads with identical 5' positions (considering orientation)
    - Allows downstream variant callers to ignore biased evidence
    - Provides statistics on the duplication rate of the library
    """
    input:
        sort_bam = '02.mapping/bwa_mem/sort_index/{sample}.sort.bam',
        sort_bam_bai = '02.mapping/bwa_mem/sort_index/{sample}.sort.bam.bai',
    output:
        duplicates_bam = temp('02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam'),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/sambamba.yaml"),
    log:
        "logs/02.mapping/duplicates_{sample}.log",
    message:
        "Running BAM MarkDuplicates for {input.sort_bam}",
    benchmark:
        "benchmarks/{sample}_bam_MarkDuplicates_benchmark.txt",
    threads: 
        config['parameter']['threads']['sambamba'],
    shell:
        """
        sambamba markdup --nthreads {threads} \
                          --show-progress \
                          {input.sort_bam} \
                          {output.duplicates_bam} 2>{log}
        """

rule Duplicates_bam_index:
    """
    Generate indices for deduplicated BAM files.

    Ensures that the BAM file containing duplicate flags is ready for 
    high-speed regional queries.
    """
    input:
        duplicates_bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
    output:
        duplicates_bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running build index for deduplicated BAM: {input.duplicates_bam}",
    log:
        "logs/02.mapping/Duplicates_bam_index_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_index_benchmark.txt",
    threads: 
        config['parameter']['threads']['samtools'],
    shell:
        """
        samtools index -@ {threads} {input.duplicates_bam} 2>{log}
        """

rule bam2cram:
    """
    Convert BAM files to CRAM format for archival and storage efficiency.

    CRAM is a highly compressed format that stores only the differences 
    between the reads and the reference genome.

    Key Benefits:
    - Achieves 40-70% storage savings compared to BAM
    - Full compatibility with modern bioinformatics toolchains
    - Lossless preservation of all read data and quality scores
    """
    input:
        duplicates_bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        duplicates_bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        cram = '02.mapping/cram/{sample}.Dup.cram',
        cram_index = '02.mapping/cram/{sample}.Dup.cram.crai',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa2.yaml"),
    message:
        "Running conversion to CRAM for: {input.duplicates_bam}",
    log:
        "logs/02.mapping/bam2cram_{sample}.log",
    benchmark:
        "benchmarks/{sample}_bam2cram_benchmark.txt",
    params:
        reference = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    threads:
        config['parameter']['threads']['bam2cram'],
    shell:
        """
        samtools view -@ {threads} -C -T {params.reference} -o {output.cram} {input.duplicates_bam}
        samtools index -@ {threads} {output.cram}
        """

rule bam_coverage:
    """
    Calculate genomic coverage and depth statistics using Mosdepth.

    Mosdepth is a lightning-fast tool for determining how many reads cover 
    each base in the genome.

    Metrics Generated:
    - Global distribution of coverage (how much of the genome is at X depth)
    - Summary of average depth per chromosome/contig
    - Quantized coverage tracks for visualization
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        dist = "02.mapping/mosdepth_coverage/{sample}.sort.Dup.global.dist.txt",
        summary = "02.mapping/mosdepth_coverage/{sample}.sort.Dup.mosdepth.summary.txt",
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/mosdepth.yaml"),
    message:
        "Running Coverage analysis for: {input.bam}",
    log:
        "logs/02.mapping/bam_coverage_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_coverage_benchmark.txt",
    params:
        prefix = "02.mapping/mosdepth_coverage/{sample}"
    threads: 
        config['parameter']["threads"]["mosdepth"],
    shell:
        """
        mosdepth -n \
                 --fast-mode \
                 -t {threads} \
                 --by 500 \
                 {params.prefix} \
                 {input.bam} &>{log}
        """

rule qualimap_qc:
    """
    Comprehensive alignment quality assessment using Qualimap.

    Provides a macroscopic view of the alignment quality, helping to detect 
    systemic issues like GC bias or library complexity problems.

    Report Highlights:
    - Mapping quality distribution
    - Insert size histogram (crucial for paired-end library validation)
    - Coverage profile across the genome
    - Per-chromosome statistics
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        qualimap_report_html = '02.mapping/qualimap_report/{sample}/qualimapReport.html',
        qualimap_report_txt = '02.mapping/qualimap_report/{sample}/genome_results.txt',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/qualimap.yaml"),
    message:
        "Running qualimap QC for: {input.bam}",
    log:
        "logs/02.mapping/qualimap_report_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_qualimap_benchmark.txt",
    params:
        outformat = config['parameter']['qualimap']["format"],
        mem = config['parameter']['qualimap']["mem"],
        genome_gff = config['bwa_mem'][config['Genome_Version']]['genome_gff'],
        prefix_dir = '02.mapping/qualimap_report/{sample}/',
    threads: 
        config['parameter']["threads"]["qualimap"],
    shell:
        """
        qualimap bamqc \
                 -nt {threads} \
                 -bam {input.bam} \
                 -gff {params.genome_gff} \
                 -outdir {params.prefix_dir} \
                 -outformat {params.outformat} \
                 --java-mem-size={params.mem} &>{log}
        """

rule samtools_flagst:
    """
    Generate alignment status summaries using Samtools flagstat.

    Provides a high-level summary of the BAM file based on bitwise flags.

    Key Statistics:
    - Total reads and percentage mapped
    - Properly paired reads (indicates library health)
    - Secondary and supplementary alignment counts
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        samtools_flagstat = '02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running flagstat for: {input.bam}",
    log:
        "logs/02.mapping/bam_dup_flagstat_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_flagstat_benchmark.txt",
    threads: 
        config['parameter']["threads"]["samtools_flagstat"],
    shell:
        """
        samtools flagstat \
                 -@ {threads} \
                 -O tsv \
                 {input.bam} > {output.samtools_flagstat} 2>{log}
        """

rule samtools_stats:
    """
    Detailed alignment metric collection using Samtools stats.

    Produces a comprehensive set of statistics that go far beyond flagstat, 
    including cycle-wise quality and GC content.

    Detailed Analytics:
    - Base quality distribution across cycles
    - Fragment length distribution for paired reads
    - Mismatches per cycle (identifies potential sequencing drift)
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        samtools_stats = '02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running stats for: {input.bam}",
    log:
        "logs/02.mapping/bam_dup_stats_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_stats_benchmark.txt",
    threads: 
        config['parameter']["threads"]["samtools_stats"],
    params:
        reference = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        samtools stats \
                 -@ {threads} \
                 --reference {params.reference} \
                 {input.bam} > {output.samtools_stats} 2>{log}
        """

rule mapping_report:
    """
    Aggregate all alignment QC reports across all samples using MultiQC.

    This rule collects results from Qualimap, Flagstat, Samtools stats, and 
    Mosdepth to create a single, interactive dashboard for project-wide QC.

    The MultiQC Report enables:
    - Comparison of mapping rates across the entire cohort
    - Detection of outliers with unusual insert sizes or coverage
    - Batch effect identification via visualization of depth and GC bias
    - Validation of the alignment stage before proceeding to variant calling
    """
    input:
        qualimap_report_html =  expand('02.mapping/qualimap_report/{sample}/qualimapReport.html',sample=samples.keys()),
        qualimap_report_txt =  expand('02.mapping/qualimap_report/{sample}/genome_results.txt',sample=samples.keys()),
        samtools_flagstat = expand('02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv',sample=samples.keys()),
        samtools_stats = expand('02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv',sample=samples.keys()),
        dist = expand("02.mapping/mosdepth_coverage/{sample}.sort.Dup.global.dist.txt",sample=samples.keys()),
        summary = expand("02.mapping/mosdepth_coverage/{sample}.sort.Dup.mosdepth.summary.txt",sample=samples.keys()),
    output:
        report = "02.mapping/mapping_report/multiqc_mapping_report.html",
    resources:
        **rule_resource(config, 'low_resource',skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/multiqc.yaml"),
    message:
        "Running MultiQC to aggregate mapping reports",
    log:
        "logs/02.mapping/multiqc_mapping_report.log",
    benchmark:
        "benchmarks/multiqc_mapping_report_benchmark.txt",
    params:
        fastqc_reports = "02.mapping/",
        report_dir = '02.mapping/mapping_report',
        report = "multiqc_mapping_report.html",
        title = "mapping_report",
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