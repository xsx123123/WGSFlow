#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
from loguru import logger
from utils.variant_chromosomes import get_chromosomes
# ----- rule ----- #

rule record_ref_metadata:
    """
    Record reference genome metadata for reproducibility.

    This rule creates a text file containing essential reference genome
    information including the genome version, FASTA file path, and sequence
    dictionary information. This metadata is crucial for ensuring reproducibility
    and proper interpretation of results across different analysis runs.
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
        samtools dict  {input.ref} | grep '^@SQ' >> {output.ref_info}
        """

rule bwa_mapping:
    input:
        r1_trimmed = "01.qc/short_read_trim/{sample}.R1.trimed.fq.gz",
        r2_trimmed = "01.qc/short_read_trim/{sample}.R2.trimed.fq.gz",
    output:
        bam = temp('../02.mapping/bwa_mem/{sample}/{sample}.raw.bam'),
    conda:
        workflow.source_path("../envs/bwa.yaml"),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    log:
        "../logs/02.mapping/bwa_mem_{sample}.log",
    message:
        "Running bwa-mem mapping on {wildcards.sample} R1 and R2",
    benchmark:
            "../benchmarks/{sample}_bwa_mem_benchmark.txt",
    params:
        index = config['bwa_mem'][config['Genome_Version']][' index'],
        pl = config['parameter']["bwa_mem"]["library_platform"],
    threads: 
        config['parameter']["threads"]["bwa_mem"],
    shell:
        """
        ( bwa-mem mem -t {threads} \
        -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:{params.pl}\tLB:{wildcards.sample}' \
        {params.index} \
        {input.r1_trimmed} {input.r2_trimmed} | samtools view -@ {threads} \
        -Sbh -o {output.bam} ) 2>{log}
        """

rule sort_index:
    input:
        bam = '02.mapping/bwa_mem/{sample}/{sample}.raw.bam',
    output:
        sort_bam =  temp('02.mapping/bwa_mem/sort_index/{sample}.sort.bam'),
        sort_bam_bai = temp('02.mapping/bwa_mem/sort_index/{sample}.sort.bam.bai'),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    message:
        "Running samtools sort & index for {wildcards.sample}",
    log:
        "logs/02.mapping/bwa_sort_index_{sample}.log",
    benchmark:
            "benchmarks/{sample}_bam_sort_index_benchmark.txt",
    threads: 
        config['parameter']["threads"]["samtools"],
    shell:
        """
        (samtools sort -@ {threads} -o {output.sort_bam} {input.bam} &&
        samtools index -@ {threads} {output.sort_bam}) 2>{log}
        """

rule bam2cram:
    """
    Convert BAM files to CRAM format for storage efficiency.

    This rule compresses the sorted BAM files into CRAM format, which typically
    achieves 40-60% smaller file sizes compared to BAM while maintaining full
    compatibility with most bioinformatics tools. The CRAM format uses reference-
    based compression, making it ideal for large-scale RNA-seq projects where
    storage costs are a concern.
    """
    input:
        sort_bam = '02.mapping/bwa_mem/sort_index/{sample}.sort.bam',
        sort_bam_bai = '02.mapping/bwa_mem/sort_index/{sample}.sort.bam.bai',
    output:
        cram = '02.mapping/cram/{sample}.cram',
        cram_index = '02.mapping/cram/{sample}.cram.crai',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/bwa2.yaml"),
    message:
        "Running bam compress of BAM : {input.bam}",
    log:
        "logs/02.mapping/bam_dup_stats_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_stats_benchmark.txt",
    threads:
        config['parameter']['threads']['bam2cram'],
    params:
        reference = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
    shell:
        """
        samtools view -@ {threads} -C -T {params.reference} -o {output.cram} {input.bam}
        samtools index  -@ {threads} {output.cram}
        """

rule sambamba_MarkDuplicates:
    input:
        sort_bam = '02.mapping/bwa_mem/sort_index/{sample}.sort.bam',
        sort_bam_bai = '02.mapping/bwa_mem/sort_index/{sample}.sort.bam.bai',
    output:
        duplicates_bam = temp('02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam'),
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    log:
        "logs/02.mapping/duplicates_{sample}.log",
    message:
        "Running BAM MarkDuplicates for {input.bam}",
    benchmark:
        "benchmarks/{sample}_bam_MarkDuplicates_benchmark.txt",
    conda:
        workflow.source_path("../envs/sambamba.yaml"),
    threads: 
        config['parameter']['threads']['sambamba'],
    shell:
        """
        sambamba markdup --nthreads {threads} \
                          --show-progress \
                          {input.bam} \
                          {output.duplicates_bam} 2>{log}
        """

rule Duplicates_bam_index:
    input:
        duplicates_bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
    output:
        duplicates_bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    resources:
        **rule_resource(config, 'medium_resource',  skip_queue_on_local=True,logger = logger),
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running build index for MarkDuplicates of BAM {input.bam}",
    log:
        "logs/02.mapping/Duplicates_bam_index_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_index_benchmark.txt",
    threads: 
        config['parameter']['threads']['samtools'],
    shell:
        """
        samtools index -@ {threads} {input.bam} 2>{log}
        """

rule bam_coverage:
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
        "Running Coverage for MarkDuplicates of BAM : {input.bam}",
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
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        qualimap_report_html = '02.mapping/qualimap_report/{sample}/qualimapReport.html',
        qualimap_report_txt = '02.mapping/qualimap_report/{sample}/genome_results.txt',
    conda:
        workflow.source_path("../envs/qualimap.yaml"),
    message:
        "Running qualimap qc for MarkDuplicates of BAM : {input.bam}",
    log:
        "logs/02.mapping/qualimap_report_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_qualimap_benchmark.txt",
    params:
        genome_gff = config['bwa_mem'][config['Genome_Version']]['genome_gff'],
        outformat = config['parameter']['qualimap']["format"],
        mem = config['parameter']['qualimap']["mem"],
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
                 --java-mem-size=16G &>{log}
        """

rule samtools_flagst:
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        samtools_flagstat = '02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv',
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running flagst for MarkDuplicates of BAM : {input.bam}",
    log:
        "logs/02.mapping/bam_dup_lagstat_{sample}.log",
    benchmark:
        "benchmarks/{sample}_Dup_bam_lagstat_benchmark.txt",
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
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        samtools_stats = '../02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv',
    conda:
        workflow.source_path("../envs/samtools.yaml"),
    message:
        "Running stats for MarkDuplicates of BAM : {input.bam}",
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
    Aggregate all mapping QC results into a comprehensive MultiQC report.

    This rule collects QC metrics from all previous mapping-related rules
    (STAR, Qualimap, Samtools, Preseq) and generates a unified HTML report
    using MultiQC. The aggregated report provides a single dashboard for
    assessing the quality of all samples in the RNA-seq experiment, making
    it easy to identify potential issues or outliers.
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
    params:
        fastqc_reports = "02.mapping/",
        report_dir = '02.mapping/mapping_report',
        report = "multiqc_mapping_report.html",
        title = "mapping_report",
    log:
        "logs/02.mapping/multiqc_mapping_report.log",
    benchmark:
        "benchmarks/multiqc_mapping_report_benchmark.txt",
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