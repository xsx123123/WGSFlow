#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule all_qc_multiqc_report:
    input:
        md5_check = "../01.qc/md5_check.tsv",
        fastqc_files_r1 = expand("../01.qc/short_read_qc_r2/{sample}_R2_fastqc.zip", sample=samples.keys()),
        fastqc_files_r2 = expand("../01.qc/short_read_qc_r2/{sample}_R2_fastqc.zip", sample=samples.keys()),
        fastp_report = expand("../01.qc/short_read_trim/{sample}.trimed.html", sample=samples.keys()),
    output:
        report = "../F.report/Data_QC_report/Data_QC_report.html"
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to raw data qc result",
    benchmark:
        "../benchmarks/all_qc_multiqc_report_benchmark.txt",
    params:
        report_dir = "../F.report/Data_QC_report/",
        fastqc_reports = "../01.qc/",
        report = "Data_QC_report.html",
        title = "Data_QC_report",
    log:
        "../logs/F.report/multiqc_trim.log",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

rule mapping_qc_multiqc_report:
    input:
        qualimap_report_html = expand('../02.mapping/qualimap_report/{sample}/qualimapReport.html', sample=samples.keys()),
        qualimap_report_txt = expand('../02.mapping/qualimap_report/{sample}/genome_results.txt', sample=samples.keys()),
        samtools_flagstat = expand('../02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv', sample=samples.keys()),
        samtools_stats = expand('../02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv', sample=samples.keys()),
    output:
        report = "../F.report/Mapping_QC_report/Mapping_QC_report.html",
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to mapping reports",
    benchmark:
        "../benchmarks/mapping_qc_multiqc_report_benchmark.txt",
    params:
        report_dir = "../F.report/Mapping_QC_report/",
        fastqc_reports = "../02.mapping/",
        report = "Mapping_QC_report.html",
        title = "Mapping_QC_report",
    log:
        "../logs/F.report/Mapping_QC_report.log",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

rule variant_stats_multiqc_report:
    input:
        fastqc_files_r1 = expand("../03.call_variant/variant_stats/{sample}.vcf.stats.txt", sample=samples.keys()),
    output:
        report = "../F.report/Mapping_QC_report/Variant_Stats_report.html",
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to bcftools variant stats",
    params:
        report_dir = "../F.report/Mapping_QC_report/",
        bcftools_reports = "../03.call_variant/variant_stats/",
        report = "Variant_Stats_report.html",
        title = "Variant_Stats_report",
    log:
        "../logs/F.report/multiqc-variant-stats.log",
    benchmark:
        "../benchmarks/F.report_multiqc-variant-stats_benchmark.txt",
    shell:
        """
        multiqc {params.bcftools_reports} \
                --force \
                --outdir {params.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

# ----- rule ----- #