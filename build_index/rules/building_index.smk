#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

rule build_bwa_index:
    """

    """
    input:
        md5_check = config["genome_version_info"]["genome_assembly"]["genome_fa"],
    output:
        bwa_index_amb = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_index.amb"),
        bwa_index_ann = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_index.ann"),
        bwa_index_bwt = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_index.bwt"),
        bwa_index_pac = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_index.pac"),
    conda:
        workflow.source_path("../envs/bwa.yaml"),
    log:
        os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_index.log"),
    message:
        "Building {new_genome_name} bwa index .... ",
    benchmark:
        os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_benchmark.txt"),
    params:
        bwa_prefix = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.bwa_index"),
    threads:
        1
    shell:
        """
        bwa index {input.md5_check} -p {params.bwa_prefix}  &> {log}
        """

rule copy_genome_annotation:
    input:
        gtf = config["genome_version_info"]["genome_annotation"]["genome_gtf"],
        gff = config["genome_version_info"]["genome_annotation"]["genome_gff"],
        go  = config["genome_version_info"]["genome_annotation"]["genome_gff"],
    output:
        gtf = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gtf"),
        gff = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gff"),
        go  = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.go"),
    log:
        os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.annotation.log"),
    message:
        "Building {new_genome_name} bwa index .... ",
    threads:
        1
    shell:
        """
        cp {input.gtf} {output.gtf} && \
        cp {input.gff} {output.gff} && \
        cp {input.go} {output.go}  & > {log}
        """

rule create_bed12:
    input:
        gtf = config["genome_version_info"]["genome_annotation"]["genome_gtf"],
    output:
        GenePred = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gtf.GenePred"),
        bed12 = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gtf.bed12"),
    log:
        os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gtf2bed12.log"),
    message:
        "Building {new_genome_name} bwa index .... ",
    conda:
        workflow.source_path("../envs/ucsc_gff.yaml"),
    threads:
        1
    shell:
        """
        gtfToGenePred -genePredExt -ignoreGroupsWithoutExons {input.gtf} {output.GenePred} && \
        genePredToBed {output.GenePred} {output.bed12} & > {log}
        """

rule create_ref_all:
    input:
        gtf = config["genome_version_info"]["genome_annotation"]["genome_gtf"],
    output:
        GenePred = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gtf.GenePred"),
        all_ref = os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.all.txt"),
    log:
        os.path.join(new_genome_path,new_genome_name,f"{new_genome_name}.gtf2bed12.log"),
    message:
        "Building {new_genome_name} bwa index .... ",
    conda:
        workflow.source_path("../envs/ucsc_gff.yaml"),
    threads:
        1
    shell:
        """
        gtfToGenePred -genePredExt -ignoreGroupsWithoutExons {input.gtf} {output.GenePred} && \
        awk 'BEGIN{OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' {output.GenePred} > {output.all_ref} & > {log}
        """
