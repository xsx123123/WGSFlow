#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import yaml

# 定义输出路径变量
# 遵循结构: {reference_path}/{species}/{genome_name}/
genome_dir = os.path.join(config["reference_path"], config["genome_version_info"]["species"], config["genome_version_info"]["genome_name"])
# 遵循结构: {reference_path}/{species}/bwa-mem_v11/
bwa_dir = os.path.join(config["reference_path"], config["genome_version_info"]["species"], "bwa-mem_v11")
genome_name = config["genome_version_info"]["genome_name"]

rule copy_fasta:
    input:
        fasta = config["genome_version_info"]["genome_assembly"]["genome_fa"]
    output:
        fasta = os.path.join(genome_dir, f"{genome_name}.fasta")
    log:
        os.path.join(genome_dir, f"{genome_name}.copy_fasta.log")
    shell:
        "cp {input.fasta} {output.fasta} &> {log}"

rule samtools_faidx:
    input:
        os.path.join(genome_dir, f"{genome_name}.fasta")
    output:
        os.path.join(genome_dir, f"{genome_name}.fasta.fai")
    conda:
        "../../envs/samtools.yaml"
    log:
        os.path.join(genome_dir, f"{genome_name}.faidx.log")
    shell:
        "samtools faidx {input} &> {log}"

rule create_sequence_dict:
    input:
        os.path.join(genome_dir, f"{genome_name}.fasta")
    output:
        os.path.join(genome_dir, f"{genome_name}.dict")
    conda:
        "../../envs/gatk4.yaml"
    log:
        os.path.join(genome_dir, f"{genome_name}.dict.log")
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output} &> {log}"

rule build_bwa_index:
    input:
        fasta = os.path.join(genome_dir, f"{genome_name}.fasta")
    output:
        amb = os.path.join(bwa_dir, f"{genome_name}.amb"),
        ann = os.path.join(bwa_dir, f"{genome_name}.ann"),
        bwt = os.path.join(bwa_dir, f"{genome_name}.bwt"),
        pac = os.path.join(bwa_dir, f"{genome_name}.pac"),
        sa  = os.path.join(bwa_dir, f"{genome_name}.sa")
    conda:
        workflow.source_path("../envs/bwa.yaml"),
    log:
        os.path.join(bwa_dir, f"{genome_name}.bwa_index.log")
    params:
        prefix = os.path.join(bwa_dir, genome_name)
    shell:
        "bwa index -p {params.prefix} {input.fasta} &> {log}"

rule copy_and_process_annotations:
    input:
        gtf = config["genome_version_info"]["genome_annotation"]["genome_gtf"],
        gff = config["genome_version_info"]["genome_annotation"]["genome_gff"],
        go  = config["genome_version_info"]["genome_annotation"]["genome_go"]
    output:
        gtf = os.path.join(genome_dir, f"{genome_name}.gtf"),
        gff = os.path.join(genome_dir, f"{genome_name}.gff"),
        go  = os.path.join(genome_dir, f"{genome_name}.go")
    log:
        os.path.join(genome_dir, f"{genome_name}.annotation.log")
    shell:
        """
        cp {input.gtf} {output.gtf} && \
        cp {input.gff} {output.gff} && \
        cp {input.go} {output.go} &> {log}
        """

rule create_bed12:
    input:
        gtf = os.path.join(genome_dir, f"{genome_name}.gtf")
    output:
        bed12 = os.path.join(genome_dir, f"{genome_name}.bed12")
    conda:
        workflow.source_path("../envs/ucsc_gff.yaml"),
    log:
        os.path.join(genome_dir, f"{genome_name}.gtf2bed12.log")
    shell:
        """
        gtfToGenePred -genePredExt -ignoreGroupsWithoutExons {input} {output}.tmp.GenePred && \
        genePredToBed {output}.tmp.GenePred {output} && \
        rm {output}.tmp.GenePred &> {log}
        """

rule create_ref_all:
    input:
        gtf = os.path.join(genome_dir, f"{genome_name}.gtf")
    output:
        ref_all = os.path.join(genome_dir, f"{genome_name}.ref_all.txt")
    conda:
        workflow.source_path("../envs/ucsc_gff.yaml"),
    log:
        os.path.join(genome_dir, f"{genome_name}.gtf2refall.log")
    shell:
        """
        gtfToGenePred -genePredExt -ignoreGroupsWithoutExons {input} {output}.tmp.GenePred && \
        awk 'BEGIN{{OFS="\\t"}} {{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' {output}.tmp.GenePred > {output} && \
        rm {output}.tmp.GenePred &> {log}
        """

rule generate_config_snippet:
    input:
        fasta = os.path.join(genome_dir, f"{genome_name}.fasta"),
        gtf   = os.path.join(genome_dir, f"{genome_name}.gtf"),
        gff   = os.path.join(genome_dir, f"{genome_name}.gff"),
        bed12 = os.path.join(genome_dir, f"{genome_name}.bed12"),
        go    = os.path.join(genome_dir, f"{genome_name}.go"),
        ref_all = os.path.join(genome_dir, f"{genome_name}.ref_all.txt")
    output:
        yaml = os.path.join(genome_dir, f"{genome_name}.reference_config.yaml")
    run:
        species = config["genome_version_info"]["species"]
        rel_genome_path = os.path.join(species, genome_name)
        rel_bwa_path = os.path.join(species, "bwa-mem_v11")
        
        snippet = {
            genome_name: {
                "index": os.path.join(rel_bwa_path, genome_name),
                "genome_fa": os.path.join(rel_genome_path, f"{genome_name}.fasta"),
                "genome_gtf": os.path.join(rel_genome_path, f"{genome_name}.gtf"),
                "genome_gff": os.path.join(rel_genome_path, f"{genome_name}.gff"),
                "bed12": os.path.join(rel_genome_path, f"{genome_name}.bed12"),
                "go_annotation": os.path.join(rel_genome_path, f"{genome_name}.go"),
                "ref_all": os.path.join(rel_genome_path, f"{genome_name}.ref_all.txt")
            }
        }
        
        with open(output.yaml, 'w') as f:
            f.write(f"# Add this to your config/reference.yaml under 'bwa_mem:'\n")
            yaml.dump(snippet, f, default_flow_style=False)
