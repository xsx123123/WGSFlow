#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - Structural Variant (SV) Detection Module

This module implements structural variant discovery using Manta, a high-performance 
SV caller optimized for germline and somatic analysis.

Key Components:
- call_SV: Discovery of large-scale genomic alterations (Deletions, Duplications, 
  Inversions, and Translocations) using the Manta workflow.

Manta integrates evidence from paired-end read spanning, split-read alignment, 
and localized assembly to provide highly accurate SV calls across a wide range 
of variant sizes.
"""

import os
from loguru import logger

# ----- rule ----- #

rule call_SV:
    """
    Detect Structural Variants (SVs) using the Manta workflow.

    Manta identifies large structural alterations by analyzing discordant 
    read pairs and split-read alignments.

    Workflow Stages:
    1. configManta.py: Configures the SV calling workflow for a specific sample.
    2. runWorkflow.py: Executes the configured workflow using local resources.

    Outputs:
    - diploidSV.vcf.gz: Final set of SV calls for diploid germline samples.
    - runWorkflow.py: The generated execution script for the Manta pipeline.
    
    Parameters:
    - reference: Path to the indexed reference genome FASTA.
    - runDir: Dedicated workspace for Manta's temporary files and results.
    """
    input:
        bam = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam',
        bam_bai = '02.mapping/bwa_mem/MarkDup/{sample}.sort.Dup.bam.bai',
    output:
        run_scripts = '05.SV/manta_SV/{sample}/runWorkflow.py',
        sv_vcf = "05.SV/manta_SV/{sample}/results/variants/diploidSV.vcf.gz",
        sv_vcf_index = "05.SV/manta_SV/{sample}/results/variants/diploidSV.vcf.gz.tbi",
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/manta.yaml"),
    log:
        "logs/05.SV/{sample}_manta.log",
    benchmark:
        "benchmarks/{sample}_manta_benchmark.txt",
    message:
        "Running Manta SV detection for {wildcards.sample}",
    params:
        reference = config['bwa_mem'][config['Genome_Version']]['genome_fa'],
        Workspaces = '05.SV/manta_SV/{sample}',
    threads:
        config['parameter']["threads"]["manta"],
    shell:
        """
        configManta.py --bam {input.bam} \
                       --referenceFasta {params.reference} \
                       --runDir {params.Workspaces} > {log} 2>&1 && \
        {params.Workspaces}/runWorkflow.py -m local -j {threads} >> {log} 2>&1
        """
# ----- END ----- #