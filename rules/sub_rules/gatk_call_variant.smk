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
        gvcf = '03.call_variant/single_calling/{sample}.g.vcf.gz',
        tbi = '03.call_variant/single_calling/{sample}.g.vcf.gz.tbi',
    resources:
        **rule_resource(config, 'medium_resource', skip_queue_on_local=True, logger=logger),
    conda:
        workflow.source_path("../envs/gatk4.yaml"),
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
        workflow.source_path("../envs/gatk4.yaml"),
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
        workflow.source_path("../envs/gatk4.yaml"),
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
# ----- rule ----- #

