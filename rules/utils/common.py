#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGSFlow Pipeline - Common Utility Functions Module

This module provides essential utility functions that are used across multiple
rules in the Whole Genome Sequencing (WGS) analysis pipeline. It handles core functionality including:

Key Components:
- Data delivery orchestration (DataDeliver function)
- Report data collection (ReportData function)
- Sample data directory resolution
- Reference index validation (BWA)
- Genome version compatibility checking

These functions ensure consistent behavior across the pipeline and provide
robust error handling for common scenarios like missing files, invalid
configurations, and path resolution issues.
"""

import os
import glob
import sys
import time
from pathlib import Path
from typing import Dict, Union, List, Callable
from rich import print as rich_print

from datadeliver import (
    qc_clean,
    mapping,
    bcftools_call_variant,
    gatk_call_variant,
    gatk_joint_call_variant,
    manta_sv_variant,
    somatic_variant_paired,
    somatic_variant_single,
)

# Import logger for consistent logging
from snakemake_logger_plugin_rich_loguru import get_analysis_logger

logger = get_analysis_logger()

# Flag to track if the QC warning has already been displayed
_qc_warning_logged = False


def DataDeliver(
    config: Dict = None, samples: Dict = None, all_contrasts: Dict = None
) -> List[str]:
    """
    Main data delivery orchestrator.
    Controls the flow of the pipeline based on 'only_qc' and specific module flags.
    """
    # Initialize config/samples if None
    config = config or {}
    samples = samples or {}

    # ---------------------------------------------------------
    # 0. 初始化基础文件列表 (MD5 Check)
    # ---------------------------------------------------------
    convert_md5_path = config.get("convert_md5", "md5_check")
    data_deliver = [
        "01.qc/md5_check.tsv",
        os.path.join("00.raw_data", convert_md5_path),
        os.path.join("00.raw_data", convert_md5_path, "raw_data_md5.json"),
    ]

    # ---------------------------------------------------------
    # 1. 定义模块类别
    # ---------------------------------------------------------

    # [A] 基础模块：QC 和 Mapping 是所有 WGS 流程的核心
    basic_modules = ["qc_clean", "mapping"]

    # [B] 下游分析模块：由配置动态加载
    downstream_modules = []

    # ---------------------------------------------------------
    # 2. 定义执行包装器 (Wrappers)
    # ---------------------------------------------------------
    def execute_qc_clean(samples, data_deliver):
        return qc_clean(samples, data_deliver)

    def execute_mapping(samples, data_deliver):
        return mapping(samples, data_deliver, config)

    # germline variant analysis (Synchronized with 07.callvraiant.smk)
    def execute_bcftools_single(samples, data_deliver):
        return bcftools_call_variant(samples, data_deliver)

    def execute_gatk_joint(samples, data_deliver):
        return gatk_joint_call_variant(samples, data_deliver)

    def execute_gatk_single(samples, data_deliver):
        return gatk_call_variant(samples, data_deliver)

    # somatic variant analysis
    def execute_samples_only_somatic_variant(samples, data_deliver):
        return somatic_variant_single(samples, data_deliver)

    def execute_paired_group_somatic_variant(_, data_deliver):
        return somatic_variant_paired(all_contrasts, data_deliver)

    def execute_manta_sv_variant(samples, data_deliver):
        return manta_sv_variant(samples, data_deliver)

    # 这里的键名必须与 config.yaml 中的 germline_call_variant_type 列表项严格一致
    module_functions: Dict[str, Callable] = {
        "qc_clean": execute_qc_clean,
        "mapping": execute_mapping,
        "bcftools_single_call_variant": execute_bcftools_single,
        "joint_call_variant": execute_gatk_joint,
        "gatk_single_call_variant": execute_gatk_single,
        "manta_sv_variant": execute_manta_sv_variant,
        "paired_group_somatic_variant": execute_paired_group_somatic_variant,
        "samples_only_somatic_variant": execute_samples_only_somatic_variant,
    }

    # 根据配置解析需要运行的下游模块
    active_callers = config.get("germline_call_variant_type", [])
    if config.get("analysis_type") == "germline" and active_callers:
        for module in active_callers:
            downstream_modules.append(module)

    elif config.get("analysis_type") == "somatic":
        if config.get("whether_paired_group"):
            downstream_modules.append("paired_group_somatic_variant")
        else:
            downstream_modules.append("samples_only_somatic_variant")

    if config.get("call_SV"):
        downstream_modules.append("manta_sv_variant")

    # ---------------------------------------------------------
    # 3. 核心逻辑控制
    # ---------------------------------------------------------
    global _qc_warning_logged

    # 默认开启基础模块
    for module in basic_modules:
        if config.get(module) is not False:
            config[module] = True

    if config.get("only_qc"):
        if not _qc_warning_logged:
            logger.warning(
                "**********************************************************************"
            )
            logger.warning(
                "   [MODE] ONLY QC ENABLED                                            "
            )
            logger.warning(
                "   - Running: Raw QC, Mapping, Coverage Stats                        "
            )
            logger.warning(
                "   - Skipping: Variants (Germline/Somatic), SV Detection             "
            )
            logger.warning(
                "**********************************************************************"
            )
            _qc_warning_logged = True

        for module in downstream_modules:
            config[module] = False
    else:
        for module in downstream_modules:
            if config.get(module) is not False:
                config[module] = True

    # ---------------------------------------------------------
    # 4. 执行并收集输出文件
    # ---------------------------------------------------------
    for module, func in module_functions.items():
        if config.get(module):
            data_deliver = func(samples, data_deliver)

    # ---------------------------------------------------------
    # 5. 调试输出
    # ---------------------------------------------------------
    if config.get("print_target"):
        rich_print("[bold green]Generated Target Files:[/bold green]")
        rich_print(data_deliver)

    return data_deliver


def ReportData(config: dict = None) -> List[str]:
    """
    Collects all files required for generating the final Quarto WGS analysis report.

    This function returns the list of manifest files, summary JSON, and HTML report
    that are needed for the final delivery and reporting stage of the pipeline.

    Args:
        config (dict): Pipeline configuration dictionary

    Returns:
        List[str]: List of report-related file paths, or empty list if reporting
                   is disabled in configuration.
    """
    if config.get("report"):
        return [
            os.path.join(config["data_deliver"], "delivery_manifest.json"),
            os.path.join(config["data_deliver"], "delivery_manifest.md5"),
            os.path.join(config["data_deliver"], "delivery_details.log"),
            os.path.join(config["data_deliver"], "report_data/project_summary.json"),
            os.path.join(
                config["data_deliver"], "report_data", "delivery_manifest.json"
            ),
            os.path.join(
                config["data_deliver"], "report_data", "delivery_manifest.md5"
            ),
            os.path.join(config["data_deliver"], "report_data", "delivery_details.log"),
            os.path.join(config["data_deliver"], "Analysis_Report/index.html"),
        ]
    else:
        return []


def get_sample_data_dir(sample_id: str = None, config: dict = None) -> str:
    """
    Resolves the directory path containing FASTQ files for a given sample ID.

    This function handles two common directory structures:
    1. Sample-specific subdirectories: raw_data/SampleID/
    2. Flat directory structure: raw_data/SampleID_R1.fq.gz

    The function uses fuzzy matching with glob patterns to handle various
    naming conventions and file extensions.

    Args:
        sample_id (str): The sample identifier to search for
        config (dict): Pipeline configuration with 'raw_data_path' key

    Returns:
        str: Absolute path to the directory containing the sample's FASTQ files

    Raises:
        ValueError: If 'raw_data_path' is missing from config
        FileNotFoundError: If no matching files or directories are found
    """

    # Ensure config has this key to prevent errors
    if "raw_data_path" not in config:
        raise ValueError("Config dictionary missing 'raw_data_path' key.")

    # Iterate through all raw data paths in the config
    for base_dir in config["raw_data_path"]:
        # --- Case A: Your previous logic (raw_data/SampleID/xxx.fq) ---
        sample_subdir = os.path.join(base_dir, sample_id)
        if os.path.isdir(sample_subdir):
            return sample_subdir

        # --- Case B: Your current ls result (raw_data/SampleID.R1.fq) ---
        # Use glob fuzzy matching: check if there are files starting with sample_id in the directory
        # pattern similar to: /data/.../00.raw_data/L1MKK1806607-a1*
        pattern = os.path.join(base_dir, f"{sample_id}*")

        # Get the list of matching files
        matching_files = glob.glob(pattern)

        # If matching files are found (and they are files not directories), it means data is at the base_dir level
        if matching_files:
            # Simple filter: ensure they are files (to avoid having a folder named SampleID_tmp by coincidence)
            # As long as one is a file, we consider it found
            if any(os.path.isfile(f) for f in matching_files):
                return base_dir

    # If the loop ends without finding anything
    raise FileNotFoundError(
        f"Could not find data directory or files for {sample_id} in {config['raw_data_path']}"
    )


def get_all_input_dirs(sample_keys=None, config: dict = None) -> list:
    """
    Aggregates unique input directories for all specified samples.

    This function iterates through all sample IDs and resolves their data
    directories using get_sample_data_dir(), then returns a deduplicated list.

    Args:
        sample_keys: Iterable of sample identifiers
        config (dict): Pipeline configuration dictionary

    Returns:
        list: Deduplicated list of input directories containing sample data
    """
    if sample_keys is None or config is None:
        return []
    dir_list = []
    for sample_id in sample_keys:
        dir_list.append(get_sample_data_dir(sample_id, config=config))

    return list(set(dir_list))


def judge_bwa_index(config: dict = None) -> bool:
    """
    Validates BWA-MEM2 index completeness by checking for required index files.

    BWA-MEM2 requires specific index files to be present for proper alignment.
    This function checks if all necessary files exist and returns True if
    the index needs to be rebuilt.

    Args:
        config (dict): Pipeline configuration with 'bwa_mem2' section

    Returns:
        bool: True if index files are missing (needs rebuild), False if complete
    """
    bwa_index = config["bwa_mem2"]["index"]
    bwa_index_files = [
        bwa_index + suffix
        for suffix in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".alt"]
    ]

    return not all(os.path.exists(f) for f in bwa_index_files)


def check_gene_version(config: dict = None, logger=None) -> None:
    """
    Validates that the configured genome version is supported by the pipeline.

    This function ensures that the user-specified genome version matches one of
    the allowed versions defined in the pipeline configuration, preventing
    runtime errors due to unsupported reference genomes.

    Args:
        config (dict): Pipeline configuration with 'Genome_Version' and
                      'can_use_genome_version' keys
        logger: Optional logger instance (falls back to unified logger if None)

    Raises:
        ValueError: If the genome version is not in the allowed list
        KeyError: If required configuration keys are missing
        TypeError: If config is not a valid dictionary
    """
    # Use the provided logger or get the unified logger
    if logger is None:
        from snakemake_logger_plugin_rich_loguru import get_analysis_logger

        logger = get_analysis_logger()

    try:
        version = config["Genome_Version"]
        allowed = config["can_use_genome_version"]

        if version not in allowed:
            logger.error(f"Version mismatch! '{version}' is not in {allowed}")
            raise ValueError(f"Unsupported genome version: {version}")

        logger.info(f"Config check passed: Genome_Version '{version}' is supported.")

    except KeyError as e:
        logger.error(f"Config structure error: Missing key {e}")
        raise
    except TypeError:
        logger.error("Config must be a valid dictionary.")
        raise


# --------------------- #
