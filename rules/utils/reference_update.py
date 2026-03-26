#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def resolve_reference_paths(config, targets, base_path=None):
    """
    解析并更新参考基因组的绝对路径,这函数可以让迁移后的分析流程自动更新到对应的应的参考基因组路径

    Args:
        config (dict): Snakemake 的配置字典
        targets: 目标基因组版本，可以是单个字符串 (e.g. "ITAG4.1")
                 也可以是字符串列表 (e.g. ["ITAG4.1", "Lsat_Salinas_v11"])
        base_path (str, optional): 参考基因组的根目录路径。
                                   如果不传，默认尝试从 config['reference_path'] 获取。
    """
    # Import the unified logger
    from snakemake_logger_plugin_rich_loguru import get_analysis_logger
    logger = get_analysis_logger()

    # get analysis reference for reference.yaml
    root_dir = base_path or config.get("reference_path")

    if not root_dir:
        return

    if isinstance(targets, str):
        genome_list = [targets]
    elif isinstance(targets, list):
        genome_list = targets
    else:
        logger.error(f"无法解析 targets 参数类型: {type(targets)}，跳过路径更新。")
        return

    # logger.info(f"正在批量更新参考基因组路径: {genome_list} (Base: {root_dir})")

    # 支持多个配置节：STAR_index 和 bwa_mem
    for config_section in ["STAR_index", "bwa_mem"]:
        index_config = config.get(config_section, {})

        for version in genome_list:
            ref_dict = index_config.get(version)

            if not ref_dict:
                if config_section == "bwa_mem":  # 只在 bwa_mem 时打印警告
                    logger.warning(f"{version} 未在 {config_section} 中定义，跳过。")
                continue

            count = 0
            for key, value in ref_dict.items():
                if isinstance(value, str):
                    if not value.startswith(root_dir) and not value.startswith("/"):
                        full_path = os.path.join(root_dir, value)
                        ref_dict[key] = full_path
                        count += 1

            if count > 0:
                # logger.success(f"[{version}] 成功更新了 {count} 个文件路径")
                pass