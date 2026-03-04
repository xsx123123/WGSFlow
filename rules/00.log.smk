#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNAFlow Pipeline - Logging and Runtime Information Module

This module handles the initialization of the logging system and captures
essential runtime information at the start of the RNA-seq analysis pipeline.

Key Responsibilities:
- Initialize the unified Rich/Loguru logger
- Log system and environment details
- Capture pipeline configuration parameters
- Provide comprehensive audit trail for reproducibility

The logging is currently commented out to avoid excessive output during
normal pipeline execution, but can be enabled for debugging purposes.
"""

import sys
import os
import platform
from datetime import datetime

try:
    from snakemake_logger_plugin_rich_loguru import get_analysis_logger
    logger = get_analysis_logger()
    logger_type = "Custom Plugin"

except ImportError:
    try:
        from loguru import logger
        logger_type = "Standard Loguru"
        logger.warning("Custom logger plugin not found. Falling back to standard loguru.")
        
    except ImportError:
        import logging 
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger("Fallback")
        logger.warning("Neither custom plugin nor loguru found. Using built-in logging.")
        logger_type = "Built-in Logging"