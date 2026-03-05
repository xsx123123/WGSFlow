#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Defensive programming: Return an empty list [] if not defined in config
active_callers = config.get('germline_call_variant_type', [])

# Only load variant calling modules if analysis_type is explicitly 'germline'
if config.get('analysis_type') == "germline":
    for caller in active_callers:
        
        # Branch 1: gatk joint call pipeline
        if caller == 'joint_call_variant':
            logger.info("🔧 Loading submodule: GATK Joint Calling")
            include "./sub_rules/gatk_joint_calling.smk"
            
        # Branch 2: single call variant pipeline by gatk 
        elif caller == 'gatk_single_call_variant':
            logger.info("🔧 Loading submodule: GATK Single Calling")
            include "./sub_rules/gatk_call_variant.smk"
            
        # Branch 3: single call variant pipeline by bcftools
        elif caller == 'bcftools_single_call_variant':
            logger.info("🔧 Loading submodule: Bcftools Single Calling")
            include "./sub_rules/bcftools_call_variant.smk"
            
        # Foolproof design: In case the user makes a typo in config.yaml
        else:
            logger.warning(f"⚠️ Unknown Call Variant mode: {caller}. Please check config.yaml!")




