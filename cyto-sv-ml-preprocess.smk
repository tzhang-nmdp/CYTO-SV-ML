"""
CYTO-SV-ML Preprocessing Pipeline (per-sample)
================================================
This Snakemake workflow processes a single sample through:
  1. SV calling via ChromoSeq (Manta + ichorCNV) and Parliament2 (Breakdancer,
     CNVnator, Delly deletion/duplication/inversion)
  2. VCF transformation and size filtering
  3. SV merging across callers using SURVIVOR
  4. SV genotyping via SVTyper
  5. Breakpoint sequence complexity analysis (SeqComplex + Komplexity)
  6. SV database annotation (12 public/internal databases)
  7. VCF info extraction and feature combination

Usage:
  conda activate cyto-sv-ml
  snakemake --cores <N> -s cyto-sv-ml-preprocess.smk --use-conda \\
      --config sample=<SAMPLE_ID> gender=<Male|Female>

Input:
  {main_dir}/in/{sample}.cram   (for ChromoSeq)
  {main_dir}/in/{sample}.bam    (for Parliament2)

Output:
  {main_dir}/out/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp
"""

import os
import sys
import pandas as pd
import numpy as np
import pathlib
import snakemake.io
from snakemake.utils import validate
from typing import Dict, Union, List

# ---- Load Configuration ----
configfile: "config.yaml"

# ---- Sample & Gender (single-sample mode) ----
SAMPLES = config['sample']
GENDERS = config['gender']

# ---- Directory Setup ----
MAIN_DIR = config['main_dir']
INPUT_DIR = config['main_dir'] + '/in'
OUTPUT_DIR = config['main_dir'] + '/out'
LOG_DIR = config['main_dir'] + '/out/log'
CYTO_SV_ML_DIR = config['cyto_sv_ml_dir']
SOFTWARE_DIR = config['cyto_sv_ml_dir'] + '/software'
DATABASE_DIR = config['cyto_sv_ml_dir'] + '/SV_database'

# ---- Docker Images ----
parliament_docker = config['parliament_docker']
chromoseq_docker = config['chromoseq_docker']

# ---- SV Caller Lists ----
parliament2_sv_callers = config['parliament2_sv_callers']
chromoseq_sv_callers = config['chromoseq_sv_callers']
all_callers = chromoseq_sv_callers + parliament2_sv_callers  # 7 callers total

# ---- SV Database List ----
SV_DB = config['sv_db']

# ---- SV Size Filter ----
size = int(config['size'])       # Minimum SV size in bp (default: 10000)
SIZE_K = round(size / 1000)      # Size in kb for file naming (default: 10)


# =============================================================================
# Target Rule: final combined feature file per sample
# =============================================================================
rule all:
    input:
        protected(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp",
            sample=SAMPLES, size_k=SIZE_K
        ))


# =============================================================================
# Rule 1: Run ChromoSeq Docker (Manta + ichorCNV)
# Calls SVs from CRAM input using the ChromoSeq clinical pipeline
# =============================================================================
rule chromoseq_sv:
    input:
        sample_cram = expand(INPUT_DIR + "/{sample}.cram", sample=SAMPLES)
    output:
        sample_vcf = protected(expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf",
            sample=SAMPLES, sv_caller=chromoseq_sv_callers
        ))
    threads: 8
    params:
        chromoseq_docker = chromoseq_docker,
        sm = SAMPLES,
        gd = GENDERS
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/run_chromoseq.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.chromoseq_docker} \
            {params.sm} {params.gd}
        """


# =============================================================================
# Rule 2: Run Parliament2 Docker (Breakdancer, CNVnator, Delly)
# Calls SVs from BAM input using the Parliament2 multi-caller pipeline
# =============================================================================
rule parliament2_sv:
    input:
        sample_bam = expand(INPUT_DIR + "/{sample}.bam", sample=SAMPLES)
    output:
        sample_vcf = protected(expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf",
            sample=SAMPLES, sv_caller=parliament2_sv_callers
        ))
    threads: 8
    params:
        parliament_docker = parliament_docker,
        sm = SAMPLES
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/run_parliament2.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.parliament_docker} {params.sm}
        """


# =============================================================================
# Rule 3: SV VCF Transformation
# - Standardizes SV IDs across callers
# - Filters SVs by minimum size
# - Splits into TRS (translocation) and nonTRS (DEL/DUP/INV) VCFs
# - Extracts simplified VCF info
# =============================================================================
rule sv_vcf_tf:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf",
            sample=SAMPLES, sv_caller=all_callers
        )
    output:
        expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf",
            sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K,
            sv_type=['trs', 'nontrs']
        ),
        temp(expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.sv_info.sim",
            sample=SAMPLES, size_k=SIZE_K, sv_caller=all_callers
        ))
    params:
        sm = SAMPLES,
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_vcf_tf.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {size}
        """


# =============================================================================
# Rule 4: SV Merging with SURVIVOR
# - Merges nonTRS SVs across callers (breakpoint distance <= 1000bp)
# - Copies TRS SVs from Manta (only caller producing BND)
# - Merges all SVs into a consolidated VCF
# - Extracts caller support (SUPP) and ID mappings
# =============================================================================
rule svmerge_qc:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf",
            sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K,
            sv_type=['trs', 'nontrs']
        )
    output:
        temp(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.{sv_type}_tf.all",
            sample=SAMPLES, size_k=SIZE_K, sv_type=['trs', 'nontrs']
        )),
        temp(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.{sv_type}",
            sample=SAMPLES, size_k=SIZE_K, sv_type=['trs', 'nontrs']
        )),
        temp(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all",
            size_k=SIZE_K, sample=SAMPLES
        )),
        protected(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping",
            size_k=SIZE_K, sample=SAMPLES
        ))
    params:
        sm = SAMPLES
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/svmerge_qc.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {SIZE_K}
        """


# =============================================================================
# Rule 5: SVTyper Genotyping
# - Runs SVTyper on each caller's VCF (TRS + nonTRS) for read-level QC
# - Extracts paired-read (PR) and split-read (SR) support
# - Combines svtyped VCFs per caller
# =============================================================================
rule svtyper_qc:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf",
            sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K,
            sv_type=['trs', 'nontrs']
        ),
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.{sv_type}_tf.all",
            sample=SAMPLES, size_k=SIZE_K, sv_type=['trs', 'nontrs']
        )
    output:
        protected(expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.{size_k}k.all.svtyped.vcf",
            sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K
        )),
        temp(expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.{size_k}k.all.svtyped.vcf.sv_info.sim",
            sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K
        )),
        protected(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.all.svtyped.vcf",
            sample=SAMPLES, size_k=SIZE_K
        )),
        temp(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.all.svtyped.vcf.sv_info.sim",
            sample=SAMPLES, size_k=SIZE_K
        ))
    params:
        sm = SAMPLES,
        sv_caller = '@'.join(str(sc) for sc in all_callers)  # '@'-delimited caller list
    conda:
        "py27.yaml"
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/svtyper_qc.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {params.sv_caller} {SIZE_K}
        """


# =============================================================================
# Rule 6: Breakpoint Sequence Complexity
# - Extracts 300bp flanking sequences around each SV breakpoint
# - Computes complexity metrics using SeqComplex (24 features) and
#   Komplexity (1 feature) for both start and end breakpoints
# =============================================================================
rule sv_seq_complex:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all",
            sample=SAMPLES, size_k=SIZE_K
        )
    output:
        temp(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.bed.bpst_bpend.kz.index_complex",
            sample=SAMPLES, size_k=SIZE_K
        ))
    params:
        sm = SAMPLES,
        py27_dir = config['py27_dir']
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_seq_complex.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {params.py27_dir} {SIZE_K}
        """


# =============================================================================
# Rule 7: SV Database Annotation
# - Annotates each SV against 12 public/internal SV databases
# - TRS SVs: matched by breakpoint distance
# - NonTRS SVs: matched by region overlap ratio
# - Produces per-SV annotation labels for downstream labeling
# =============================================================================
rule sv_database_ann:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.{sv_type}",
            sample=SAMPLES, size_k=SIZE_K, sv_type=['trs', 'nontrs']
        ),
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping",
            size_k=SIZE_K, sample=SAMPLES
        )
    output:
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping.all_anno",
            sample=SAMPLES, size_k=SIZE_K
        )
    params:
        sm = SAMPLES,
        sv_db = '@'.join(str(sd) for sd in SV_DB),  # '@'-delimited database list
        py27_dir = config['py27_dir']
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_database_ann.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {params.sv_db} \
            {params.py27_dir} {SIZE_K}
        """


# =============================================================================
# Rule 8: SV VCF Info Extraction
# - Extracts read-level features (PR, SR, GT, CN, BND_DEPTH, etc.)
#   from svtyped VCFs and original caller VCFs
# - Maps features back to consolidated SV IDs
# =============================================================================
rule sv_info_extract:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.{size_k}k.all.svtyped.vcf.sv_info.sim",
            sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K
        ),
        expand(
            OUTPUT_DIR + "/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.sv_info.sim",
            sample=SAMPLES, size_k=SIZE_K, sv_caller=all_callers
        ),
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping",
            size_k=SIZE_K, sample=SAMPLES
        )
    output:
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping.all_info",
            sample=SAMPLES, size_k=SIZE_K
        )
    params:
        sm = SAMPLES,
        sv_caller = '@'.join(str(sc) for sc in all_callers)
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_info_extract.sh \
            {MAIN_DIR} {params.sm} {params.sv_caller} {SIZE_K}
        """


# =============================================================================
# Rule 9: Combine All SV Features
# - Joins database annotations, VCF info features, sequence complexity,
#   and caller support count into a single feature matrix per sample
# - This is the final output used for cohort-level modeling
# =============================================================================
rule sv_all_combine:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping.{feature}",
            sample=SAMPLES, size_k=SIZE_K, feature=['all_anno', 'all_info']
        ),
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.bed.bpst_bpend.kz.index_complex",
            sample=SAMPLES, size_k=SIZE_K
        )
    output:
        report(expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp",
            sample=SAMPLES, size_k=SIZE_K
        ))
    params:
        sm = SAMPLES
    shell:
        """
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_all_combine.sh \
            {MAIN_DIR} {params.sm} {SIZE_K}
        """
