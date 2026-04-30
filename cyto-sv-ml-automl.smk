"""
CYTO-SV-ML AutoML Modeling Pipeline (cohort-level)
====================================================
This Snakemake workflow combines all per-sample SV feature matrices and
trains XGBoost classification models to classify SVs into:
  -1: True Artifact (TA)
   1: True Germline (TG)
   2: True Somatic (TS)

Two separate models are trained:
  - TRS model: for translocations / BND SVs
  - NonTRS model: for DEL, DUP, INV SVs

Usage:
  conda activate cyto-sv-ml
  snakemake --cores <N> -s cyto-sv-ml-automl.smk \\
      --config cohort_name=<COHORT_NAME>

Input:
  Per-sample feature files from the preprocess pipeline:
  {main_dir}/out/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp

Output:
  - Combined cohort SV data
  - SV summary plots (TRS/nonTRS)
  - Model performance metrics (precision, recall, F1, accuracy)
  - SHAP feature importance plots
  - Confusion matrices and AUC-ROC curves
  - Per-SV predictions for all samples
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

# ---- Cohort Setup ----
cohort_name = config['cohort_name']

# Read sample list (tab-separated: sample_id, sex)
samples_information = pd.read_csv(
    config['sample_list'], sep='\t', header=None, index_col=False
)
samples_information.columns = ['id', 'sex']
SAMPLES = list(samples_information['id'])
GENDERS = list(samples_information['sex'])

# '@'-delimited sample vector for shell scripts
SAMPLES_vector = '@'.join(str(sm) for sm in SAMPLES)

# ---- Directory Setup ----
MAIN_DIR = config['main_dir']
INPUT_DIR = config['main_dir'] + '/in'
OUTPUT_DIR = config['main_dir'] + '/out'
LOG_DIR = config['main_dir'] + '/out/log'
CYTO_SV_ML_DIR = config['cyto_sv_ml_dir']
SOFTWARE_DIR = config['cyto_sv_ml_dir'] + '/software'
DATABASE_DIR = config['cyto_sv_ml_dir'] + '/SV_database'

# ---- Docker & Caller Config ----
parliament_docker = config['parliament_docker']
chromoseq_docker = config['chromoseq_docker']
parliament2_sv_callers = config['parliament2_sv_callers']
chromoseq_sv_callers = config['chromoseq_sv_callers']
all_callers = chromoseq_sv_callers + parliament2_sv_callers

# ---- SV Size Filter ----
size = int(config['size'])
SIZE_K = round(size / 1000)


# =============================================================================
# Target Rule: final model metrics for both TRS and nonTRS
# =============================================================================
rule all:
    input:
        protected(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_{sv_type}_sv_ml_metrics_sub.csv",
            cohort_name=cohort_name, sv_type=['trs', 'nontrs']
        ))


# =============================================================================
# Checkpoint: Wait for all per-sample preprocessing to complete
# This checkpoint ensures all sample feature files exist before combining
# =============================================================================
checkpoint all_sample_sv_ready:
    input:
        expand(
            OUTPUT_DIR + "/{sample}/{sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp",
            sample=SAMPLES, size_k=SIZE_K
        )
    output:
        pathlib.Path(OUTPUT_DIR + "/log_files/sample_sv_ready.out")
    run:
        shell('echo {SAMPLES_vector} >> {output}')


def check_sample_file(*wildcards):
    """Callback to verify checkpoint completion before proceeding."""
    return checkpoints.all_sample_sv_ready.get().output


# =============================================================================
# Rule: Combine All Sample SV Data
# Concatenates per-sample feature matrices into a single cohort-level file,
# appending sample_id as the last column
# =============================================================================
rule all_sample_sv_combine:
    input:
        check_sample_file
    output:
        protected(expand(
            OUTPUT_DIR + "/{cohort_name}/{cohort_name}.sv.all.combine_all",
            cohort_name=cohort_name
        ))
    shell:
        """
        cat {input} && \
        bash {CYTO_SV_ML_DIR}/Pipeline_script/all_sample_sv_combine.sh \
            {MAIN_DIR} {cohort_name} {input}
        """


# =============================================================================
# Rule: CYTO-SV-ML Model Training & Evaluation
# 1. Data transformation: splits into TRS/nonTRS, applies database labeling,
#    computes derived features (read ratios, CI ranges, etc.)
# 2. Model training: XGBoost via mljar-supervised AutoML with k-fold
#    sub-oversampling to handle class imbalance
# 3. Evaluation: confusion matrix, AUC-ROC, SHAP importance, metrics CSV
# =============================================================================
rule cyto_sv_ml:
    input:
        expand(
            OUTPUT_DIR + "/{cohort_name}/{cohort_name}.sv.all.combine_all",
            cohort_name=cohort_name
        )
    output:
        # Summary plots
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_sv_summary_plot.pdf",
            cohort_name=cohort_name, k=0),
            category="sv data summary", subcategory="data",
            labels={"data name": "sv type distribution", "sv type": "trs", "data type": "plot"}
        ),
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_sv_summary_plot.pdf",
            cohort_name=cohort_name, k=0),
            category="sv data summary", subcategory="data",
            labels={"data name": "sv type distribution", "sv type": "nontrs", "data type": "plot"}
        ),
        # Model metrics tables
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_sv_ml_metrics_sub.csv",
            cohort_name=cohort_name),
            category="sv model summary", subcategory="model",
            labels={"data name": "model performance metrics", "sv type": "trs", "data type": "table"}
        ),
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_sv_ml_metrics_sub.csv",
            cohort_name=cohort_name),
            category="sv model summary", subcategory="model",
            labels={"data name": "model performance metrics", "sv type": "nontrs", "data type": "table"}
        ),
        # SHAP feature importance plots
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_{k}_ts_EXP/learner_fold_0_shap_summary.png",
            cohort_name=cohort_name, k=0),
            category="sv model summary", subcategory="model",
            labels={"data name": "shap feature importance", "sv type": "trs", "data type": "plot"}
        ),
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_{k}_ts_EXP/learner_fold_0_shap_summary.png",
            cohort_name=cohort_name, k=0),
            category="sv model summary", subcategory="model",
            labels={"data name": "shap feature importance", "sv type": "nontrs", "data type": "plot"}
        ),
        # Confusion matrices
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_{k}_ts_model_confusion_matrix.pdf",
            cohort_name=cohort_name, k=0),
            category="sv model summary", subcategory="model",
            labels={"data name": "model confusion matrix", "sv type": "trs", "data type": "plot"}
        ),
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_{k}_ts_model_confusion_matrix.pdf",
            cohort_name=cohort_name, k=0),
            category="sv model summary", subcategory="model",
            labels={"data name": "model confusion matrix", "sv type": "nontrs", "data type": "plot"}
        ),
        # AUC-ROC curves
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_{k}_ts_model_aucroc_curve.pdf",
            cohort_name=cohort_name, k=0),
            category="sv model summary", subcategory="model",
            labels={"data name": "model aucroc curve", "sv type": "trs", "data type": "plot"}
        ),
        report(expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_{k}_ts_model_aucroc_curve.pdf",
            cohort_name=cohort_name, k=0),
            category="sv model summary", subcategory="model",
            labels={"data name": "model aucroc curve", "sv type": "nontrs", "data type": "plot"}
        )
    params:
        kfs = config['kfolds'],
        trs_sv_cutoff = config['trs_sv_cutoff'],
        nontrs_sv_cutoff = config['nontrs_sv_cutoff'],
        sv_feature_metrics_index = config['sv_feature_metrics_index_file'],
        py39_dir = config['py39_dir']
    shell:
        """
        sudo mkdir -p {OUTPUT_DIR}/{cohort_name}/cyto_sv_ml && \
        {params.py39_dir}/python \
            {CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML_transformation.py \
            -s {cohort_name} \
            -o {OUTPUT_DIR}/{cohort_name} \
            -t {params.trs_sv_cutoff} \
            -c {params.nontrs_sv_cutoff} && \
        {params.py39_dir}/python \
            {CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML_modelling.py \
            -s {cohort_name} \
            -o {OUTPUT_DIR}/{cohort_name} \
            -x {params.sv_feature_metrics_index} \
            -k {params.kfs}
        """
