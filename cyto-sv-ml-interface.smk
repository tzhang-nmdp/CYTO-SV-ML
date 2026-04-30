"""
CYTO-SV-ML Interface Pipeline (cohort-level)
==============================================
This Snakemake workflow builds an R Shiny web portal for interactive
exploration of SV classification results. It:
  1. Loads the trained AutoML models (last k-fold)
  2. Predicts SV classes for all samples
  3. Prepares data for the R Shiny app
  4. Builds a Docker image containing the Shiny app

Usage:
  conda activate cyto-sv-ml
  snakemake --cores <N> -s cyto-sv-ml-interface.smk \\
      --config cohort_name=<COHORT_NAME>

  # Then run the Docker container:
  sudo docker run -d -p 8000:80 cyto-sv-ml-app:<COHORT_NAME>
  # Open http://localhost:8000/ in a browser

Note: ~2 GB disk space is required to build and run the Docker image.
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

# Read sample list
samples_information = pd.read_csv(
    config['sample_list'], sep='\t', header=None, index_col=False
)
samples_information.columns = ['id', 'sex']
SAMPLES = list(samples_information['id'])
GENDERS = list(samples_information['sex'])
SAMPLES_vector = '@'.join(str(sm) for sm in SAMPLES)
print(SAMPLES_vector)

# ---- Directory Setup ----
MAIN_DIR = config['main_dir']
INPUT_DIR = config['main_dir'] + '/in'
OUTPUT_DIR = config['main_dir'] + '/out'
LOG_DIR = config['main_dir'] + '/out/log'
CYTO_SV_ML_DIR = config['cyto_sv_ml_dir']
SOFTWARE_DIR = config['cyto_sv_ml_dir'] + '/software'
DATABASE_DIR = config['cyto_sv_ml_dir'] + '/SV_database'
CYTO_BAND_FILE = config['cyto_band_file']


# =============================================================================
# Rule: Build R Shiny Docker Interface
# - Uses the last k-fold model (k=kfolds-1) for predictions
# - Runs interface_data.py to generate predictions for all SVs
# - Runs interface_data.R to prepare RData for the Shiny app
# - Builds Docker image with the Shiny app
# =============================================================================
rule interface_docker:
    input:
        # Requires all k-fold model directories for both TRS and nonTRS
        expand(
            OUTPUT_DIR + "/{cohort_name}/cyto_sv_ml/{cohort_name}_{sv_type}_{k}_ts_EXP",
            cohort_name=cohort_name,
            k=[0, 1, 2, 3, 4],
            sv_type=['trs', 'nontrs']
        )
    output:
        # Combined prediction file for all SVs
        expand(
            OUTPUT_DIR + "/{cohort_name}/{cohort_name}.all_pred",
            cohort_name=cohort_name
        ),
        # RData file for the Shiny app
        CYTO_SV_ML_DIR + "/docker/cyto-sv-ml/data/cyto_sv_ml.RData"
    params:
        kfolds = config['kfolds']
    shell:
        """
        echo {input} && \
        bash {CYTO_SV_ML_DIR}/Pipeline_script/interface_docker.sh \
            {MAIN_DIR} {CYTO_SV_ML_DIR} {cohort_name} \
            {CYTO_BAND_FILE} {params.kfolds}
        """
