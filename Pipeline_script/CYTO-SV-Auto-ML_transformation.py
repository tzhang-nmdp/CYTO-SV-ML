"""
CYTO-SV-Auto-ML Data Transformation
=====================================
Transforms the combined cohort SV data into model-ready feature matrices.

For each SV type (TRS and nonTRS):
  1. Splits combined data by SV type
  2. Applies database-based labeling thresholds
  3. Transforms read statistics and confidence intervals
  4. Filters low-quality SVs
  5. Computes derived features (read ratios, CI ranges, etc.)
  6. Assigns benchmark labels (-1=artifact, 0=unlabeled, 1=germline, 2=somatic)
  7. Generates summary plots
  8. Outputs two files per SV type:
     - Reduced feature matrix (for modeling)
     - Full data matrix (for prediction on all SVs)

Usage:
  python CYTO-SV-Auto-ML_transformation.py \\
      -s <cohort_name> -o <output_dir> \\
      -t <trs_sv_cutoff> -c <nontrs_sv_cutoff>

Arguments:
  -s : Cohort name (used for file naming)
  -o : Output directory path
  -t : TRS breakpoint distance cutoff (default: 1000 bp)
  -c : NonTRS overlap ratio cutoff (default: 0.9)
"""

import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import subprocess
import random
import pickle
import warnings
import getopt
import statistics
from matplotlib.pyplot import *
from pandas import read_csv
from time import time
from numpy import asarray, sqrt, argmax
from statistics import median
from itertools import *
import sv_dataframe_transform
from sklearn import svm, metrics, datasets, preprocessing
from sklearn.preprocessing import (
    LabelEncoder, label_binarize, StandardScaler,
    OrdinalEncoder, OneHotEncoder
)
from sklearn.ensemble import (
    GradientBoostingRegressor, GradientBoostingClassifier,
    AdaBoostClassifier
)
from sklearn.metrics import (
    mean_squared_error, median_absolute_error, roc_curve, auc,
    accuracy_score, precision_score, recall_score,
    precision_recall_curve, confusion_matrix
)
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer, TransformedTargetRegressor
from sklearn.pipeline import make_pipeline
from supervised.automl import AutoML

warnings.filterwarnings("ignore")

# =============================================================================
# Parse Command-Line Arguments
# =============================================================================
wd = sys.path[0]

# Default cutoffs
trs_sv_cutoff = 1000    # Max breakpoint distance for TRS database matching
nontrs_sv_cutoff = 0.9  # Min overlap ratio for nonTRS database matching

opts, args = getopt.getopt(sys.argv[1:], "s:o:t:c:")
for op, value in opts:
    if op == "-s":
        cohort_name = str(value)
    elif op == "-o":
        outdir = str(value)
    elif op == "-t":
        trs_sv_cutoff = int(value)
    elif op == "-c":
        nontrs_sv_cutoff = float(value)

# SV types to process
sv_type_vector = ['trs', 'nontrs']

# =============================================================================
# Load Combined Cohort SV Data
# =============================================================================
sv_data = pd.read_csv(
    os.path.join(outdir, f'{cohort_name}.sv.all.combine_all'),
    sep='\t', header=0, index_col=None, keep_default_na=False
)
print(sv_data.columns)

# =============================================================================
# Process Each SV Type
# =============================================================================
for sv_type in sv_type_vector:

    # --- Split by SV type ---
    if sv_type == 'trs':
        # TRS: keep BND/TRA SVs (exclude DEL, DUP, INV, INS)
        sv_data_01 = sv_data[
            ~sv_data['sv_type'].isin(['DEL', 'DUP', 'INV', 'INS'])
        ].copy()
        sv_data_01.to_csv(
            os.path.join(outdir, f'{cohort_name}.sv.all.combine_all_trs_o'),
            sep='\t', header=True, index=None
        )

        # Transform: apply database labeling, compute features, filter QC
        sv_data_tf = sv_dataframe_transform.trs_sv_data_transform(
            sv_data_01, trs_sv_cutoff
        )
        sv_data_1 = sv_data_tf[0]   # Reduced feature matrix (for modeling)
        sv_data_1.to_csv(
            os.path.join(outdir, f'{cohort_name}.sv.all.combine_all_trs'),
            sep='\t', header=True, index=None
        )
        sv_data_10 = sv_data_tf[1]  # Full data matrix (for prediction)
        sv_data_10.to_csv(
            os.path.join(outdir, f'{cohort_name}.sv.all.combine_all_trs_all'),
            sep='\t', header=True, index=None
        )
    else:
        # NonTRS: keep DEL, DUP, INV, INS SVs
        sv_data_01 = sv_data[
            sv_data['sv_type'].isin(['DEL', 'DUP', 'INV', 'INS'])
        ].copy()
        sv_data_01.to_csv(
            os.path.join(outdir, f'{cohort_name}.sv.all.combine_all_nontrs_o'),
            sep='\t', header=True, index=None
        )

        # Transform: apply database labeling, compute features, filter QC
        sv_data_tf = sv_dataframe_transform.nontrs_sv_data_transform(
            sv_data_01, nontrs_sv_cutoff
        )
        sv_data_1 = sv_data_tf[0]   # Reduced feature matrix
        sv_data_1.to_csv(
            os.path.join(outdir, f'{cohort_name}.sv.all.combine_all_nontrs'),
            sep='\t', header=True, index=None
        )
        sv_data_10 = sv_data_tf[1]  # Full data matrix
        sv_data_10.to_csv(
            os.path.join(outdir, f'{cohort_name}.sv.all.combine_all_nontrs_all'),
            sep='\t', header=True, index=None
        )

    # --- Generate SV data summary plot ---
    sv_summary_plot = sv_dataframe_transform.sv_data_summary_plot(sv_data_1)
    sv_summary_plot.savefig(
        os.path.join(
            outdir, 'cyto_sv_ml',
            f'{cohort_name}_{sv_type}_sv_summary_plot.pdf'
        ),
        transparent=True
    )
