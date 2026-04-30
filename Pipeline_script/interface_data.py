"""
interface_data.py
==================
Generates SV predictions for the R Shiny web portal.

Loads the trained AutoML models (last k-fold) and preprocessing pipelines,
predicts SV classes for all SVs, then prepares a combined output file
with predictions, read statistics, and database annotations.

Steps:
  1. Load trained TRS and nonTRS AutoML models + preprocessing pipelines
  2. Read transformed SV feature matrices
  3. Apply preprocessing and predict SV classes
  4. Rename columns to human-readable names for the web portal
  5. Compute derived read statistics (total ref/alt reads, read diff)
  6. Annotate which databases each SV matches
  7. Combine TRS + nonTRS predictions into a single file
  8. Compute max prediction probability across classes

Usage:
  python interface_data.py -t <trs_model_prefix> -n <nontrs_model_prefix> -i <data_prefix>

Arguments:
  -t : TRS model prefix (reads {prefix}_ts_EXP/ and {prefix}_tf.pickle)
  -n : nonTRS model prefix (same pattern)
  -i : Data file prefix (reads {prefix}.sv.all.combine_all_trs etc.)

Output:
  {data_prefix}.trs_pred     - TRS predictions with metadata
  {data_prefix}.nontrs_pred  - nonTRS predictions with metadata
  {data_prefix}.all_pred     - Combined predictions for Shiny app
"""

import sys, getopt, os
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys, subprocess, random, pickle, warnings, getopt, statistics
from matplotlib.pyplot import *
from pandas import read_csv
from time import time
from numpy import asarray, sqrt, argmax
from statistics import median
from itertools import *
import sv_dataframe_transform
from sklearn import svm, metrics, datasets, preprocessing
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler, OrdinalEncoder, OneHotEncoder 
from sklearn.ensemble import GradientBoostingRegressor,GradientBoostingClassifier, AdaBoostClassifier
from sklearn.metrics import mean_squared_error, median_absolute_error, roc_curve, auc, accuracy_score, precision_score, recall_score, precision_recall_curve, confusion_matrix
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer, TransformedTargetRegressor
from sklearn.pipeline import make_pipeline
from supervised.automl import AutoML 
warnings.filterwarnings("ignore")

############################################################################################################################################################# 

#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"t:n:i:")
inFile = ""

for op, value in opts:
	if op == "-t":
	    trs_dir = str(value)
	if op == "-n":
	    nontrs_dir = str(value)
	if op == "-i":
	    inFile =  str(value)

if inFile == "":
	print("-i invalid")
	sys.exit()
############################################################################################################################################################# 
      
# load the optimal model
automl_trs = AutoML(mode="Explain", n_jobs= 6, results_path=trs_dir+'_ts_EXP')
automl_nontrs = AutoML(mode="Explain", n_jobs= 6, results_path=nontrs_dir+'_ts_EXP')

# load the training data tranformer
model_trs=pickle.load(open(trs_dir+'_tf.pickle', 'rb'))
model_nontrs=pickle.load(open(nontrs_dir+'_tf.pickle', 'rb')) 

# read in all WGS SV data for model prediction
X_trs=pd.read_csv(inFile+'.sv.all.combine_all_trs',sep='\t',header=0, index_col=None, keep_default_na=False)
X_nontrs=pd.read_csv(inFile+'.sv.all.combine_all_nontrs',sep='\t',header=0, index_col=None, keep_default_na=False)
X_trs_all=pd.read_csv(inFile+'.sv.all.combine_all_trs_all',sep='\t',header=0, index_col=None, keep_default_na=False)
X_nontrs_all=pd.read_csv(inFile+'.sv.all.combine_all_nontrs_all',sep='\t',header=0, index_col=None, keep_default_na=False)

# transform all the original SV data
s_trs=model_trs.transform(X_trs)
s_nontrs=model_nontrs.transform(X_nontrs)

# predict all the transformed SV data
trs_label=automl_trs.predict_all(s_trs)
nontrs_label=automl_nontrs.predict_all(s_nontrs)
trs_label=trs_label.rename(columns={'label':'predict_label'})
nontrs_label=nontrs_label.rename(columns={'label':'predict_label'})
print(X_trs_all.shape)
print(X_nontrs_all.shape)
print(trs_label.shape)
print(nontrs_label.shape)

# attach the prediction label to trs/nontrs the original SV data
X_trs_all=pd.concat([X_trs_all,trs_label],axis=1)
X_nontrs_all=pd.concat([X_nontrs_all,nontrs_label],axis=1)

# re-organize in trs/nontrs data
col_rename={'SVTYPE':'svtype','prediction_-1':'prediction_TA','prediction_1':'prediction_TG','prediction_2':'prediction_TS','sv_chr': 'CHROM', 'sv_chr2':'CHR2', 'sv_bp_st':'POS', 'sv_bp_end':'END', 'sv_bp_st_ci_range':'cipos_range', 'sv_bp_end_ci_range':'ciend_range', 'sv_bp_st_ci0':'ci_pos0', 'sv_bp_st_ci1':'ci_pos1', 'sv_bp_end_ci0':'sv_bp_end_POS', 'sv_bp_end_ci1':'sv_bp_end_END', 'sv_read_ratio':'read_ratio', 'PR_ref':'sv_PR_ref', 'SR_ref':'sv_SR_ref','PR_alt':'sv_PR_alt', 'SR_alt':'sv_SR_alt', 'uwstl_s':'chromoseq_s', '1000_g':'1000G_g', '1000_gall':'1000G_all_g', 'cosmic_s':'COMSIC_s', 'cytoatlas_s':'CytoAtlas_s', 'gnomad_g2ab':'gnomAD_nonsingleton_g', 'gnomad_gall':'gnomAD_all_g', 'gnomad_g':'gnomAD_g', 'gnomad_qc':'gnomAD_qc', 'dgv_g':'DGV_g', 'control_g':'donor_g', 'control_gall':'donor_all_g'}
X_trs_all=X_trs_all.rename(columns=col_rename)
X_nontrs_all=X_nontrs_all.rename(columns=col_rename)
X_trs_all=X_trs_all.assign(sv_read_r=0,sv_read_a=0,sv_read_diff=0,sv_database=0)
X_trs_all.loc[:,'sv_read_r']=X_trs_all.loc[:,['sv_PR_ref','sv_SR_ref']].apply(lambda x: sum(x), axis=1) 
X_trs_all.loc[:,'sv_read_a']=X_trs_all.loc[:,['sv_PR_alt','sv_SR_alt']].apply(lambda x: sum(x), axis=1) 
X_trs_all.loc[:,'sv_read_diff']=X_trs_all.loc[:,['sv_read_r','sv_read_a']].apply(lambda x: x[1]-x[0], axis=1) 

X_nontrs_all=X_nontrs_all.assign(sv_read_r=0,sv_read_a=0,sv_read_diff=0,sv_database=0)
X_nontrs_all.loc[:,'sv_read_r']=X_nontrs_all.loc[:,['sv_PR_ref','sv_SR_ref']].apply(lambda x: sum(x), axis=1) 
X_nontrs_all.loc[:,'sv_read_a']=X_nontrs_all.loc[:,['sv_PR_alt','sv_SR_alt']].apply(lambda x: sum(x), axis=1) 
X_nontrs_all.loc[:,'sv_read_diff']=X_nontrs_all.loc[:,['sv_read_r','sv_read_a']].apply(lambda x: x[1]-x[0], axis=1) 

# label the sv benchmark database
def sv_db_an(sv_data,sv_type):
    if sv_type=='trs':
        sv_db_idx=np.where(sv_data=='YVID')
    elif sv_type=='nontrs':
        sv_data[sv_data=='Not_in_database']=0
        sv_db_idx=np.where(sv_data.astype(float)>=0.9)
    if sv_db_idx!=[]:
        return '&'.join(str(w) for w in sv_data.index[sv_db_idx])
    else:
        return 'NA'
col_database=['chromoseq_s', '1000G_g', '1000G_all_g', 'COMSIC_s', 'CytoAtlas_s', 'gnomAD_nonsingleton_g', 'gnomAD_all_g', 'gnomAD_g', 'gnomAD_qc', 'DGV_g', 'donor_g', 'donor_all_g']
X_trs_all['sv_database']=X_trs_all.loc[:,col_database].apply(lambda x: sv_db_an(x,'trs'),axis=1) 
X_nontrs_all['sv_database']=X_nontrs_all.loc[:,col_database].apply(lambda x: sv_db_an(x,'nontrs'),axis=1) 

# mocked unqiue id
tmp_id=['sv_id','sample_id']
X_trs_all['tmp_id']=0
X_nontrs_all['tmp_id']=0
X_trs_all.loc[:,'tmp_id']=X_trs_all.loc[:,tmp_id].apply(lambda x: ':'.join(str(w) for w in x), axis=1)      
X_nontrs_all.loc[:,'tmp_id']=X_nontrs_all.loc[:,tmp_id].apply(lambda x: ':'.join(str(w) for w in x), axis=1)

# save the temp files
X_trs_all.to_csv(inFile+'.trs_pred',sep='\t',index=False,header=True)
X_nontrs_all.to_csv(inFile+'.nontrs_pred',sep='\t',index=False,header=True)

# re-organize trs/nontrs data
col_keep=['sv_type', 'sv_chr1', 'sv_chr2', 'sv_read_r', 'sv_read_a','sv_PR_ref','sv_SR_ref','sv_PR_alt','sv_SR_alt', 'read_ratio', 'read_ratio_diff', 'sv_start_bp', 'sv_end_bp', 'ci_pos0', 'ci_pos1', 'sv_bp_end_POS', 'sv_bp_end_END', 'cipos_range', 'ciend_range', 'sv_bp_st_cc_v1', 'sv_bp_end_cc_v1', 'sv_database', 'predict_label', 'prediction_TA', 'prediction_TG', 'prediction_TS', 'label'] 
col_drop=[str(w) for w in X_trs_all.columns if w not in col_keep]
print(col_drop)
X_trs_all=X_trs_all.drop(col_drop,axis=1)
print(X_trs_all.columns)
col_drop=[str(w) for w in X_nontrs_all.columns if w not in col_keep]
print(col_drop)
X_nontrs_all=X_nontrs_all.drop(col_drop,axis=1)
print(X_nontrs_all.columns)

# combined trs/nontrs data
X_all=pd.concat([X_trs_all,X_nontrs_all],axis=0)
X_all=X_all.rename(columns={'sv_chr1': 'sv_chr','sv_start_bp':'sv_bp_st','sv_end_bp':'sv_bp_end', 'read_ratio':'sv_read_ratio','read_ratio_diff':'sv_read_ratio_diff', 'ci_pos0':'sv_bp_st_ci0', 'ci_pos1':'sv_bp_st_ci1', 'sv_bp_end_POS':'sv_bp_end_ci0', 'sv_bp_end_END':'sv_bp_end_ci1','cipos_range':'sv_bp_st_ci_range', 'ciend_range':'sv_bp_end_ci_range'})  

# attach the prediction label to all the original SV data
X_all['predict_max']=0
X_all['predict_max']=X_all.loc[:,['prediction_TA', 'prediction_TG', 'prediction_TS']].apply(lambda x: max(x), axis=1) 
print(X_all.iloc[0:5,:])
print(X_all.columns)

# save all the SV data
X_all.to_csv(inFile+'.all_pred',sep='\t',index=False,header=True)
