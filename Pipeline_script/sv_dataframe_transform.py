"""
sv_dataframe_transform.py
=========================
Core data transformation module for CYTO-SV-ML.

Provides three functions:
  - trs_sv_data_transform:   Transform TRS (translocation/BND) SV data
  - nontrs_sv_data_transform: Transform nonTRS (DEL/DUP/INV) SV data
  - sv_data_summary_plot:    Generate stacked bar plot of SV class distribution

Transformation steps (both TRS and nonTRS):
  1. Shuffle data randomly (seed=0 for reproducibility)
  2. Convert database annotations to binary labels (YVID = matched)
  3. Parse confidence intervals (CIPOS/CIEND) into numeric columns
  4. Harmonize read support fields (PR/SR) across different caller formats
  5. Filter low-quality SVs
  6. Compute derived features (read ratios, CI ranges, depth ratios)
  7. Assign benchmark labels: -1=artifact, 0=unlabeled, 1=germline, 2=somatic
  8. Select model features and return (reduced_df, full_df)

Labeling logic:
  - Germline (1): matched in 1000G, gnomAD (AF>=0.05), or DGV
  - Artifact (-1): matched in gnomAD QC-fail, normal donor controls, or centromere
  - Somatic (2): matched in COSMIC/CytoAtlas/ChromoSeq, NOT in germline/control DBs
  - Unlabeled (0): no exclusive match to any category
  - Conflicts (matched in multiple categories) are left unlabeled (0)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from itertools import repeat
import os, sys, getopt


def trs_sv_data_transform(sv_data_trs, trs_sv_cutoff):
    """Transform TRS (translocation/BND) SV data for ML modeling.

    Args:
        sv_data_trs (pd.DataFrame): Raw TRS SV data with all features.
        trs_sv_cutoff (int): Max breakpoint distance (bp) for database matching.
            SVs with distance <= cutoff are labeled 'YVID' (Yes, Validated In Database).

    Returns:
        tuple: (sv_data_trs_2, sv_data_trs)
            - sv_data_trs_2: Reduced feature matrix (model features + label)
            - sv_data_trs: Full data matrix with all columns + label
    """
    # Shuffle data randomly for unbiased training
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs.shape[0])), k=sv_data_trs.shape[0]) 
    sv_data_trs=sv_data_trs.reset_index(drop=True)
    sv_data_trs=sv_data_trs.reindex(new_index)
    sv_data_trs=sv_data_trs.reset_index(drop=True)   
    print(sv_data_trs.shape) 
    
    # transform database label
    sv_db_list=['1000_gall','1000_g','gnomad_gall','gnomad_g','control_gall', 'control_g', 'cytoatlas_s','cosmic_s','gnomad_qc', 'dgv_g', 'centromere_qc','uwstl_s' ]
    for sv_db in sv_db_list:
        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db].isnull()) ],sv_db]=trs_sv_cutoff*10000
        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db]=='NAN')],sv_db]=trs_sv_cutoff*10000
        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db]=='Not_in_database')],sv_db]=trs_sv_cutoff*10000   
        
    sv_data_trs.loc[sv_data_trs['1000_gall'].astype(int)<=trs_sv_cutoff,'1000_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['1000_g'].astype(int)<=trs_sv_cutoff,'1000_g']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_gall'].astype(int)<=trs_sv_cutoff,'gnomad_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_g'].astype(int)<=trs_sv_cutoff,'gnomad_g']="YVID"
    sv_data_trs.loc[sv_data_trs['control_gall'].astype(int)<=trs_sv_cutoff,'control_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['control_g'].astype(int)<=trs_sv_cutoff,'control_g']="YVID"
    sv_data_trs.loc[sv_data_trs['cytoatlas_s'].astype(int)<=trs_sv_cutoff,'cytoatlas_s']="YVID"
    sv_data_trs.loc[sv_data_trs['cosmic_s'].astype(int)<=trs_sv_cutoff,'cosmic_s']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_qc'].astype(int)<=trs_sv_cutoff,'gnomad_qc']="YVID"
    sv_data_trs.loc[sv_data_trs['dgv_g'].astype(int)<=trs_sv_cutoff,'dgv_g']="YVID"
    sv_data_trs.loc[sv_data_trs['centromere_qc'].astype(int)<=trs_sv_cutoff,'centromere_qc']="YVID"
    sv_data_trs.loc[sv_data_trs['uwstl_s'].astype(int)<=trs_sv_cutoff,'uwstl_s']="YVID"
    
    # transform sv_info variables
    sv_data_trs['ci_pos0']=0
    sv_data_trs['ci_pos1']=0
    sv_data_trs['ci_end']=0
    sv_data_trs['PR_ref']=0
    sv_data_trs['PR_alt']=0
    sv_data_trs['SR_ref']=0
    sv_data_trs['SR_alt']=0
    sv_data_trs['ci_pos0']=sv_data_trs['CIPOS'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['ci_pos1']=sv_data_trs['CIPOS'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_trs['sv_bp_end_POS']=sv_data_trs['CIEND'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['sv_bp_end_END']=sv_data_trs['CIEND'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0) 

    # read statistics transformation
    sv_data_trs.loc[(sv_data_trs['RP'].isnull()) | (sv_data_trs['RP']=='NA') | (sv_data_trs['RP']=='.'),'RP']=0  
    sv_data_trs.loc[(sv_data_trs['AP'].isnull()) | (sv_data_trs['AP']=='NA') | (sv_data_trs['AP']=='.'),'AP']=0  
    sv_data_trs.loc[(sv_data_trs['RS'].isnull()) | (sv_data_trs['RS']=='NA') | (sv_data_trs['RS']=='.'),'RS']=0          
    sv_data_trs.loc[(sv_data_trs['AS'].isnull()) | (sv_data_trs['AS']=='NA') | (sv_data_trs['AS']=='.'),'AS']=0    
    sv_data_trs.loc[(sv_data_trs['DR'].isnull()) | (sv_data_trs['DR']=='NA') | (sv_data_trs['DR']=='.'),'DR']=0  
    sv_data_trs.loc[(sv_data_trs['DV'].isnull()) | (sv_data_trs['DV']=='NA') | (sv_data_trs['DV']=='.'),'DV']=0  
    sv_data_trs.loc[(sv_data_trs['RR'].isnull()) | (sv_data_trs['RR']=='NA') | (sv_data_trs['RR']=='.'),'RR']=0          
    sv_data_trs.loc[(sv_data_trs['RV'].isnull()) | (sv_data_trs['RV']=='NA') | (sv_data_trs['RV']=='.'),'RV']=0  
  
    sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['AP']!=0),'PR']=sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['AP']!=0),['RP','AP']].apply(lambda x: str(x['RP'])+','+str(x['AP']), axis=1)
    sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(',')) | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['AS']!=0),'SR']=sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(',')) | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['AS']!=0),['RS','AS']].apply(lambda x: str(x['RS'])+','+str(x['AS']), axis=1)             
    sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['DV']!=0),'PR']=sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['DV']!=0),['DR','DV']].apply(lambda x: str(x['DR'])+','+str(x['DV']), axis=1)
    sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(','))  | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['RV']!=0),'SR']=sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(','))  | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['RV']!=0),['RR','RV']].apply(lambda x: str(x['RR'])+','+str(x['RV']), axis=1)  
    sv_data_trs.loc[~sv_data_trs['PR'].str.contains(','),'PR']='0,0'
    sv_data_trs.loc[~sv_data_trs['SR'].str.contains(','),'SR']='0,0'         
         
    sv_data_trs['PR_ref']=sv_data_trs['PR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['PR_alt']=sv_data_trs['PR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_trs['SR_ref']=sv_data_trs['SR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['SR_alt']=sv_data_trs['SR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)        

    # replace na
    sv_data_trs=sv_data_trs.fillna(0)
    sv_data_trs.replace('', 0, inplace=True)    
    sv_data_trs.replace('.', 0, inplace=True) 
    
    # filter low quality sv
    sv_data_trs=sv_data_trs[(sv_data_trs['PR_alt']>=1) & (sv_data_trs['SR_alt']>=1) & ( sv_data_trs['BND_DEPTH'].astype(int) >=1) & ( sv_data_trs['MATE_BND_DEPTH'].astype(int) >=1) & (sv_data_trs['sv_bp_end_POS']!='NAN') & (sv_data_trs['sv_bp_end_POS']!='chrY')]

    # generate sv_info variables (length=bp_st-bp_end; range=ci_st-ci_end; ratio=read_alt/read_ref; diff=read_alt-read_ref)
    sv_data_trs['cipos_range']=sv_data_trs['ci_pos0'].astype(int) -sv_data_trs['ci_pos1'].astype(int)
    sv_data_trs['ciend_range']=sv_data_trs['sv_bp_end_POS'].astype(int) -sv_data_trs['sv_bp_end_END'].astype(int)
     
    sv_data_trs['PR_read_ratio']=abs(sv_data_trs['PR_alt'])/(abs(sv_data_trs['PR_ref'])+abs(sv_data_trs['PR_alt']))
    sv_data_trs['SR_read_ratio']=abs(sv_data_trs['SR_alt'])/(abs(sv_data_trs['SR_ref'])+abs(sv_data_trs['SR_alt']))
    sv_data_trs['read_ratio']=abs(sv_data_trs['PR_alt'])+abs(abs(sv_data_trs['SR_alt']))/(abs(sv_data_trs['SR_ref'])+abs(sv_data_trs['SR_alt'])+abs(sv_data_trs['PR_ref'])+abs(sv_data_trs['PR_alt']))
    sv_data_trs['bnd_depth_ratio']=abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))/(abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))+abs(sv_data_trs['BND_DEPTH'].astype(float)))
    sv_data_trs['PR_read_ratio_diff']=(abs(sv_data_trs['PR_alt'])-abs(sv_data_trs['PR_ref']))/(abs(sv_data_trs['PR_alt'])+abs(sv_data_trs['PR_ref']))
    sv_data_trs['SR_read_ratio_diff']=(abs(sv_data_trs['SR_alt'])-abs(sv_data_trs['SR_ref']))/(abs(sv_data_trs['SR_alt'])+ abs(sv_data_trs['SR_ref']))
    sv_data_trs['read_ratio_diff']=(abs(sv_data_trs['PR_alt'])+ abs(sv_data_trs['SR_alt'])-abs(sv_data_trs['PR_ref'])-abs(sv_data_trs['SR_ref']))/(abs(sv_data_trs['PR_alt'])+ abs(sv_data_trs['SR_alt'])+abs(sv_data_trs['PR_ref'])+abs(sv_data_trs['SR_ref']))
    sv_data_trs['bnd_depth_ratio_diff']=(abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))-abs(sv_data_trs['BND_DEPTH'].astype(float)))/(abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))+abs(sv_data_trs['BND_DEPTH'].astype(float)))
    print(sv_data_trs.shape)   
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs.shape[0])), k=sv_data_trs.shape[0]) 
    sv_data_trs=sv_data_trs.reset_index(drop=True)
    sv_data_trs=sv_data_trs.reindex(new_index)
    sv_data_trs=sv_data_trs.reset_index(drop=True)     
    # replace na
    sv_data_trs=sv_data_trs.fillna(0)    
    # replace infinite value
    sv_data_trs.replace([np.inf, -np.inf], 0, inplace=True) 
    
    # label benchmark sv 
    sv_data_trs.insert(sv_data_trs.shape[1], "label", list(repeat(0,sv_data_trs.shape[0])), True)  
    tg_idx0=np.where(((sv_data_trs['1000_g']=="YVID") | (sv_data_trs['gnomad_g']=="YVID") | (sv_data_trs['dgv_g']=="YVID")))
    ta_idx0=np.where(((sv_data_trs['gnomad_qc']=="YVID") | (sv_data_trs['control_g']=="YVID") | (sv_data_trs['centromere_qc']=="YVID")))    
    ts_idx0=np.where((((sv_data_trs['cytoatlas_s']=="YVID") | (sv_data_trs['cosmic_s']=="YVID") | (sv_data_trs['uwstl_s']=="YVID")) & (sv_data_trs['control_gall']!="YVID"))) 
    print(tg_idx0[0].shape)
    print(ta_idx0[0].shape)    
    print(ts_idx0[0].shape)       
    tg_idx=[i for i in tg_idx0[0]  if i not in ta_idx0[0] and i not in ts_idx0[0] ]
    ta_idx=[i for i in ta_idx0[0]  if i not in tg_idx0[0] and i not in ts_idx0[0] ] 
    ts_idx=[i for i in ts_idx0[0]  if i not in tg_idx0[0] and i not in ta_idx0[0] ]      
    print(len(tg_idx))
    print(len(ta_idx))
    print(len(ts_idx))    
    sv_data_trs.loc[np.r_[tg_idx],'label']=1   
    sv_data_trs.loc[np.r_[ta_idx],'label']=-1  
    sv_data_trs.loc[np.r_[ts_idx],'label']=2    
    print(sv_data_trs['label'].value_counts())    
    
    # make summary plot / table
    sv_data_trs['sv_start_bp0']=0
    sv_data_trs['sv_start_bp1']=0
    sv_data_trs['sv_start_bp0']=sv_data_trs['sv_start_bp'].astype(int)-sv_data_trs['ci_pos0'].astype(int)
    sv_data_trs['sv_start_bp1']=sv_data_trs['sv_start_bp'].astype(int)+sv_data_trs['ci_pos1'].astype(int)
    sv_data_trs['sv_chr12']=sv_data_trs['sv_chr1']+":"+sv_data_trs['sv_chr2']
    sv_data_trs.replace('NAN', 0, inplace=True) 
    sv_data_trs.replace('-', 0, inplace=True) 
    sv_data_trs.replace('', 0, inplace=True) 
    
    # drop columns wih rebundant variables
    sv_data_trs_2=sv_data_trs.copy()
    sv_data_trs_2=sv_data_trs_2.loc[:,['sv_type','sv_bp_end_cc1', 'sv_bp_end_cc_v1', 'sv_bp_end_cc_v2',  'sv_bp_end_cc_v3', 'sv_bp_end_cc_v4', 'sv_bp_end_cc_v5', 'sv_bp_end_cc_v6', 'sv_bp_end_cc_v7', 'sv_bp_end_cc_v8', 'sv_bp_end_cc_v9', 'sv_bp_end_cc_v10', 'sv_bp_end_cc_v11',  'sv_bp_end_cc_v12', 'sv_bp_end_cc_v13', 'sv_bp_end_cc_v14',   'sv_bp_end_cc_v15', 'sv_bp_end_cc_v16', 'sv_bp_end_cc_v17', 'sv_bp_end_cc_v18', 'sv_bp_end_cc_v19', 'sv_bp_end_cc_v20', 'sv_bp_end_cc_v21', 'sv_bp_end_cc_v22', 'sv_bp_end_cc_v23', 'sv_bp_end_cc_v24', 'sv_bp_st_cc1', 'sv_bp_st_cc_v1', 'sv_bp_st_cc_v2',  'sv_bp_st_cc_v3', 'sv_bp_st_cc_v4', 'sv_bp_st_cc_v5', 'sv_bp_st_cc_v6',  'sv_bp_st_cc_v7', 'sv_bp_st_cc_v8', 'sv_bp_st_cc_v9', 'sv_bp_st_cc_v10', 'sv_bp_st_cc_v11', 'sv_bp_st_cc_v12', 'sv_bp_st_cc_v13', 'sv_bp_st_cc_v14', 'sv_bp_st_cc_v15', 'sv_bp_st_cc_v16', 'sv_bp_st_cc_v17', 'sv_bp_st_cc_v18', 'sv_bp_st_cc_v19', 'sv_bp_st_cc_v20', 'sv_bp_st_cc_v21', 'sv_bp_st_cc_v22', 'sv_bp_st_cc_v23', 'sv_bp_st_cc_v24', 'cipos_range', 'ciend_range',  'PR_read_ratio', 'SR_read_ratio', 'read_ratio', 'bnd_depth_ratio',  'PR_read_ratio_diff', 'SR_read_ratio_diff', 'read_ratio_diff', 'bnd_depth_ratio_diff', 'label']]
    sv_data_trs_2.columns
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs_2.shape[0])), k=sv_data_trs_2.shape[0]) 
    sv_data_trs_2=sv_data_trs_2.reset_index(drop=True)
    sv_data_trs_2=sv_data_trs_2.reindex(new_index)
    sv_data_trs_2=sv_data_trs_2.reset_index(drop=True)
    return sv_data_trs_2,sv_data_trs

def nontrs_sv_data_transform(sv_data_nontrs, nontrs_sv_cutoff):
    """Transform nonTRS (DEL/DUP/INV) SV data for ML modeling.

    Args:
        sv_data_nontrs (pd.DataFrame): Raw nonTRS SV data with all features.
        nontrs_sv_cutoff (float): Min overlap ratio for database matching.
            SVs with overlap >= cutoff are labeled as matching the database.

    Returns:
        tuple: (sv_data_nontrs_2, sv_data_nontrs)
            - sv_data_nontrs_2: Reduced feature matrix (model features + label)
            - sv_data_nontrs: Full data matrix with all columns + label

    Note:
        Additional filters applied for nonTRS:
        - sv_length >= 1,000,000 bp (large SVs only for benchmarking)
        - Genotype must be 0/1 (heterozygous)
    """
    # Shuffle data randomly for unbiased training
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs.shape[0])), k=sv_data_nontrs.shape[0]) 
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)
    sv_data_nontrs=sv_data_nontrs.reindex(new_index)
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)   
    print(sv_data_nontrs.shape) 
    
    # transform database label
    sv_data_nontrs.loc[sv_data_nontrs['1000_gall']=='Not_in_database','1000_gall']=0
    sv_data_nontrs.loc[sv_data_nontrs['1000_g']=='Not_in_database','1000_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['gnomad_gall']=='Not_in_database','gnomad_gall']=0
    sv_data_nontrs.loc[sv_data_nontrs['gnomad_g']=='Not_in_database','gnomad_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['control_gall']=='Not_in_database','control_gall']=0
    sv_data_nontrs.loc[sv_data_nontrs['control_g']=='Not_in_database','control_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['cytoatlas_s']=='Not_in_database','cytoatlas_s']=0
    sv_data_nontrs.loc[sv_data_nontrs['cosmic_s']=='Not_in_database','cosmic_s']=0
    sv_data_nontrs.loc[sv_data_nontrs['gnomad_qc']=='Not_in_database','gnomad_qc']=0
    sv_data_nontrs.loc[sv_data_nontrs['dgv_g']=='Not_in_database','dgv_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['centromere_qc']=='Not_in_database','centromere_qc']=0
    sv_data_nontrs.loc[sv_data_nontrs['uwstl_s']=='Not_in_database','uwstl_s']=0
    
    # transform sv_info variables
    sv_data_nontrs['ci_pos0']=0
    sv_data_nontrs['ci_pos1']=0
    sv_data_nontrs['sv_bp_end_POS']=0
    sv_data_nontrs['sv_bp_end_END']=0
    sv_data_nontrs['PR_ref']=0
    sv_data_nontrs['PR_alt']=0
    sv_data_nontrs['SR_ref']=0
    sv_data_nontrs['SR_alt']=0
    sv_data_nontrs['ci_pos0']=sv_data_nontrs['CIPOS'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['ci_pos1']=sv_data_nontrs['CIPOS'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_nontrs['sv_bp_end_POS']=sv_data_nontrs['CIEND'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['sv_bp_end_END']=sv_data_nontrs['CIEND'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
#     print(sv_data_nontrs.shape) 
    
    # read statistics transformation
    sv_data_nontrs.loc[(sv_data_nontrs['RP'].isnull()) | (sv_data_nontrs['RP']=='NA') | (sv_data_nontrs['RP']=='.'),'RP']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['AP'].isnull()) | (sv_data_nontrs['AP']=='NA') | (sv_data_nontrs['AP']=='.'),'AP']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['RS'].isnull()) | (sv_data_nontrs['RS']=='NA') | (sv_data_nontrs['RS']=='.'),'RS']=0          
    sv_data_nontrs.loc[(sv_data_nontrs['AS'].isnull()) | (sv_data_nontrs['AS']=='NA') | (sv_data_nontrs['AS']=='.'),'AS']=0    
    sv_data_nontrs.loc[(sv_data_nontrs['DR'].isnull()) | (sv_data_nontrs['DR']=='NA') | (sv_data_nontrs['DR']=='.'),'DR']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['DV'].isnull()) | (sv_data_nontrs['DV']=='NA') | (sv_data_nontrs['DV']=='.'),'DV']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['RR'].isnull()) | (sv_data_nontrs['RR']=='NA') | (sv_data_nontrs['RR']=='.'),'RR']=0          
    sv_data_nontrs.loc[(sv_data_nontrs['RV'].isnull()) | (sv_data_nontrs['RV']=='NA') | (sv_data_nontrs['RV']=='.'),'RV']=0  
 
    sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['AP']!=0),'PR']=sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['AP']!=0),['RP','AP']].apply(lambda x: str(x['RP'])+','+str(x['AP']), axis=1)
    sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(',')) | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['AS']!=0),'SR']=sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(',')) | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['AS']!=0),['RS','AS']].apply(lambda x: str(x['RS'])+','+str(x['AS']), axis=1)             
    sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['DV']!=0),'PR']=sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['DV']!=0),['DR','DV']].apply(lambda x: str(x['DR'])+','+str(x['DV']), axis=1)
    sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(','))  | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['RV']!=0),'SR']=sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(','))  | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['RV']!=0),['RR','RV']].apply(lambda x: str(x['RR'])+','+str(x['RV']), axis=1)  
    sv_data_nontrs.loc[~sv_data_nontrs['PR'].str.contains(','),'PR']='0,0'
    sv_data_nontrs.loc[~sv_data_nontrs['SR'].str.contains(','),'SR']='0,0'       

         
    sv_data_nontrs['PR_ref']=sv_data_nontrs['PR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['PR_alt']=sv_data_nontrs['PR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_nontrs['SR_ref']=sv_data_nontrs['SR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['SR_alt']=sv_data_nontrs['SR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)        

    # replace na
    sv_data_nontrs=sv_data_nontrs.fillna(0)
    sv_data_nontrs.replace('', 0, inplace=True)    
    sv_data_nontrs.replace('.', 0, inplace=True)  
    
    # filter low quality sv
    sv_data_nontrs['sv_length']=sv_data_nontrs['sv_start_bp'].astype(int)-sv_data_nontrs['sv_end_bp'].astype(int)
    sv_data_nontrs=sv_data_nontrs[sv_data_nontrs['CIEND']!='NA'] 
  
    sv_data_nontrs=sv_data_nontrs[ (sv_data_nontrs['CN']!='2') | (abs(sv_data_nontrs['PR_alt'].astype(int))>1) | (abs(sv_data_nontrs['SR_alt'].astype(int))>0) | (abs(sv_data_nontrs['PR_ref'].astype(int))>0) |  (abs(sv_data_nontrs['SR_ref'].astype(int))>0) | (sv_data_nontrs['SUPP'].astype(int)>1)]
   
    # generate sv_info variables (length=bp_st-bp_end; range=ci_st-ci_end; ratio=read_alt/read_ref; diff=read_alt-read_ref)
    sv_data_nontrs['cipos_range']=sv_data_nontrs['ci_pos0'].astype(int) -sv_data_nontrs['ci_pos1'].astype(int)
    sv_data_nontrs['ciend_range']=sv_data_nontrs['sv_bp_end_POS'].astype(int) -sv_data_nontrs['sv_bp_end_END'].astype(int)
     
    sv_data_nontrs['PR_read_ratio']=abs(sv_data_nontrs['PR_alt'])/(abs(sv_data_nontrs['PR_ref'])+abs(sv_data_nontrs['PR_alt']))
    sv_data_nontrs['SR_read_ratio']=abs(sv_data_nontrs['SR_alt'])/(abs(sv_data_nontrs['SR_ref'])+abs(sv_data_nontrs['SR_alt']))
    sv_data_nontrs['read_ratio']=abs(sv_data_nontrs['PR_alt'])+abs(abs(sv_data_nontrs['SR_alt']))/(abs(sv_data_nontrs['SR_ref'])+abs(sv_data_nontrs['SR_alt'])+abs(sv_data_nontrs['PR_ref'])+abs(sv_data_nontrs['PR_alt']))
    sv_data_nontrs['PR_read_ratio_diff']=(abs(sv_data_nontrs['PR_alt'])-abs(sv_data_nontrs['PR_ref']))/(abs(sv_data_nontrs['PR_alt'])+abs(sv_data_nontrs['PR_ref']))
    sv_data_nontrs['SR_read_ratio_diff']=(abs(sv_data_nontrs['SR_alt'])-abs(sv_data_nontrs['SR_ref']))/(abs(sv_data_nontrs['SR_alt'])+ abs(sv_data_nontrs['SR_ref']))
    sv_data_nontrs['read_ratio_diff']=(abs(sv_data_nontrs['PR_alt'])+ abs(sv_data_nontrs['SR_alt'])-abs(sv_data_nontrs['PR_ref'])-abs(sv_data_nontrs['SR_ref']))/(abs(sv_data_nontrs['PR_alt'])+ abs(sv_data_nontrs['SR_alt'])+abs(sv_data_nontrs['PR_ref'])+abs(sv_data_nontrs['SR_ref']))
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs.shape[0])), k=sv_data_nontrs.shape[0]) 
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)
    sv_data_nontrs=sv_data_nontrs.reindex(new_index)
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)       
    # replace na
    sv_data_nontrs=sv_data_nontrs.fillna(0)
    # replace infinite value
    sv_data_nontrs.replace([np.inf, -np.inf], 0, inplace=True)
    print(sv_data_nontrs.shape)  
        
    # label benchmark sv 
    sv_data_nontrs=sv_data_nontrs[abs(sv_data_nontrs['sv_length'])>=1000000]   
    print(sv_data_nontrs.shape)     
    sv_data_nontrs=sv_data_nontrs[(sv_data_nontrs['GT']=='0/1')]
    print(sv_data_nontrs.shape)  
 
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs.shape[0])), k=sv_data_nontrs.shape[0]) 
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)
    sv_data_nontrs=sv_data_nontrs.reindex(new_index)
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)     
    sv_data_nontrs.insert(sv_data_nontrs.shape[1], "label", list(repeat(0,sv_data_nontrs.shape[0])), True)    
    tg_idx0=np.where(((sv_data_nontrs['1000_g'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['gnomad_g'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['dgv_g'].astype(float)>nontrs_sv_cutoff)))
    ta_idx0=np.where(((sv_data_nontrs['gnomad_qc'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['control_g'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['centromere_qc'].astype(float)>nontrs_sv_cutoff)))    
    ts_idx0=np.where((((sv_data_nontrs['cosmic_s'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['uwstl_s'].astype(float)>nontrs_sv_cutoff))  & (sv_data_nontrs['dgv_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['1000_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['gnomad_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['control_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['1000_gall'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['gnomad_gall'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['control_gall'].astype(float)<nontrs_sv_cutoff))) 
    print(tg_idx0[0].shape)
    print(ta_idx0[0].shape)    
    print(ts_idx0[0].shape)       
    tg_idx=[i for i in tg_idx0[0]  if i not in ta_idx0[0] and i not in ts_idx0[0] ]
    ta_idx=[i for i in ta_idx0[0]  if i not in tg_idx0[0] and i not in ts_idx0[0] ] 
    ts_idx=[i for i in ts_idx0[0]  if i not in tg_idx0[0] and i not in ta_idx0[0] ]     
    print(len(tg_idx))
    print(len(ta_idx))  
    print(len(ts_idx))      
    sv_data_nontrs.loc[np.r_[tg_idx],'label']=1    
    sv_data_nontrs.loc[np.r_[ta_idx],'label']=-1    
    sv_data_nontrs.loc[np.r_[ts_idx],'label']=2    


    print(sv_data_nontrs['label'].value_counts())      
    sv_data_nontrs=sv_data_nontrs.rename(columns={'QUAL':'svtyper_score','SUPP':'sv_caller_supp', 'CN':'sv_copy_number'})  
    sv_data_nontrs.replace('NAN', 0, inplace=True) 
    sv_data_nontrs.replace('NA', 0, inplace=True)
    sv_data_nontrs.replace('-', 0, inplace=True)  
    sv_data_nontrs.replace('', 0, inplace=True)  
    
    # drop columns wih rebundant variables
    sv_data_nontrs_2=sv_data_nontrs.copy()
    sv_data_nontrs_2=sv_data_nontrs_2.loc[:,['sv_type', 'svtyper_score', 'sv_bp_end_cc1', 'sv_bp_end_cc_v1','sv_bp_end_cc_v2', 'sv_bp_end_cc_v3', 'sv_bp_end_cc_v4','sv_bp_end_cc_v5', 'sv_bp_end_cc_v6', 'sv_bp_end_cc_v7', 'sv_bp_end_cc_v8', 'sv_bp_end_cc_v9', 'sv_bp_end_cc_v10', 'sv_bp_end_cc_v11', 'sv_bp_end_cc_v12', 'sv_bp_end_cc_v13', 'sv_bp_end_cc_v14', 'sv_bp_end_cc_v15', 'sv_bp_end_cc_v16','sv_bp_end_cc_v17', 'sv_bp_end_cc_v18', 'sv_bp_end_cc_v19', 'sv_bp_end_cc_v20', 'sv_bp_end_cc_v21', 'sv_bp_end_cc_v22', 'sv_bp_end_cc_v23', 'sv_bp_end_cc_v24', 'sv_bp_st_cc1', 'sv_bp_st_cc_v1', 'sv_bp_st_cc_v2', 'sv_bp_st_cc_v3', 'sv_bp_st_cc_v4', 'sv_bp_st_cc_v5', 'sv_bp_st_cc_v6', 'sv_bp_st_cc_v7', 'sv_bp_st_cc_v8', 'sv_bp_st_cc_v9', 'sv_bp_st_cc_v10', 'sv_bp_st_cc_v11', 'sv_bp_st_cc_v12', 'sv_bp_st_cc_v13', 'sv_bp_st_cc_v14', 'sv_bp_st_cc_v15', 'sv_bp_st_cc_v16', 'sv_bp_st_cc_v17', 'sv_bp_st_cc_v18', 'sv_bp_st_cc_v19', 'sv_bp_st_cc_v20', 'sv_bp_st_cc_v21', 'sv_bp_st_cc_v22', 'sv_bp_st_cc_v23', 'sv_bp_st_cc_v24', 'sv_copy_number', 'sv_caller_supp', 'cipos_range', 'ciend_range', 'PR_read_ratio', 'SR_read_ratio', 'read_ratio', 'PR_read_ratio_diff',  'SR_read_ratio_diff', 'read_ratio_diff','label']]
    sv_data_nontrs_2.columns
    
    # re-arrange the data 
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs_2.shape[0])), k=sv_data_nontrs_2.shape[0]) 
    sv_data_nontrs_2=sv_data_nontrs_2.reset_index(drop=True)
    sv_data_nontrs_2=sv_data_nontrs_2.reindex(new_index)
    sv_data_nontrs_2=sv_data_nontrs_2.reset_index(drop=True)    
    return sv_data_nontrs_2,sv_data_nontrs

def sv_data_summary_plot(sv_data):
    """Generate a stacked bar plot showing SV class distribution by SV type.

    Args:
        sv_data (pd.DataFrame): Transformed SV data with 'sv_type' and 'label' columns.

    Returns:
        matplotlib.figure.Figure: Stacked bar chart with counts per class per SV type.
            Classes: -1 (artifact), 1 (germline), 2 (somatic)
    """
    # Count labels by SV type
    sv_datat2=sv_data.copy()
    sv_datat2.insert(sv_datat2.shape[1], "count", list(repeat(1,sv_datat2.shape[0])), True)
    sv_datat2=sv_datat2.loc[:,['sv_type','label','count']]
    print(sv_datat2.groupby(['sv_type','label']).sum())
    SA=sv_datat2.loc[sv_datat2['label']==-1,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    UL=sv_datat2.loc[sv_datat2['label']==0,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    TG=sv_datat2.loc[sv_datat2['label']==1,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    TC=sv_datat2.loc[sv_datat2['label']==2,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    index=sv_datat2['sv_type'].unique()
    df1 = pd.DataFrame({'class -1': SA,'class 1': TG,'class 2': TC}, index=index)
    print(df1.iloc[0,:].sum())    
    ax = df1.plot.bar(stacked=True, figsize=(5, 7),width = 0.2,legend=None)
    print(df1)
    for idx, sv in enumerate(df1.index):
        df2= df1[df1.index==sv]
        d0=0
        for i in range(df2.shape[1]):
            data=df2.iloc[0,i]
            ax.text(x=idx-0.1 , y = d0+data/2*pow(-2,i), s=f"{data}" , fontdict=dict(fontsize=10))
            d0=d0+data
    ax.set_yscale('linear')
    ax.ticklabel_format(style='plain', axis='y')
    ax.set_ylabel("SV number")
    ax.set_ylim([0,pow(10,len(str(df1.iloc[0,:].sum())))*1.1])
    fig = ax.figure 
    return fig
