"""
sv_db_tf.py
============
Merges SV database annotation results into the cumulative annotation file.

For each database, joins the annotation label (overlap ratio or breakpoint
distance) to the running _anno file by SV ID.

Usage:
  python sv_db_tf.py <sv_file_prefix> <database_name> [force_mode]

Arguments:
  sv_file_prefix - Path prefix (reads {prefix}.{database_name} and {prefix}_anno)
  database_name  - Name of the SV database (e.g., '1000_g', 'cosmic_s_1000')
  force_mode     - 'f' to force overwrite existing column

Output:
  Updates {sv_file_prefix}_anno in place
"""

import sys, os, re
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

# read the sv input and database vcf
SV_database_name=str(sys.argv[2])
out_vcf=str(sys.argv[1])+'_anno' 
in_vcf=pd.read_csv(out_vcf,sep='\t',header=0,keep_default_na=False)  
db_vcf=pd.read_csv(str(sys.argv[1])+'.'+str(sys.argv[2]),sep='\t',header=None,names=range(10))    
try:
    id_mode=str(sys.argv[3])
except:
    id_mode='n'

# simplify the sv database vcf
db_vcf=db_vcf.rename(columns={db_vcf.columns[6]: SV_database_name})

if (db_vcf.iloc[:,4]=='TRA').any():
    db_vcf.loc[db_vcf[SV_database_name]!='Not_in_database',SV_database_name]=db_vcf.iloc[:,6]
else:
    db_vcf.loc[db_vcf[SV_database_name]!='Not_in_database',SV_database_name]=db_vcf.iloc[:,8]
db_vcf_sim=db_vcf.iloc[:, [5,6]]

db_vcf_sim.columns=['sv_id',SV_database_name]


# merge the sv database vcf
if SV_database_name not in in_vcf.columns:
    in_db_vcf=pd.merge(in_vcf, db_vcf_sim, on=['sv_id'], how='inner')

    # label NA database
#    if  in_db_vcf[SV_database_name].isnull().any():
#    if  db_vcf_sim['db_label'].str.contains('NAN').any():        
    in_db_vcf.loc[in_db_vcf[SV_database_name].isnull(),SV_database_name]='Not_in_database'   
elif id_mode=='f':
    in_vcf=in_vcf.drop(SV_database_name, axis=1) 
    in_db_vcf=pd.merge(in_vcf, db_vcf_sim, on=['sv_id'], how='inner')
    in_db_vcf=in_db_vcf.rename(columns={in_db_vcf.columns[in_db_vcf.shape[1]-1]:SV_database_name})

    # label NA database
#    if  in_db_vcf[SV_database_name].isnull().any():
#    if  db_vcf_sim['db_label'].str.contains('NAN').any():        
    in_db_vcf.loc[in_db_vcf[SV_database_name].isnull(),SV_database_name]='Not_in_database'     
else:
    in_db_vcf=in_vcf
# export to csv
in_db_vcf.to_csv(out_vcf,sep='\t',index=False,header=True)
