"""
sv_bnd_database_mapping.py
==========================
Annotates TRS (BND/translocation) SVs against a TRS SV database using
breakpoint distance matching.

For each query TRS SV, finds database SVs where both breakpoints are
within the specified distance threshold. Reports the minimum distance.

Usage:
  python sv_bnd_database_mapping.py <database.trs> <query_trs_file> <db_name_distance>

Arguments:
  database.trs      - TRS database file (chr1, bp1, bp2, chr2, type)
  query_trs_file    - Query TRS SV file with header
  db_name_distance  - Database name with distance suffix (e.g., '1000_g_1000')
                      Distance is parsed from the last '_' token

Output:
  {query_trs_file}.{db_name_distance}
"""

import sys, os, re
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# read the sv input data and database
ref_data=pd.read_csv(sys.argv[1],sep='\t',header=None)
in_vcf=pd.read_csv(sys.argv[2],sep='\t',header=0,keep_default_na=False)
out_vcf=str(sys.argv[2])+"."+str(sys.argv[3])
bp_dis=int(str(sys.argv[3]).split('_')[2])

in_vcf['database']='Not_in_database'
for i in range(in_vcf.shape[0]):
    
    # define sv chromosome coordinates
    chr1=str(in_vcf.iloc[i,0])
    bp1=int(in_vcf.iloc[i,1])
    chr2=str(in_vcf.iloc[i,3])   
    bp2=int(in_vcf.iloc[i,2])
    
    # extract mapping records 
    bnd_dict=ref_data[(ref_data.iloc[:,0]==chr1) & (ref_data.iloc[:,3]==chr2) & (abs(ref_data.iloc[:,1].astype(int)-bp1)<=bp_dis) & (abs(ref_data.iloc[:,2].astype(int)-bp2)<=bp_dis)] 
    if bnd_dict.empty or bnd_dict.shape[0]==0:
        continue
    else:
        bnd_dict.iloc[:,1]=abs(bnd_dict.iloc[:,1].astype(int)-bp1)
        bnd_dict.iloc[:,2]=abs(bnd_dict.iloc[:,2].astype(int)-bp2) 
        bnd_dict['info']=bnd_dict.apply(lambda x: max(x[1:3]),axis=1)
        in_vcf.iloc[i,in_vcf.shape[1]-1]=min(bnd_dict['info']) 

# export to csv
in_vcf.to_csv(out_vcf,sep='\t',index=False,header=False)
