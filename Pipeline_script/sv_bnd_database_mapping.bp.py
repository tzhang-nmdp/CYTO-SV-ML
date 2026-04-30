"""
sv_bnd_database_mapping.bp.py
==============================
Annotates TRS (BND/translocation) SVs against a TRS SV database using
breakpoint distance matching WITH confidence interval (CI) support.

Unlike sv_bnd_database_mapping.py, this version uses databases that include
CI ranges (columns: chr1, bp1_start, bp1_end, chr2, bp2_start, bp2_end).
A match occurs when the query breakpoint falls within the CI range of a
database SV (expanded by bp_dis).

Usage:
  python sv_bnd_database_mapping.bp.py <database.bp.trs> <query_trs_file> <db_name_distance>

Arguments:
  database.bp.trs   - TRS database with CI info (6+ columns)
  query_trs_file    - Query TRS SV file with header
  db_name_distance  - Database name with distance suffix (e.g., '1000_g_1000')

Output:
  {query_trs_file}.{db_name_distance}
"""

import sys, os, re
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

ref_data=pd.read_csv(sys.argv[1],sep='\t',header=None)#,skiprows=1)
in_vcf=pd.read_csv(sys.argv[2],sep='\t',header=0,keep_default_na=False)
out_vcf=str(sys.argv[2])+"."+str(sys.argv[3])
bp_dis=int(str(sys.argv[3]).split('_')[2])

in_vcf['database']='Not_in_database'
for i in range(in_vcf.shape[0]):
    chr1=str(in_vcf.iloc[i,0])
    bp1=int(in_vcf.iloc[i,1])
    chr2=str(in_vcf.iloc[i,3])   
    bp2=int(in_vcf.iloc[i,2])
    bnd_dict=ref_data[(ref_data.iloc[:,0]==chr1) & (ref_data.iloc[:,3]==chr2) & (ref_data.iloc[:,1].astype(int)<=(bp1+bp_dis))  & (ref_data.iloc[:,2].astype(int)>=(bp1-bp_dis)) & (ref_data.iloc[:,4].astype(int)<=(bp2+bp_dis)) & (ref_data.iloc[:,5].astype(int)>=(bp2-bp_dis))] 

    if bnd_dict.empty or bnd_dict.shape[0]==0:
        continue
    else:       
        bnd_dict.iloc[:,1]=bnd_dict.iloc[:,1].astype(int)-bp1
        bnd_dict.iloc[:,2]=bp1-bnd_dict.iloc[:,2].astype(int) 
        bnd_dict.iloc[:,4]=bnd_dict.iloc[:,4].astype(int)-bp2
        bnd_dict.iloc[:,5]=bp2-bnd_dict.iloc[:,5].astype(int) 
        bnd_dict['info']=bnd_dict.apply(lambda x: max(0,x[1],x[2], x[4], x[5]),axis=1)  
        #print(bnd_dict)        
        in_vcf.iloc[i,in_vcf.shape[1]-1]=min(bnd_dict['info']) 
in_vcf.to_csv(out_vcf,sep='\t',index=False,header=False)
