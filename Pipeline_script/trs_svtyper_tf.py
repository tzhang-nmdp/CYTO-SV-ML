"""
trs_svtyper_tf.py
==================
Prepares TRS (BND) VCF records for SVTyper compatibility.

SVTyper requires MATEID to be properly formatted. This script swaps the
last digit of MATEID (0<->1) to create proper mate pair references.

Usage:
  python trs_svtyper_tf.py <input_trs.vcf>

Output:
  {input_trs.vcf}.tmp
"""

import sys, os, re, collections

in_vcf = open(sys.argv[1], 'r')
out_vcf = open(str(sys.argv[1]) + '.tmp', 'w')

def trs_svtyper_tf(line):
    info_dict={}
    item=line.strip().split('\t')
    
    # info column extraction
    info=item[7]
    info_list=info.split(';') 
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
        else:
            info_dict[inf]=inf
    
    # correct mate_id
    mate_id=info_dict['MATEID']   
    mate_id_key=mate_id.split(':')[-1]
    mate_id_cr=mate_id.split(':')[0:-1]  
 #   print(mate_id)
    if int(mate_id_key)==0:
      item[2]=':'.join(str(m) for m in mate_id_cr)+':1'
    elif int(mate_id_key)==1:
      item[2]=':'.join(str(m) for m in mate_id_cr)+':0' 
    else:
      print ("error!!!")
    line='\t'.join(str(l) for l in item)+'\n'
    return line
  
for line in in_vcf:
    if line.startswith('#'):
        out_vcf.write(line)             
    else:
        out_vcf.write(trs_svtyper_tf(line))   
        n=1
