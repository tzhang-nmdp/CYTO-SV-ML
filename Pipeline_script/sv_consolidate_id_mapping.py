"""
sv_consolidate_id_mapping.py
=============================
Creates a mapping between consolidated SV IDs (from SURVIVOR merge) and
the original per-caller SV IDs.

Parses the FORMAT:ID field from the merged VCF to extract individual caller
IDs and maps them to the consolidated chr:pos:end:chr2:type format.

Usage:
  python sv_consolidate_id_mapping.py <merged.vcf>

Output:
  {merged.vcf}.sv_id_mapping (tab-separated: consolidated_id, caller_id)
"""

import sys, os, re

in_vcf = open(sys.argv[1], 'r')
out_vcf = open(str(sys.argv[1]) + '.sv_id_mapping', 'w')

def sv_id_tf(line):
    item=line.strip().split('\t')
    sv_chr2=item[0]
    sv_end=item[1]
    sv_type=re.sub('<|>','',item[4])
    
    #info column extraction
    info=item[7]
    info_list=info.split(';') 
    info_dict={}    
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
        else:
            info_dict[inf]=inf  
            
    # correct sv info        
    for inf in info_list:        
        if 'CHR2' in inf:
            sv_chr2=inf.split('=')[1]
        if 'END' in inf and 'CI' not in inf:
            sv_end=inf.split('=')[1]   
        if 'SVTYPE' in inf:
            sv_type=inf.split('=')[1]      
    if sv_type=='TRA':
        sv_type='BND'
    if info_dict['SVTYPE']=='TRA':
        info_dict['SVTYPE']='BND'        
    
    # consolidate sv_id column extraction
    if 'chr' in item[2]:
        consolidate_sv_id=item[2]
    else:
        consolidate_sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+ sv_type +':'+item[2]    
 
    # inv sv_id column extraction
    sv_id_list=[]    
    id_index_list=item[8].split(':')  # "format list"  
    id_index=id_index_list.index('ID') # id forward index in "format list"
    id_len=len(id_index_list)-id_index-1 # id backward index in "format list"
    for i in range(9,len(item)):
        geno_list=item[i].split(':')           
        inv_sv_id=':'.join(str(w) for w in geno_list[id_index:(len(geno_list)-id_len)])  # id = forward - backward in "format list"  
        if "CNVnator" not in inv_sv_id and "_" in inv_sv_id:
            inv_sv_id=re.sub('_',':',inv_sv_id)  # correct id (like Manta id)
        if "CNVnator" in inv_sv_id and ":" not in inv_sv_id:
            inv_sv_id=re.sub('_',':',inv_sv_id)  # correct id (like CNVnator id)    
            inv_sv_id_list=inv_sv_id.split(':')
            inv_sv_id=':'.join(str(w) for w in inv_sv_id_list[:-3])+':'+inv_sv_id_list[-3]+':'+inv_sv_id_list[-2]+':'+inv_sv_id_list[-1]
        if (inv_sv_id=='NaN' or inv_sv_id=='.') and geno_list[-1]!='NAN':
            inv_sv_co=re.sub('_|-',':',geno_list[-1].split(',')[0]).split(':')
            inv_sv_id=':'.join(str(w) for w in [inv_sv_co[i] for i in [0,1,3,2]])+':'+sv_type
        if inv_sv_id!='NaN' and inv_sv_id!='.' and 'chr' not in inv_sv_id:                   
            inv_sv_co=re.sub('_|-',':',geno_list[-1].split(',')[0]).split(':')
            inv_sv_id=':'.join(str(w) for w in [inv_sv_co[i] for i in [0,1,3,2]])+':'+sv_type+':'+inv_sv_id
        # if inv_sv_id!='NaN':         
        #     print(inv_sv_id)
        if inv_sv_id!='NaN' and inv_sv_id!='.' and 'chr' in inv_sv_id:                   
            sv_id_list.append(consolidate_sv_id+'\t'+inv_sv_id)
    # print out mapping consolidate sv_id and inv sv_id
    line='\n'.join(str(sv_id) for sv_id in sv_id_list)+'\n'        
    return line
  
for line in in_vcf:
    item=line.strip().split('\t')
    if line.startswith('#CHROM'):
        out_vcf.write("sv_id\tID\n")           
    elif line.startswith('#'):
        continue       
    else:
        out_vcf.write(sv_id_tf(line))   
