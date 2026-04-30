"""
sv_info_tf_sim.py
==================
Extracts key VCF INFO and FORMAT fields into a simplified tab-separated file.

Default fields extracted:
  END, CHR2, CIPOS, CIEND, SVTYPE, BND_DEPTH, MATE_BND_DEPTH,
  GT, CN, PR, SR, DR, DV, RR, RV, RP, AP, RS, AS, ASC, PE

Usage:
  python sv_info_tf_sim.py <input.vcf> [keep_fields]

Arguments:
  input.vcf    - Input VCF file
  keep_fields  - Optional '|'-delimited list of fields to extract

Output:
  {input.vcf}.sv_info.sim
"""
import sys, os, re, collections
import numpy as np

in_vcf=open(sys.argv[1],'r')
out_vcf=open(str(sys.argv[1])+'.sv_info.sim','w')
try:
    keep_list=str(sys.argv[2]).split('|')
except:
#    keep_list=['END','CHR2','CIPOS','CIEND','SVTYPE','BND_DEPTH','MATE_BND_DEPTH','GT','CN','RP','AP','RS','AS','ASC']    
    keep_list=['END','CHR2','CIPOS','CIEND','SVTYPE','BND_DEPTH','MATE_BND_DEPTH','GT','CN','PR','SR','DR','DV','RR','RV','RP','AP','RS','AS','ASC','PE']
    
def sv_info_tf(line):
    info_dict={}
    info_key_list=[]
    item=line.strip().split('\t')
    sv_chr2=item[0]
    sv_end=item[1]
    sv_type=re.sub('<|>','',item[4])
    
    #info column extraction
    info=item[7]
    info_list=info.split(';') 
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
            info_key_list.append(inf.split('=')[0])
                
    # correct sv info        
    if 'SVTYPE' not in info_list:
        if 'DEL|DUP|INV|INS|TRA|BND' in item[2]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[2])[0]
        elif 'DEL|DUP|INV|INS|TRA|BND' in item[4]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[4])[0]
            
    alt_info=re.sub('\[|\]',':',item[4])     
    if 'CHR2' not in info_key_list:
        if info_dict['SVTYPE']=='BND':
            for ai in alt_info.split(':'):
                if re.findall('chr',ai):
                    info_dict['CHR2']=ai
        else:
            info_dict['CHR2']=item[0]   
            
    if 'END' not in info_key_list:
        if info_dict['SVTYPE']=='BND':
            for ai in alt_info.split(':'):
                if ai.isdigit():
                    info_dict['END']=ai
        else:
            info_dict['END']=item[1]
            
    if len(idx)!=0:
        keep_var.update(zip(idx_k,[item[x] for x in idx]))
        
    # geno column extraction
    geno_index=item[8].split(':')     
    geno_list=item[9].split(':') 
    if 'ID' in geno_index and len(geno_index)<len(geno_list):
        id_index=geno_index.index('ID') 
        id_len=len(geno_index)-id_index-1    
        inv_sv_id=':'.join(str(w) for w in geno_list[id_index:(len(geno_list)-id_len)])  
        for i in range(0,len(geno_index)):
            if i<id_index:
                info_dict[geno_index[i]]=geno_list[i]  
            elif i==id_index:
                info_dict[geno_index[i]]=inv_sv_id
            else:
                info_dict[geno_index[i]]=geno_list[i+id_len]    
    else:
        for i in range(0,len(geno_list)):
            info_dict[geno_index[i]]=geno_list[i]      
            
    # print out dict into line 
    key_list=[]
    value_list=[]
    key_all=[]
    for key,value in info_dict.items():
        key_all.append(str(key))
    for kp in keep_list:
        if kp not in key_all:
            info_dict[kp]='NA'
    info_dict = collections.OrderedDict(sorted(info_dict.items()))    
    for key,value in info_dict.items():
        if key in keep_list:
            key_list.append(str(key))
            value_list.append(str(value))  
    
    # correct info
    for inf in info_list:        
        if 'CHR2' in inf:
            sv_chr2=inf.split('=')[1]
        if 'END' in inf and 'CI' not in inf:
            sv_end=inf.split('=')[1]   
        if 'SVTYPE' in inf:
            sv_type=inf.split('=')[1]       
            
    # fix sv_id column extraction
    if 'chr' in item[2] or re.findall('DEL|DUP|INV|INS|TRA|BND',item[2]) or re.findall('del|dup|inv|ins|tra|bnd',item[2]) or ':' in item[2]:
        sv_id=item[2]
    else:
        sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+sv_type
    item[2]=sv_id
    
    if n==0:
        line='\t'.join(str(ch) for ch in col_head)+'\t'+'\t'.join(str(key) for key in key_list)+'\n'+'\t'.join(str(itm) for itm in item[0:6])+'\t'+'\t'.join(str(value) for value in value_list)+'\n'
    else:
        line='\t'.join(str(itm) for itm in item[0:6])+'\t'+'\t'.join(str(value) for value in value_list)+'\n'        
    return line

n=0  
for line in in_vcf:
    item=line.strip().split('\t')
    if line.startswith('#CHROM'):
        if '\t' not in line:
            item=line.strip().split(' ')        
        col_head=item[0:6]
        col_head[0]=re.sub("#", "", col_head[0])
        idx= [x for x, v in enumerate (item) if v in keep_list]
        idx_k=[item[x] for x in idx]           
        continue    
    elif line.startswith('#'):
        continue           
    else:
        out_vcf.write(sv_info_tf(line))   
        n=1
