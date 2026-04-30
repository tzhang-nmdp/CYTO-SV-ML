"""
sv_id_tf.py
===========
Standardizes SV IDs in VCF files to a consistent format:
  chr:pos:end:chr2:svtype[:original_id]

This ensures all callers use the same ID convention for downstream merging.

Usage:
  python sv_id_tf.py <input.vcf> [id_mode]

Modes:
  'a' (default): auto-detect; keep existing chr-based IDs, add format to others
  'c': force chr:pos:end:chr2:type:ID format for non-chr IDs
  'f'/'sf': always force full format

Output:
  {input.vcf}.re_id
"""

import sys, os, re

in_vcf = open(sys.argv[1], 'r')
out_vcf = open(str(sys.argv[1]) + '.re_id', 'w')
try:
    id_mode = str(sys.argv[2])
except:
    id_mode = 'a'

def sv_id_tf(line):
    item=line.strip().split('\t')
    sv_chr2=item[0]
    sv_end=item[1]
    sv_type=re.sub('<|>','',item[4])
    
    #collect info from column extraction
    info=item[7]
    info_list=info.split(';') 
    info_dict={}    
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
        else:
            info_dict[inf]=inf+"&no&id"       

    # correct sv info        
    if 'SVTYPE=' not in info_list:
        if 'DEL|DUP|INV|INS|TRA|BND' in item[2]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[2])[0]
        elif 'DEL|DUP|INV|INS|TRA|BND' in item[4]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[4])[0]
            
    alt_info=re.sub('\[|\]',':',item[4])     
    if 'CHR2=' not in info_list:
        if info_dict['SVTYPE']=='BND' or info_dict['SVTYPE']=='TRA':
            for ai in alt_info.split(':'):
                if re.findall('chr',ai):
                    info_dict['CHR2']=ai
        else:
            info_dict['CHR2']=item[0]   
            
    if 'END=' not in info:
        if info_dict['SVTYPE']=='BND' or info_dict['SVTYPE']=='TRA':
            for ai in alt_info.split(':'):
                if ai.isdigit():
                    info_dict['END']=ai
        else:
            info_dict['END']=item[1]
            
    sv_end=info_dict['END']   
    sv_chr2=info_dict['CHR2']
    sv_type=info_dict['SVTYPE']  
 #   print(str(item[1])+":"+str(sv_end))        
    # fix sv_id column extraction
    if id_mode=='f' or id_mode=='sf':
        sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+sv_type+":"+item[2]
    elif id_mode=='c' and 'chr' not in item[2]:
        sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+sv_type+":"+item[2]        
    elif 'chr' in item[2] or re.findall('DEL|DUP|INV|INS|TRA|BND',item[2]) or re.findall('del|dup|inv|ins|tra|bnd',item[2]) or ':' in item[2]:
        sv_id=item[2]
    else:
        sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+sv_type+":"+item[2]
    item[2]=sv_id
    
    # re-compose info columns
    info_new=""
    m=0
    for key,value in info_dict.items():
        if "&no&id" in key and m==0 and key!="CSQ":
            info_new=value
            m=1
        elif "&no&id" in key and m==1 and key!="CSQ":
            info_new=info_new+";"+value
            m=1
        elif "&no&id" not in key and m==0 and key!="CSQ":
            info_new=str(key)+"="+str(value) 
            m=1
        elif "&no&id" not in key and m==1 and key!="CSQ":
            info_new=info_new+";"+str(key)+"="+str(value) 
            m=1   
#    print(info_new)
    if ';CSQ=' in item[7] and id_mode!='sf': 
        item[7]=info_new+';CSQ='+info_dict['CSQ']
    else:
        item[7]=info_new
    
    # print out line
    line='\t'.join(str(w) for w in item)+'\n'        
    return line
  
for line in in_vcf:
    item=line.strip().split('\t')
    n=0
    if line.startswith('#'):
        out_vcf.write(line)        
    else:
        out_vcf.write(sv_id_tf(line))   
 
