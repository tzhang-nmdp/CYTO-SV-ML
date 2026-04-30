"""
sv_size.py
==========
Filters SVs by size from a VCF file.

Usage:
  python sv_size.py <input.vcf> <size_threshold> <direction>

Arguments:
  input.vcf      - Input VCF file
  size_threshold - Size threshold in bp
  direction      - 'down': keep SVs >= threshold (filter small)
                   'up': keep SVs <= threshold (filter large)

Output:
  Filtered VCF lines printed to stdout
"""

import sys, os, re

in_vcf = open(sys.argv[1], 'r')
size = int(sys.argv[2])
cut_direction = str(sys.argv[3])

def vcf_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    if item[0].startswith('#'):
        line=' '.join(item)
    elif len(item)<10:
        line='NONE' 
    else:
        alt_list=item[4]        
        info_list=item[7]
        info=info_list.split(';')
        for inf in info:
            if "=" in inf:
                sv_dict[inf.split('=')[0]]=inf.split('=')[1]
        if 'SVTYPE' in sv_dict:
            sv_type=sv_dict['SVTYPE']
        else:
            sv_type=re.sub('<|>','',item[4])
        if sv_type=='BND' and ']' in alt_list:
            alt=alt_list.split(']')[1]
            sv_len=size
            sv_dict['END']=alt.split(':')[1]
            sv_chr2=alt.split(':')[0]            
        elif sv_type=='BND' and '[' in alt_list:
            alt=alt_list.split('[')[1]             
            sv_len=size
            sv_dict['END']=alt.split(':')[1]
            sv_chr2=alt.split(':')[0]
        #elif sv_type=='INS':
        #    sv_len=sv_dict['INSLEN']
        #    sv_dict['END']=item[1]
        #    sv_chr2=item[0]             
        else:
            sv_len=abs(int(item[1])-int(sv_dict['END']))+1 
            sv_chr2=item[0]       
        
        # check sv size for filtering
        if cut_direction=='down' and sv_len<size:
            line='NONE'
        elif cut_direction=='up' and sv_len>size:
            line='NONE' 
        else:
            line='\t'.join(item)
    return line  

for line in in_vcf:
    new_line=vcf_sv_sim(line)
    if not new_line.startswith('NONE'):
        print(new_line)