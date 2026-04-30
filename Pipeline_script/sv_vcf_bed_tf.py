"""
sv_vcf_bed_tf.py
=================
Converts a consolidated SV VCF to BED format for breakpoint extraction.

Output format (tab-separated):
  sv_chr  sv_start_bp  sv_end_bp  sv_chr2  sv_type  sv_id

Usage:
  python sv_vcf_bed_tf.py <input.vcf>

Output:
  {input.vcf}.bed
"""

import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_bed=open(str(sys.argv[1])+'.bed','w')

def vcf_bed_tf(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';') 
    sv_chr2=item[0]
    sv_bp_end=item[1]
    if 'DEL|DUP|INV|INS|TRA|BND' in item[2]:
        sv_type=re.findall('DEL|DUP|INV|INS|TRA|BND',item[2])[0]
    for inf in info:
        if inf.startswith('END='):
            sv_bp_end=inf.split('=')[1]
        if inf.startswith('CHR2='):
            sv_chr2=inf.split('=')[1]    
        if inf.startswith('SVTYPE='):
            sv_type=inf.split('=')[1] 
    if 'chr' not in item[2]:       
        item[2]=str(item[0])+':'+str(item[1])+':'+str(sv_bp_end)+':'+str(sv_chr2)+':'+str(sv_type)+':'+str(item[2])  
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(sv_bp_end)+'\t'+str(sv_chr2)+'\t'+str(sv_type)+'\t'+ str(item[2])+'\n'
    return line  
  
for line in in_vcf:
    item=line.strip().split('\t')    
    if line.startswith('#'):
        continue
    else:
        out_bed.write(vcf_bed_tf(line))        
