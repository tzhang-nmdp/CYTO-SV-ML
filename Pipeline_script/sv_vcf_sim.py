"""
sv_vcf_sim.py
==============
Splits a consolidated SV VCF into simplified TRS and nonTRS bed-like files.

TRS output format:
  sv_chr  sv_start_bp  sv_end_bp  sv_chr2  sv_type  sv_id

NonTRS output format:
  sv_chr  sv_start_bp  sv_end_bp  sv_type  sv_id

Usage:
  python sv_vcf_sim.py <input.vcf> [id_mode]

Output:
  {input.vcf}.trs     - TRS (BND/TRA) SVs
  {input.vcf}.nontrs  - nonTRS (DEL/DUP/INV) SVs
"""

import sys,os,re
in_vcf=open(sys.argv[1],'r')
trs_out_vcf=open(str(sys.argv[1])+'.trs','w')
non_trs_out_vcf=open(str(sys.argv[1])+'.nontrs','w')
try:
    id_mode=str(sys.argv[2])
except:
    id_mode='a'
    
def trs_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';')  
    for inf in info:
        if inf.startswith('SVTYPE'):
          sv_type=inf.split('=')[1]            
        if inf.startswith('END'):
          sv_bp_end=inf.split('=')[1]
        if inf.startswith('CHR2'):
          sv_chr2=inf.split('=')[1] 
    if 'chr' not in item[2] and id_mode!='s':
        line=str(item[0])+'\t'+str(item[1])+'\t'+str(sv_bp_end)+'\t'+str(sv_chr2)+'\t'+sv_type+'\t'+str(item[0])+':'+str(item[1])+':'+str(sv_bp_end)+':'+str(sv_chr2)+':BND:'+str(item[2])+'\n'
    else:
            line=str(item[0])+'\t'+str(item[1])+'\t'+str(sv_bp_end)+'\t'+str(sv_chr2)+'\t'+sv_type+'\t'+item[2]+'\n'
    return line  
  
def nontrs_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';')  
    for inf in info:
        if inf.startswith('SVTYPE'):
          sv_type=inf.split('=')[1]          
        if inf.startswith('END'):
          sv_bp_end=inf.split('=')[1]
    if 'chr' not in item[2] and id_mode!='s':       
        line=str(item[0])+'\t'+str(item[1])+'\t'+str(sv_bp_end)+'\t'+sv_type+'\t'+str(item[0])+':'+str(item[1])+':'+str(sv_bp_end)+':'+str(item[0])+':'+str(sv_type)+':'+str(item[2])+'\n'
    else:
            line=str(item[0])+'\t'+str(item[1])+'\t'+str(sv_bp_end)+'\t'+sv_type+'\t'+item[2]+'\n'
    return line  
  
for line in in_vcf:
    item=line.strip().split('\t')    
    if line.startswith('##'):
        continue
    elif item[0].startswith('#CHROM'):      
        trs_out_vcf.write('sv_chr\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type\tsv_id\n')
        non_trs_out_vcf.write('sv_chr\tsv_start_bp\tsv_end_bp\tsv_type\tsv_id\n')
    elif 'SVTYPE=TRA' in item[7] or 'SVTYPE=BND' in item[7] or 'MantaBND' in item[7]:
        trs_out_vcf.write(trs_sv_sim(line))        
    else:
        non_trs_out_vcf.write(nontrs_sv_sim(line))       
