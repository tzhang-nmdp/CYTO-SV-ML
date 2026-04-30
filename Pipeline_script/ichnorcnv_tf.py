"""
ichnorcnv_tf.py
================
Transforms ichorCNV segment output (.segs.txt) into VCF format.

ichorCNV produces copy number segments (GAIN/AMP/HETD/HOMD) which are
converted to standard SV VCF records (DUP/DEL) for downstream processing.

Usage:
  python ichnorcnv_tf.py <input_segs.txt> <output.vcf>

Input format (tab-separated, with header):
  ID  chr  start  end  ...  copy_number  log_RD  call_type(GAIN/AMP/HETD/HOMD/NEUT)

Output:
  VCF 4.1 format with SVTYPE=DEL or DUP, GT, CN fields
"""

import sys
import os
import re

in_vcf = open(sys.argv[1], 'r')
out_vcf = open(sys.argv[2], 'w')

# Extract sample ID from input filename
sm_id = (str(sys.argv[1]).split('/')[-1]).split('.')[0]


def ichnorcnv_tf(line):
    """Convert a single ichorCNV segment line to VCF format.

    Mapping:
      GAIN/AMP -> DUP (duplication), positive SVLEN
      HETD     -> DEL (deletion), genotype 0/1, negative SVLEN
      HOMD     -> DEL (deletion), genotype 1/1, negative SVLEN
      NEUT     -> skipped (filtered in main loop)
    """
    item = line.strip().split('\t')

    # Header line: convert to VCF #CHROM header
    if item[0].startswith('ID'):
        line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sm_id + '\n'
        return line

    # Data line: parse segment coordinates and call type
    genotype = '0/1'
    sv_end = item[3]
    sv_type = 'NA'

    if item[11] == 'GAIN' or item[11] == 'AMP':
        sv_type = 'DUP'
        sv_len = int(sv_end) - int(item[2])
    elif item[11] == 'HETD':
        sv_type = 'DEL'
        sv_len = int(item[2]) - int(sv_end)
    elif item[11] == 'HOMD':
        sv_type = 'DEL'
        sv_len = int(item[2]) - int(sv_end)
        genotype = '1/1'  # Homozygous deletion

    # Build VCF fields
    sv_id = ':'.join([item[0], item[1], item[2]])
    sv_format = 'GT:CN'
    sv_cnv = item[10]       # Copy number
    log_RD = item[9]        # Log read depth ratio

    sv_info = (
        'END=' + str(sv_end) +
        ';SVTYPE=' + sv_type +
        ';SVLEN=' + str(sv_len) +
        ';IMPRECISE' +
        ';natorRD=' + str(log_RD)
    )

    line = '\t'.join(str(w) for w in [
        item[1], item[2], sv_id, '.', sv_type, '.', '.',
        sv_info, sv_format, genotype + ':' + str(sv_cnv)
    ]) + '\n'
    return line


# Write VCF header
out_vcf.write(
    '##fileformat=VCFv4.1\n'
    '##fileDate=20221013\n'
    '##reference=1000GenomesPhase3_decoy-GRCh37\n'
    '##source=CNVnator\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n'
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
    '##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">\n'
    '##INFO=<ID=natorP1,Number=1,Type=Float,Description="e-val by t-test">\n'
    '##INFO=<ID=natorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">\n'
    '##INFO=<ID=natorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">\n'
    '##INFO=<ID=natorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">\n'
    '##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">\n'
    '##INFO=<ID=natorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">\n'
    '##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">\n'
    '##ALT=<ID=DEL,Description="Deletion">\n'
    '##ALT=<ID=DUP,Description="Duplication">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n'
    '##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-ends that support the event">\n'
)

# Process each segment line (skip NEUT = neutral copy number)
for line in in_vcf:
    item = line.strip().split('\t')
    if item[11] != 'NEUT':
        out_vcf.write(ichnorcnv_tf(line))
