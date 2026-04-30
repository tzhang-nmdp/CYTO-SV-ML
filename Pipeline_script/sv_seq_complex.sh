#!/bin/bash
###############################################################################
# sv_seq_complex.sh
#
# Computes breakpoint sequence complexity features for each SV.
# Extracts 300bp flanking sequences around start and end breakpoints,
# then runs SeqComplex (24 metrics) and Komplexity (1 metric) on each.
#
# Steps:
#   1. Convert consolidated VCF to BED format
#   2. Extract 300bp windows around start (bpst) and end (bpend) breakpoints
#   3. Extract FASTA sequences using bedtools getfasta
#   4. Remove low-complexity (poly-N) sequences
#   5. Run SeqComplex (profileComplexSeq.pl) for 24 complexity features
#   6. Run Komplexity (kz) for 1 complexity feature
#   7. Consolidate start + end breakpoint complexity into one file
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - sample: sample ID
#   $4 - py27_dir: path to Python 2.7 bin directory (for bedtools)
#   $5 - size_k: SV size threshold in kb (e.g., 10)
#
# Output:
#   {sample}.{size_k}k.sv.all.bed.bpst_bpend.kz.index_complex
#   (53 columns: 1 ID + 26 end-bp features + 26 start-bp features)
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
py27_dir=$4
size_k=$5

# ---- Step 1: Convert VCF to BED format ----
echo "# extract SV coordinate information" && date
python ${cyto_sv_ml_dir}/Pipeline_script/sv_vcf_bed_tf.py \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all

# ---- Step 2: Extract 300bp windows around breakpoints ----
# Start breakpoint (bpst): 150bp upstream + 150bp downstream
# Clamp to chromosome boundaries using hg38_chromosome_size.txt
awk 'FNR==NR{a[$1];b[$1]=$2;next}{
    c=b[$1]-150;
    if (($2>=150)&&($2<=c)) {$2=$2-150; $3=$2+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6}
    else if ($2>c) {$2=c-150;$3=c+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6}
    else if ($2<150){$2=1;$3=300; print $1"\t"$2"\t"$3"\t"$5"\t"$6}
}' ${cyto_sv_ml_dir}/reference/hg38_chromosome_size.txt \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed \
    | sed 's% %\t%g' \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst

# End breakpoint (bpend): swap chr2 as primary, extract 300bp window
awk 'FNR==NR{a[$1];b[$1]=$2;next}{
    $1=$4; $4=$1;
    c=b[$1]-150;
    if ($3>=c) {$3=c+150;$2=c-150; print $1"\t"$2"\t"$3"\t"$5"\t"$6}
    else {$2=$3-150; $3=$3+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6}
}' ${cyto_sv_ml_dir}/reference/hg38_chromosome_size.txt \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed \
    | sed 's% %\t%g' \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpend

# ---- Steps 3-6: Process each breakpoint (start and end) ----
export PATH=${main_dir}/software/SeqComplex:$PATH
cd ${cyto_sv_ml_dir}/software/SeqComplex

for bp in bpst bpend
do
    echo ${bp}

    # Step 3: Extract FASTA sequences from reference genome
    echo "# make bed file for SV breakpoints" && date
    ${py27_dir}/bedtools getfasta \
        -fi ${cyto_sv_ml_dir}/reference/hg38/hs38.fasta \
        -bed ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp} \
        -fo ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out

    # Step 4: Remove low-complexity sequences (long poly-N stretches)
    # Identify lines with >137 consecutive N's
    grep NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out -n \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.cr

    cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp

    # Mark low-complexity sequences and their headers for removal
    for i in $(awk '{print $1}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.cr | cut -d ":" -f 1)
    do
        echo $i
        i1=$((i-1))
        awk -v a="$i1" -v b="$i" '{
            if ((FNR==a)||(FNR==b)) {print $0"\tlc"}
            else {print $0}
        }' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp \
            > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc
        cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp
    done

    # Keep only non-low-complexity sequences
    awk '($NF!="lc"){print $0}' \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc

    # Steps 5-6: Run complexity analysis tools
    echo "# run SeqComplex and KZ for SV breakpoints" && date

    # SeqComplex: 24 sequence complexity metrics per sequence
    perl ${cyto_sv_ml_dir}/software/SeqComplex/profileComplexSeq.pl \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc

    # Komplexity (kz): 1 complexity metric per sequence
    kz --fasta < ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz

    # Join kz scores back to BED coordinates
    awk 'FNR==NR{a[$1]; b[$1]=$0; next}{
        c=$1":"$2"-"$3;
        if (c in a) {print $0"\t"b[c]}
    }' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp} \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz.index

    # Join SeqComplex scores to kz-indexed file
    awk 'FNR==NR{a[$1]; b[$1]=$0; next}($6 in a) {print $0"\t"b[$6]}' \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.complex \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz.index \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz.index_complex
done

# ---- Step 7: Consolidate start + end breakpoint complexity ----
echo "# consolidate SV breakpoint complexity information" && date

# Join end-bp and start-bp complexity by SV ID (column 5)
awk 'FNR==NR{a[$5];b[$5]=$0;next} ($5 in a){print $0"\t"b[$5]}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst.fa.out.kz.index_complex \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpend.fa.out.kz.index_complex \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.fa.out.kz.index_complex

# Extract final columns and add header
# Output: ID + 26 end-bp features (kz + SeqComplex) + 26 start-bp features
awk 'FNR==NR{a[$5];b[$5]=$0;next} ($6 in a){print $0"\t"b[$6]}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.fa.out.kz.index_complex \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed \
    | sed 's% %\t%g' \
    | cut -d$'\t' -f 6,14-15,17-40,48-49,51- \
    | awk '{
        if (FNR==1) {
            print "ID\tsv_bp_end_cc0\tsv_bp_end_cc1\tsv_bp_end_cc_v1\tsv_bp_end_cc_v2\tsv_bp_end_cc_v3\tsv_bp_end_cc_v4\tsv_bp_end_cc_v5\tsv_bp_end_cc_v6\tsv_bp_end_cc_v7\tsv_bp_end_cc_v8\tsv_bp_end_cc_v9\tsv_bp_end_cc_v10\tsv_bp_end_cc_v11\tsv_bp_end_cc_v12\tsv_bp_end_cc_v13\tsv_bp_end_cc_v14\tsv_bp_end_cc_v15\tsv_bp_end_cc_v16\tsv_bp_end_cc_v17\tsv_bp_end_cc_v18\tsv_bp_end_cc_v19\tsv_bp_end_cc_v20\tsv_bp_end_cc_v21\tsv_bp_end_cc_v22\tsv_bp_end_cc_v23\tsv_bp_end_cc_v24\tsv_bp_st_cc0\tsv_bp_st_cc1\tsv_bp_st_cc_v1\tsv_bp_st_cc_v2\tsv_bp_st_cc_v3\tsv_bp_st_cc_v4\tsv_bp_st_cc_v5\tsv_bp_st_cc_v6\tsv_bp_st_cc_v7\tsv_bp_st_cc_v8\tsv_bp_st_cc_v9\tsv_bp_st_cc_v10\tsv_bp_st_cc_v11\tsv_bp_st_cc_v12\tsv_bp_st_cc_v13\tsv_bp_st_cc_v14\tsv_bp_st_cc_v15\tsv_bp_st_cc_v16\tsv_bp_st_cc_v17\tsv_bp_st_cc_v18\tsv_bp_st_cc_v19\tsv_bp_st_cc_v20\tsv_bp_st_cc_v21\tsv_bp_st_cc_v22\tsv_bp_st_cc_v23\tsv_bp_st_cc_v24\n"$0
        } else {print $0}
    }' > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.kz.index_complex
