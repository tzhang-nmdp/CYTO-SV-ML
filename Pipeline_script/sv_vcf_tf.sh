#!/bin/bash
###############################################################################
# sv_vcf_tf.sh
#
# Transforms raw SV VCF files from each caller:
#   1. Standardizes SV IDs using sv_id_tf.py
#   2. Deduplicates SVs within each caller using SURVIVOR merge
#   3. Filters SVs by minimum size using sv_size.py
#   4. Splits into TRS (BND) and nonTRS (DEL/DUP/INV) using sv_vcf_tf.py
#   5. Extracts simplified VCF info using sv_info_tf_sim.py
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - sample: sample ID
#   $4 - size: minimum SV size in bp (e.g., 10000)
#
# Output per caller:
#   {sample}.{caller}.vcf.{size_k}k.trs_tf      - TRS SVs (transformed)
#   {sample}.{caller}.vcf.{size_k}k.nontrs_tf    - nonTRS SVs (transformed)
#   {sample}.{caller}.vcf.{size_k}k.sv_info.sim  - Simplified VCF info
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
size=$4
size_k=$((size/1000))

echo "# Prepare SV VCF files with size restriction and extract sv_vcf_info" && date

# Process each SV caller's VCF
for sv_caller in breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta
do
    echo ${sv_caller}

    # Step 1: Standardize SV IDs (add chr:pos:end:chr2:type format)
    python ${cyto_sv_ml_dir}/Pipeline_script/sv_id_tf.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf c

    # Step 2: Deduplicate SVs within the caller using SURVIVOR
    ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.list
    SURVIVOR merge \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.list \
        1000 1 1 0 0 10 \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s

    # Map deduplicated IDs back to original IDs
    awk '($1!~"#"){split($3,b,":");print b[1]":"b[2]":"b[3]":"b[4]":"b[5]"\t"$3}' \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s.id

    awk 'FNR==NR{a[$1];c[$1]=$2;next}{split($3,b,":"); e=b[1]":"b[2]":"b[3]":"b[4]":"b[5]; if (($1!~"#")&&(e in a)) {$3=c[e]; print $0}}' \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s.id \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id \
        | sed 's% %\t%g' \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.s

    # Reconstruct VCF with header + deduplicated body
    awk '($1~"#"){print $0}' \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.hd
    cat ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.hd \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.s \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id

    # Step 3: Filter by minimum SV size
    echo ${sv_caller} "SV size tf" && date
    python ${cyto_sv_ml_dir}/Pipeline_script/sv_size.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id \
        $size down \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k

    # Step 4: Split into TRS and nonTRS VCFs
    echo ${sv_caller} "SV type tf" && date
    python ${cyto_sv_ml_dir}/Pipeline_script/sv_vcf_tf.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k

    # Step 5: Extract simplified VCF info
    echo ${sv_caller} "SV info tf" && date
    python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
done
