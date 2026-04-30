#!/bin/bash
###############################################################################
# sv_all_combine.sh
#
# Combines all SV feature sources into a single per-sample feature matrix:
#   1. Joins database annotations with VCF info features (by SV ID)
#   2. Appends breakpoint sequence complexity features
#   3. Appends caller support count (SUPP) from merged nonTRS VCF
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - sample: sample ID
#   $3 - size_k: SV size threshold in kb (e.g., 10)
#
# Output:
#   {sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp
#   (Final feature matrix with all annotations, info, complexity, and support)
###############################################################################

main_dir=$1
sample=$2
size_k=$3

echo ${sample}
echo "# Combine SV annotation, VCF info, complexity, and caller support" && date

# Step 1: Join database annotations (all_anno) with VCF info (all_info)
# Match on column 2 (SV ID) from both files
awk 'FNR==NR{a[$2];b[$2]=$0;next} ($2 in a) {print $0"\t"b[$2]}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_anno \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info

# Step 2: Append breakpoint sequence complexity features
# Match on column 1 (complexity file) to column 2 (combined file)
# If no complexity data found, fill with NAN placeholders (53 columns)
awk 'FNR==NR{a[$1];b[$1]=$0;next} {
    if ($2 in a) {print $0"\t"b[$2]}
    else {print $0"\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN"}
}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.kz.index_complex \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info.all_complex

# Step 3: Append caller support count (SUPP)
# Match SUPP from merged nonTRS VCF; TRS (MantaBND) SVs get SUPP=1
awk 'FNR==NR{a[$2];b[$2]=$3;next} {
    if (FNR==1) {print $0"\tSUPP"}
    else if ($2 in a) {print $0"\t"b[$2]}
    else if ($2~"MantaBND") {print $0"\t1"}
    else {print $1}
}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.sv_info.sim \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info.all_complex \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info.all_complex.supp
