#!/bin/bash
###############################################################################
# sv_info_extract.sh
#
# Extracts read-level features from SVTyper-genotyped VCFs and maps them
# back to consolidated SV IDs. Iterates over all callers and joins their
# info columns into a single wide-format file.
#
# Features extracted include: PR, SR, GT, CN, BND_DEPTH, MATE_BND_DEPTH,
# CIPOS, CIEND, SVTYPE, END, CHR2, etc.
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - sample: sample ID
#   $3 - sv_caller_vector: '@'-delimited list of SV callers
#   $4 - size_k: SV size threshold in kb (e.g., 10)
#
# Output:
#   {sample}.{size_k}k.sv.all.sv_id_mapping.all_info
###############################################################################

main_dir=$1
sample=$2
sv_caller_vector=$3
size_k=$4

echo "# starting SV info extraction" && date

# Start with the SV ID mapping file as the base
cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t

# Parse '@'-delimited caller list
echo ${sv_caller_vector} | sed "s%@%\n%g" > ${main_dir}/out/${sample}/sv_caller_vector.tmp2
sc_ln=$(wc -l ${main_dir}/out/${sample}/sv_caller_vector.tmp2 | awk '{print $1}')

# Iteratively join each caller's info to the mapping file
n=0
for i in $(seq 1 $sc_ln)
do
    sv_caller=$(awk -v a="$i" '(FNR==a){print $1}' ${main_dir}/out/${sample}/sv_caller_vector.tmp2)

    # Try SVTyper-genotyped info first
    if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim ]; then
        echo ${sv_caller} "svtype ok" && date
        if [ "${n}" == 0 ]; then
            # First caller: join including header
            awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' \
                ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim \
                ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
                > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t
            n=1
        else
            # Subsequent callers: skip header row
            awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' \
                ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim \
                ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
                > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t
        fi
    else
        # Fallback: use original (non-svtyped) VCF info
        echo ${sv_caller} "svtype pass" && date
        if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.sv_info.sim ]; then
            cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.sv_info.sim \
                ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim.tmp
            if [ "${n}" == 0 ]; then
                awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' \
                    ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.svtyped.vcf.sv_info.sim.tmp \
                    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
                    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t
                n=1
            else
                awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' \
                    ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.svtyped.vcf.sv_info.sim.tmp \
                    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
                    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t
            fi
        else
            echo ${sv_caller} "all pass" && date
        fi
    fi

    # Update the running file for the next iteration
    cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t
done

# ---- Validate: check all rows have the same number of columns ----
echo "#check the complete status of sv_info extraction" && date
awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t | awk '!a[$1]++'
info_check=$(awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t | awk '!a[$1]++' | wc -l)

if (( $info_check > 1 )); then
    # Some SVs are missing info columns — pad with placeholder dots
    echo "sv_info missing !!!"
    awk '{if (NF==2){print $0}}' \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info_missing

    awk '{if (NF==2){print $0"\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."} else{print $0}}' \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
        > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info

    # Verify padding fixed the issue
    info_check2=$(awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++' | wc -l)
    if (( $info_check2 > 1 )); then
        awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++'
        echo "sv_info2 missing !!!"
        exit
    fi
else
    # All rows have consistent column count
    cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info
fi

# Clean up temporary files
rm ${main_dir}/out/${sample}/${sample}.*tmp*
rm ${main_dir}/out/${sample}/sv_caller_vector.tmp2
