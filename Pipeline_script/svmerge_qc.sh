#!/bin/bash
###############################################################################
# svmerge_qc.sh
#
# Merges SVs across all callers using SURVIVOR and extracts metadata:
#   1. Merges nonTRS SVs across callers (breakpoint distance <= 1000bp)
#   2. Copies TRS SVs from Manta (only caller producing BND calls)
#   3. Merges ALL SVs (TRS + nonTRS) into a consolidated VCF
#   4. Extracts caller support count (SUPP) per SV
#   5. Creates SV ID mapping (consolidated ID <-> per-caller ID)
#   6. Generates simplified SV representations (TRS/nonTRS bed-like files)
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - sample: sample ID
#   $4 - size_k: SV size threshold in kb (e.g., 10)
#
# Output:
#   {sample}.{size_k}k.sv.all              - Consolidated VCF (all SVs)
#   {sample}.{size_k}k.sv.all.sv_id_mapping - ID mapping file
#   {sample}.{size_k}k.sv.all.trs           - Simplified TRS SV file
#   {sample}.{size_k}k.sv.all.nontrs        - Simplified nonTRS SV file
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
size_k=$4

# Step 1: Merge nonTRS SVs across callers using SURVIVOR
echo "# Merge all non/trs SV" && date
ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf.${size_k}k.nontrs_tf \
    | grep -v svtype \
    > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.nontrs_tf.list

# SURVIVOR merge: max distance=1000bp, min callers=1, same type=1, same strand=0, min size=0, min size=10
SURVIVOR merge \
    ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.nontrs_tf.list \
    1000 1 1 0 0 10 \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all

# TRS SVs: only Manta produces BND calls, so just copy
cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf.${size_k}k.trs_tf \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all

# Step 2: Merge ALL SVs (TRS + nonTRS) into one consolidated VCF
echo "# Merge all SV" && date
ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf.${size_k}k.*trs_tf \
    | grep -v svtype \
    > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.all.list

SURVIVOR merge \
    ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.all.list \
    1000 1 1 0 0 10 \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all

# Step 3: Extract caller support (SUPP) from merged nonTRS VCF
echo "# SV caller SUPP info extraction" && date
python ${cyto_sv_ml_dir}/Pipeline_script/sv_consolidate_info_tf_sim.py \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all SUPP

# Step 4: Create SV ID mapping (consolidated <-> per-caller IDs)
python ${cyto_sv_ml_dir}/Pipeline_script/sv_consolidate_id_mapping.py \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all

# Step 5: Generate simplified SV representations
echo "# SV vcf simplified transformation" && date

# Extract simplified VCF info (END, CHR2, CIPOS, CIEND, etc.)
python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all

# Split into TRS and nonTRS bed-like files
python ${cyto_sv_ml_dir}/Pipeline_script/sv_vcf_sim.py \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all
