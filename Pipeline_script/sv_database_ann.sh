#!/bin/bash
###############################################################################
# sv_database_ann.sh
#
# Annotates SVs against public and internal SV databases.
# For each database, computes overlap (nonTRS) or breakpoint distance (TRS).
#
# Steps:
#   1. Initialize annotation files with SV coordinates
#   2. Annotate against sample-specific uwstl_s database (ChromoSeq calls)
#   3. Annotate against each public SV database (1000G, gnomAD, COSMIC, etc.)
#   4. Consolidate all annotations into a single file per sample
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - sample: sample ID
#   $4 - sv_db_vector: '@'-delimited list of SV database names
#   $5 - py27_dir: path to Python 2.7 bin directory
#   $6 - size_k: SV size threshold in kb (e.g., 10)
#
# Output:
#   {sample}.{size_k}k.sv.all.sv_id_mapping.all_anno
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
sv_db_vector=$4
py27_dir=$5
size_k=$6

# Parse '@'-delimited database list
echo ${sv_db_vector} | sed "s%@%\n%g" > ${main_dir}/out/${sample}/sv_db_vector.tmp
sd_ln=$(wc -l ${main_dir}/out/${sample}/sv_db_vector.tmp | awk '{print $1}')

# ---- Step 1: Initialize annotation files ----
# Create header + coordinate columns for TRS and nonTRS annotation files
echo "# initiate the touch of SV anno files" && date
awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs_anno

awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $5"\t"$1"\t"$2"\t"$3"\t"$1"\t"$4}}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs_anno

# ---- Step 2: Annotate against uwstl_s (ChromoSeq sample-specific SVs) ----
echo "# SV database annotation label" && date
for SV_database_name in uwstl_s
do
    echo ${SV_database_name} "ok" && date

    # NonTRS: use overlap-based mapping if database file exists
    if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${SV_database_name}.nontrs.gz ]; then
        ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_database_mapping.py \
            -i ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs \
            -t ${main_dir}/out/${sample}/sv_caller_results/${sample}.${SV_database_name}.nontrs.gz \
            -d 1000 -p 0.7 \
            -o ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs.${SV_database_name}
    else
        # No database file: fill with NAN
        awk 'FNR!=1{$4=$1"\t"$4; print $0"\tNAN"}' \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs \
            | sed 's% %\t%g' \
            > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs.${SV_database_name}
    fi

    # TRS: use breakpoint-distance mapping if database file exists
    if [ -s ${main_dir}/out/${sample}/${sample}/sv_caller_results/${sample}.${SV_database_name}.trs ]; then
        ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_bnd_database_mapping.py \
            ${main_dir}/out/${sample}/sv_caller_results/${sample}.${SV_database_name}.trs \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs \
            ${SV_database_name}_1000
    else
        awk 'FNR!=1{print $0"\tNAN"}' \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs \
            | sed 's% %\t%g' \
            > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs.${SV_database_name}_1000
    fi

    # Consolidate annotation into the _anno file
    ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_db_tf.py \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs ${SV_database_name} f
    ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_db_tf.py \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs ${SV_database_name}_1000 f
done

# ---- Step 3: Annotate against each public SV database ----
for i in $(seq 1 $sd_ln)
do
    SV_database_name=$(awk -v a="$i" '(FNR==a){print $1}' ${main_dir}/out/${sample}/sv_db_vector.tmp)
    echo ${SV_database_name} "ok" && date

    # NonTRS annotation: overlap-based mapping via tabix
    if [ -s ${cyto_sv_ml_dir}/SV_database/${SV_database_name}.nontrs.gz ]; then
        ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_database_mapping.py \
            -i ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs \
            -t ${cyto_sv_ml_dir}/SV_database/${SV_database_name}.nontrs.gz \
            -d 1000 -p 0.7 \
            -o ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs.${SV_database_name}
    else
        awk 'FNR!=1{$4=$1"\t"$4; print $0"\tNAN"}' \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs \
            | sed 's% %\t%g' \
            > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs.${SV_database_name}
    fi

    # TRS annotation: breakpoint-distance mapping
    # Use .bp.trs (with CI info) if available, otherwise .trs (simple format)
    if [ -s ${cyto_sv_ml_dir}/SV_database/${SV_database_name}.bp.trs ]; then
        ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_bnd_database_mapping.bp.py \
            ${cyto_sv_ml_dir}/SV_database/${SV_database_name}.bp.trs \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs \
            ${SV_database_name}_1000
    elif [ -s ${cyto_sv_ml_dir}/SV_database/${SV_database_name}.trs ]; then
        ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_bnd_database_mapping.py \
            ${cyto_sv_ml_dir}/SV_database/${SV_database_name}.trs \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs \
            ${SV_database_name}_1000
    else
        awk 'FNR!=1{print $0"\tNAN"}' \
            ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs \
            | sed 's% %\t%g' \
            > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs.${SV_database_name}_1000
    fi

    # Consolidate annotation into the _anno file
    ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_db_tf.py \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs ${SV_database_name} f
    ${py27_dir}/python ${cyto_sv_ml_dir}/Pipeline_script/sv_db_tf.py \
        ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs ${SV_database_name}_1000 f
done

# ---- Step 4: Consolidate all annotations ----
# Merge TRS and nonTRS annotations, then join with SV ID mapping
echo "# prepare SV DB annotation file" && date
cat ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.trs_anno \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.nontrs_anno \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno

awk 'FNR==NR{a[$1];b[$1]=$0;next} ($1 in a) {print $0"\t"b[$1]}' \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno \
    ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping \
    > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_anno

# Clean up
rm -rf ${main_dir}/out/${sample}/sv_db_vector.tmp
