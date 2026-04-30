#!/bin/bash
###############################################################################
# svtyper_qc.sh
#
# Runs SVTyper on each caller's VCF to extract read-level genotype evidence.
# SVTyper computes paired-read (PR) and split-read (SR) support for each SV.
#
# Steps per caller:
#   1. Transform TRS VCF for SVTyper compatibility (MATEID fix)
#   2. Run SVTyper on TRS and nonTRS VCFs separately
#   3. Re-ID TRS SVs with chr:pos:end:chr2:BND format
#   4. Combine TRS + nonTRS svtyped VCFs
#   5. Standardize SV IDs and extract simplified info
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - sample: sample ID
#   $4 - sv_caller_vector: '@'-delimited list of SV callers
#   $5 - size_k: SV size threshold in kb (e.g., 10)
#
# Output per caller:
#   {sample}.{caller}.{size_k}k.all.svtyped.vcf          - SVTyper genotyped VCF
#   {sample}.{caller}.{size_k}k.all.svtyped.vcf.sv_info.sim - Simplified info
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
sv_caller_vector=$4
size_k=$5

# SVTyper command (uses 8 cores for split-read extraction)
svtyper="svtyper-sso --core 8"

# Parse '@'-delimited caller list
echo ${sv_caller_vector} | sed "s%@%\n%g" > sv_caller_vector.tmp
sc_ln=$(wc -l sv_caller_vector.tmp | awk '{print $1}')

echo "# run svtyper for all sv callers" && date
for i in $(seq 1 $sc_ln)
do
    sv_caller=$(awk -v a="$i" '(FNR==a){print $1}' sv_caller_vector.tmp)
    echo ${sv_caller}

    # ---- Step 1: Prepare TRS VCF for SVTyper ----
    # SVTyper requires MATEID to be properly formatted for BND records
    python ${cyto_sv_ml_dir}/Pipeline_script/trs_svtyper_tf.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.trs_tf

    # ---- Step 2: Run SVTyper on TRS and nonTRS separately ----
    # TRS (translocation) SVs
    ${svtyper} --max_reads 100000 \
        -i ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.trs_tf.tmp \
        -B ${main_dir}/in/${sample}.bam \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf

    # NonTRS (DEL/DUP/INV) SVs
    ${svtyper} --max_reads 100000 \
        -i ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.nontrs_tf \
        -B ${main_dir}/in/${sample}.bam \
        > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.nontrs_tf.svtyped.vcf

    # ---- Step 3: Re-ID TRS SVs with standardized BND format ----
    if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf ]; then
        # Add chr:pos:end:chr2:BND:ID format to TRS SVs
        awk '{if ($1~"#") {print $0} else if ($5~"\\["){split($5,a,"["); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0} else if ($5~"\\]"){split($5,a,"]"); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0}}' \
            ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf \
            > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.re_id

        # Extract body only (no header) for concatenation
        awk '($1!~"#"){if ($5~"\\["){split($5,a,"["); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0} else if ($5~"\\]"){split($5,a,"]"); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0}}' \
            ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf \
            | sed 's% %\t%g' \
            > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp
    fi

    # ---- Step 4: Combine TRS + nonTRS svtyped VCFs ----
    if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp ]; then
        # Concatenate nonTRS (with header) + TRS (body only)
        cat ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.nontrs_tf.svtyped.vcf \
            ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp \
            > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf
        sudo rm -rf ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp
    else
        # No TRS SVs: use nonTRS only
        cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.nontrs_tf.svtyped.vcf \
            ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf
    fi

    # ---- Step 5: Standardize IDs and extract simplified info ----
    python ${cyto_sv_ml_dir}/Pipeline_script/sv_id_tf.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf c
    sudo mv ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.re_id \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf

    python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py \
        ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf
done

# Clean up
sudo rm -rf sv_caller_vector.tmp

# Remove BAM files after SVTyper is done (saves disk space)
sudo rm -rf ${main_dir}/in/${sample}.bam*
