#!/bin/bash
###############################################################################
# all_sample_sv_combine.sh
# 
# Combines per-sample SV feature matrices into a single cohort-level file.
# Adds a 'sample_id' column to identify which sample each SV belongs to.
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - sample_all: cohort name (used for output naming)
#   $3 - SAMPLES_vector: file containing '@'-delimited sample IDs
#
# Output:
#   {main_dir}/out/{sample_all}/{sample_all}.sv.all.combine_all
###############################################################################

main_dir=$1
sample_all=$2
SAMPLES_vector=$3

echo "# Combine all sample SV data" && date

# Create cohort output directory
sudo rm -rf ${main_dir}/out/${sample_all}/
sudo mkdir ${main_dir}/out/${sample_all}/
sudo chmod 777 -R ${main_dir}/out/${sample_all}/

# Parse '@'-delimited sample list into one sample per line
sed "s%@%\n%g" ${SAMPLES_vector} > ${SAMPLES_vector}.tmp
sm_ln=$(wc -l ${SAMPLES_vector}.tmp | awk '{print $1}')

# Iterate over each sample and concatenate feature files
for i in $(seq 1 $sm_ln)
do
    sample=$(awk -v a="$i" '(FNR==a){print $1}' ${SAMPLES_vector}.tmp)
    sample_file="${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex.supp"
    echo ${sample} ${sample_file}

    if (( i==1 )); then
        # First sample: include header row with 'sample_id' column appended
        awk -v a="$sample" '{
            if (FNR==1) {print $0"\tsample_id"}
            else {print $0"\t"a}
        }' ${sample_file} > ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all.tmp
    else
        # Subsequent samples: skip header, append with sample_id
        awk -v a="$sample" '(FNR>1){print $0"\t"a}' \
            ${sample_file} >> ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all.tmp
    fi
done

# Verify the combined file exists
ls -a ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all.tmp

# Replace spaces with tabs for consistent formatting
sed 's% %\t%g' ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all.tmp \
    > ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all

# Clean up temporary files
sudo rm -rf ${SAMPLES_vector}.tmp
sudo rm -rf ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all.tmp
