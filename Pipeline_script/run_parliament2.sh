#!/bin/bash
###############################################################################
# run_parliament2.sh
#
# Runs the Parliament2 multi-caller SV pipeline via Docker.
# Parliament2 integrates: Breakdancer, CNVnator, Delly (del/dup/inv).
#
# Steps:
#   1. Copy BAM and reference files to expected locations
#   2. Run Parliament2 Docker container with selected callers
#   3. Clean up temporary reference copies
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - parliament_docker: Docker image name for Parliament2
#   $4 - sample: sample ID
#
# Input:
#   {main_dir}/in/{sample}.bam
#   {main_dir}/in/{sample}.bam.bai
#
# Output:
#   {sample}.breakdancer.vcf
#   {sample}.cnvnator.vcf
#   {sample}.delly.deletion.vcf
#   {sample}.delly.duplication.vcf
#   {sample}.delly.inversion.vcf
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
parliament_docker=$3
sample=$4

echo "# run parliament docker" && date

# Copy BAM to expected input location (Parliament2 expects specific paths)
cp ${main_dir}/in/${sample}.bam ${main_dir}/in/input.bam
cp ${main_dir}/in/${sample}.bam.bai ${main_dir}/in/input.bam.bai

# Copy hg38 reference to input directory (required by Parliament2)
sudo mkdir -p ${main_dir}/in/reference/hg38/
sudo chmod -R 777 ${main_dir}/in/reference/hg38/
cp ${cyto_sv_ml_dir}/reference/hg38/hs38.fasta ${main_dir}/in/reference/hg38/hs38.fasta
cp ${cyto_sv_ml_dir}/reference/hg38/hs38.fasta.fai ${main_dir}/in/reference/hg38/hs38.fasta.fai

# Run Parliament2 Docker with selected SV callers
# --filter_short_contigs: skip alt/random contigs
# Enabled callers: breakdancer, cnvnator, delly (deletion, inversion, duplication)
# Note: Lumpy is NOT enabled (requires specific library prep metadata)
sudo docker run --rm --privileged \
    -v ${main_dir}/in/:/home/dnanexus/in \
    -v ${main_dir}/out/${sample}/:/home/dnanexus/out \
    ${parliament_docker} \
    --bam ${sample}.bam \
    --bai ${sample}.bam.bai \
    -r reference/hg38/hs38.fasta \
    --fai reference/hg38/hs38.fasta.fai \
    --prefix ${sample} \
    --filter_short_contigs \
    --breakdancer \
    --cnvnator \
    --delly_deletion \
    --delly_inversion \
    --delly_duplication

# Clean up temporary reference copy
sudo rm -rf ${main_dir}/in/reference/
