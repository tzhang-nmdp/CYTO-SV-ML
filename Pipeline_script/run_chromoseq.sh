#!/bin/bash
###############################################################################
# run_chromoseq.sh
#
# Runs the ChromoSeq clinical WGS pipeline via Docker to call SVs using
# Manta (translocations/BND) and ichorCNV (copy number / DEL/DUP).
#
# Steps:
#   1. Create output directories
#   2. Generate ChromoSeq input JSON with sample-specific values
#   3. Run ChromoSeq Docker container (Cromwell + WDL workflow)
#   4. Extract ChromoSeq-called SVs as uwstl_s database (for annotation)
#   5. Separate Manta VCF and transform ichorCNV segments to VCF format
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - chromoseq_docker: Docker image name for ChromoSeq
#   $4 - sample: sample ID
#   $5 - gender: sample sex (Male/Female)
#
# Input:
#   {main_dir}/in/{sample}.cram
#
# Output:
#   {sample}.manta.vcf      - Manta SV calls (BND + DEL/DUP/INV)
#   {sample}.ichnorcnv.vcf   - ichorCNV copy number calls (DEL/DUP)
#   {sample}.uwstl_s.trs     - ChromoSeq TRS SVs (for database annotation)
#   {sample}.uwstl_s.nontrs.gz - ChromoSeq nonTRS SVs (bgzipped + indexed)
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
chromoseq_docker=$3
sample=$4
gender=$5

echo $sample $gender

# ---- Step 1: Create output directories ----
sudo mkdir -p ${main_dir}/out/${sample}/sv_caller_results/
sudo chmod 777 -R ${main_dir}/out/${sample}/sv_caller_results/

# ---- Step 2: Generate ChromoSeq input JSON ----
# Replace placeholder sample ID and gender in the template JSON
echo "# set up the configure file for chromoseq" && date
sed "s%XXXXXX%${sample}%g" \
    ${cyto_sv_ml_dir}/software/docker-basespace_chromoseq/lsf/inputs.json \
    | sed "s%Male%${gender}%g" \
    > ${cyto_sv_ml_dir}/software/docker-basespace_chromoseq/lsf/inputs.json.tmp

# ---- Step 3: Run ChromoSeq Docker ----
# Mounts main_dir as /scratch and cyto_sv_ml_dir as /cyto_sv_ml_dir
# Runs Cromwell WDL workflow inside the container
echo "# run chromoseq docker" && date
sudo docker run --rm --privileged \
    -v ${main_dir}/:/scratch \
    -v ${cyto_sv_ml_dir}/:/cyto_sv_ml_dir \
    --entrypoint /bin/sh ${chromoseq_docker} -c \
    '/usr/bin/java \
        -Dconfig.file=/cyto_sv_ml_dir/software/docker-basespace_chromoseq/lsf/application.new.conf \
        -jar /opt/cromwell-36.jar run -t wdl \
        -i /cyto_sv_ml_dir/software/docker-basespace_chromoseq/lsf/inputs.json.tmp \
        /cyto_sv_ml_dir/software/docker-basespace_chromoseq/workflow_files/Chromoseq.v17.wdl'

# ---- Step 4: Extract ChromoSeq SVs as uwstl_s database ----
# These are used later as a sample-specific SV database for annotation
echo "# prepare the SV db uwstl_s file for chromoseq pipeline" && date

# NonTRS SVs (DEL/DUP/INV/INS): extract from chromoseq.txt, bgzip + index
awk '($1=="DUP")||($1=="DEL")||($1=="INS")||($1=="INV"){
    print $2"\t"$3"\t"$5"\t"$1"\t"$2":"$3":"$5":"$1
}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.chromoseq.txt \
    > ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.nontrs

sort -k 1,1 -k 2,3n \
    ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.nontrs \
    | bgzip -f \
    > ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.nontrs.gz
tabix -p vcf ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.nontrs.gz

# TRS SVs (BND): extract from chromoseq.txt
awk '($1=="BND"){
    print $2"\t"$3"\t"$5"\t"$4"\t"$2":"$3":"$5":"$4":"$1
}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.chromoseq.txt \
    > ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.trs

# ---- Step 5: Separate Manta VCF and transform ichorCNV ----
echo "# prepare the SV vcf files (manta + ichnorcnv)" && date

# Decompress annotated VCF from ChromoSeq output
gunzip -f ${main_dir}/out/${sample}/sv_caller_results/${sample}.svs_annotated.vcf.gz

# Extract Manta calls (exclude ichorCNV lines marked with "ICHOR:")
grep -v ICHOR: \
    ${main_dir}/out/${sample}/sv_caller_results/${sample}.svs_annotated.vcf \
    > ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf

# Transform ichorCNV segment file (.segs.txt) into VCF format
python ${cyto_sv_ml_dir}/Pipeline_script/ichnorcnv_tf.py \
    ${main_dir}/out/${sample}/sv_caller_results/${sample}.segs.txt \
    ${main_dir}/out/${sample}/sv_caller_results/${sample}.ichnorcnv.vcf
