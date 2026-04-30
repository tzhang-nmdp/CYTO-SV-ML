#!/bin/bash
###############################################################################
# interface_docker.sh
#
# Prepares data and builds a Docker image for the R Shiny web portal.
#
# Steps:
#   1. Run interface_data.py to predict SV classes using the last k-fold model
#   2. Run interface_data.R to prepare RData for the Shiny app
#   3. Build Docker image containing the Shiny app
#
# Arguments:
#   $1 - main_dir: root working directory
#   $2 - cyto_sv_ml_dir: CYTO-SV-ML repository path
#   $3 - sample_all: cohort name
#   $4 - cyto_band_file: path to cytoband dictionary CSV
#   $5 - k: number of k-folds (uses last fold k-1 for final model)
#
# Output:
#   {sample_all}.all_pred              - Combined TRS + nonTRS predictions
#   cyto_sv_ml.RData                   - R data file for Shiny app
#   Docker image: cyto-sv-ml-app:{sample_all}
###############################################################################

main_dir=$1
cyto_sv_ml_dir=$2
sample_all=$3
cyto_band_file=$4
k=$5

# Use the last k-fold model (0-indexed: k-1)
k=$((k-1))

echo ${main_dir} ${cyto_sv_ml_dir} ${sample_all}

# ---- Step 1: Generate SV predictions ----
# Uses the last k-fold TRS and nonTRS models to predict all SVs
echo "# prepare the data file for cyto-sv-ml application" && date
python ${cyto_sv_ml_dir}/Pipeline_script/interface_data.py \
    -t ${main_dir}'/out/'${sample_all}'/cyto_sv_ml/'${sample_all}'_trs_'${k} \
    -n ${main_dir}'/out/'${sample_all}'/cyto_sv_ml/'${sample_all}'_nontrs_'${k} \
    -i ${main_dir}'/out/'${sample_all}'/'${sample_all}

# ---- Step 2: Prepare RData for Shiny app ----
# Combines predictions with cytoband coordinates for genome visualization
Rscript ${cyto_sv_ml_dir}/Pipeline_script/interface_data.R \
    -i ${main_dir}'/out/'${sample_all}'/'${sample_all} \
    -o ${cyto_sv_ml_dir}"/docker/cyto-sv-ml/data/cyto_sv_ml.RData" \
    -c ${cyto_band_file}

# ---- Step 3: Build Docker image ----
cd ${cyto_sv_ml_dir}"/docker"
echo "# build R-shiny based user interface in docker container" && date
sudo docker build -t cyto-sv-ml-app:${sample_all} .

# To run the Docker image locally:
# sudo docker run -d -p 8000:80 cyto-sv-ml-app:${sample_all}
# Then open http://localhost:8000/ in a browser
