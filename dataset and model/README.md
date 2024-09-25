# Instructions for training, validation and prediction of SV classification using CYTO-SV-ML AUTOML models

### 1. Prepare CYTO-SV-ML environment and dataset
1.1 Install CYTO-SV-ML environment according to the main page

1.2 Create output folder: ${OUTPUT_DIR}/${cohort_name}/cyto_sv_ml/

1.3 Unzip and copy the consolidated datasets into the folder as: 
${OUTPUT_DIR}/${cohort_name}".sv.all.combine_all_trs" 
${OUTPUT_DIR}/${cohort_name}".sv.all.combine_all_nontrs"
${OUTPUT_DIR}/${cohort_name}".sv.all.combine_all_trs_all" 
${OUTPUT_DIR}/${cohort_name}".sv.all.combine_all_nontrs_all"

### 2. Run CYTO-SV-ML AUTOML model training and validation
```
python {CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML_modelling.py -s ${cohort_name} -o ${OUTPUT_DIR}/${cohort_name} -x ${params.sv_feature_metrics_index} -k ${params.kfs} 
```
Note: This is an example dataset (a subset of > 5GB full dataset) and full data information will be released in our publication:

CYTO-SV-ML: a machine learning framework for discovery and classification of cytogenetic structural variants using whole genome sequencing data (under review).

### 3. Run CYTO-SV-ML AUTOML model prediction
Copy pre-trained model to ${OUTPUT_DIR}/  
```
# load pre-trained model

# load test dataset

# run prediction


```
Note: The scripts in step 2 and 3 has been tested with python version 3.9.6
