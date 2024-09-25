# Instructions for training, validation and prediction of SV classification using CYTO-SV-ML AUTOML models

### 1. Prepare CYTO-SV-ML environment and dataset
1.1 Install CYTO-SV-ML environment according to the main page

1.2 Create output folder: ${OUTPUT_DIR}/${cohort_name}/cyto_sv_ml/

1.3 Unzip and copy the consolidated datasets into the folder as: 
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_trs
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_nontrs
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_trs_all
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_nontrs_all

1.4 Copy pre-trained model and test data to into the folder as:
${OUTPUT_DIR}/pre_trained_model.pickle
${OUTPUT_DIR}/test_data.csv

### 2. Run CYTO-SV-ML AUTOML model training and validation
```
python {CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML_modelling.py -s ${cohort_name} -o ${OUTPUT_DIR}/${cohort_name} -x ${params.sv_feature_metrics_index} -k ${params.kfs} 
```
Note: This is a subset of > 5GB dataset and full data information will be released in our publication:

CYTO-SV-ML: a machine learning framework for discovery and classification of cytogenetic structural variants using whole genome sequencing data (under review).

### 3. Run CYTO-SV-ML AUTOML model prediction
 
```
from supervised.automl import AutoML
from sklearn import metrics, datasets
from sklearn.metrics import roc_curve, auc
import pandas as pd
import pickle

# load pre-trained model
pre_trained_trs_model_folder=${OUTPUT_DIR}"/pre_trained_trs_model_folder"
pre_trained_automl_nontrs_model_folder=${OUTPUT_DIR}"/pre_trained_nontrs_model_folder"
pre_trained_automl_trs_model = AutoML(mode="Explain", algorithms=["Xgboost"], n_jobs= 6, results_path=pre_trained_trs_model_folder)
pre_trained_automl_nontrs_model = AutoML(mode="Explain", algorithms=["Xgboost"], n_jobs= 6, results_path=pre_trained_automl_nontrs_model_folder)

# load test dataset
test_trs_sv_data_file=${OUTPUT_DIR}"/test_data.csv"
test_trs_sv_data = pd.read_csv(test_sv_data_file,sep="\t", header=0, index_col=None, keep_default_na=False)
test_trs_sv_data_file=${OUTPUT_DIR}"/test_data.csv"
test_trs_sv_data = pd.read_csv(test_sv_data_file,sep="\t", header=0, index_col=None, keep_default_na=False)

# run sv classification prediction
predictions = pre_trained_automl_trs_model.predict_all(sv_data_10_tf)
predictions = pre_trained_automl_nontrs_model.predict_all(sv_data_10_tf)

# run model performance evaluation


```
Note: The scripts in step 2 and 3 has been tested with python version 3.9.6
