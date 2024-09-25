# Instructions for training, validation and prediction of SV classification using CYTO-SV-ML AUTOML models

### 1. Prepare CYTO-SV-ML environment and data
```
1.1 Install CYTO-SV-ML environment according to the instruction at main page

1.2 Create output folder: ${OUTPUT_DIR}/${cohort_name}/cyto_sv_ml/

1.3 Unzip and copy the consolidated dataset into the folder as: 
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_trs
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_nontrs
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_trs_all
${OUTPUT_DIR}/${cohort_name}.sv.all.combine_all_nontrs_all

1.4 Unzip and copy pre-trained model and test data into the folder as:
${OUTPUT_DIR}/pre_trained_trs_model_folder
${OUTPUT_DIR}/pre_trained_nontrs_model_folder
${OUTPUT_DIR}/test_trs_sv_data_file
${OUTPUT_DIR}/test_nontrs_sv_data_file
```

### 2. Run CYTO-SV-ML AUTOML model training and validation
```
python {CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML_modelling.py -s ${cohort_name} -o ${OUTPUT_DIR}/${cohort_name} -x ${params.sv_feature_metrics_index} -k ${params.kfs} 
```
Note: This is a subset of > 5GB dataset and full data information will be available in our publication:

CYTO-SV-ML: a machine learning framework for discovery and classification of cytogenetic structural variants using whole genome sequencing data (under review).


### 3. Run CYTO-SV-ML AUTOML model prediction
 
```
from supervised.automl import AutoML
from sklearn import metrics, datasets
from sklearn.metrics import roc_curve, auc
import pandas as pd

# load pre-trained model
pre_trained_trs_model_folder=${OUTPUT_DIR}"/pre_trained_trs_model_folder"
pre_trained_automl_nontrs_model_folder=${OUTPUT_DIR}"/pre_trained_nontrs_model_folder"
pre_trained_automl_trs_model = AutoML(mode="Explain", algorithms=["Xgboost"], results_path=pre_trained_trs_model_folder)
pre_trained_automl_nontrs_model = AutoML(mode="Explain", algorithms=["Xgboost"], results_path=pre_trained_automl_nontrs_model_folder)

# load test dataset
test_trs_sv_data_file=${OUTPUT_DIR}"/test_trs_sv_data_file"
test_trs_sv_data_file=${OUTPUT_DIR}"/test_trs_sv_data_file"
test_trs_sv_data = pd.read_csv(test_trs_sv_data_file,sep="\t", header=0, index_col=None, keep_default_na=False)
test_nontrs_sv_data = pd.read_csv(test_trs_sv_data_file,sep="\t", header=0, index_col=None, keep_default_na=False)

# run sv classification prediction
test_trs_sv_predictions = pre_trained_automl_trs_model.predict_all(test_trs_sv_data)
test_nontrs_sv_predictions = pre_trained_automl_nontrs_model.predict_all(test_nontrs_sv_data)

```
Note: The scripts in step 2 and 3 have been tested with python version 3.9.6
