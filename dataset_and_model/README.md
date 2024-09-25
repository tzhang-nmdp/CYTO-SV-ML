# Instructions for training, validation and prediction of SV classification using CYTO-SV-ML AUTOML models

### 1. Prepare CYTO-SV-ML environment and data
```
1.1 Install CYTO-SV-ML environment according to the instruction at main page

1.2 Create output folder: ${OUTPUT_DIR}/${cohort_name}/cyto_sv_ml/

1.3 Unzip and copy the consolidated dataset into the folder as: 
${OUTPUT_DIR}/test.sv.all.combine_all_trs
${OUTPUT_DIR}/test.sv.all.combine_all_nontrs
${OUTPUT_DIR}/test.sv.all.combine_all_trs_all
${OUTPUT_DIR}/test.sv.all.combine_all_nontrs_all

1.4 Unzip and copy pre-trained model and test data into the folder as:
${OUTPUT_DIR}/pre_trained_trs_model_folder
${OUTPUT_DIR}/pre_trained_nontrs_model_folder
${OUTPUT_DIR}/test_trs_sv_data_file
${OUTPUT_DIR}/test_nontrs_sv_data_file
```

### 2. Run CYTO-SV-ML AUTOML model training and validation
```
python ${CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML_modelling.py -s test -o ${OUTPUT_DIR}/test -x  ${CYTO_SV_ML_DIR}/sv_feature_metrics.csv -k 5 
```
Note: This is a subset of > * GB dataset and full data information will be available in our publication:

CYTO-SV-ML: a machine learning framework for discovery and classification of cytogenetic structural variants using whole genome sequencing data (under review).


### 3. Run CYTO-SV-ML AUTOML model prediction
 
```
from supervised.automl import AutoML
import pandas as pd
import pickle

# load pre-trained model
pre_trained_trs_model_folder=${OUTPUT_DIR}"/pre_trained_trs_model_folder"
pre_trained_automl_nontrs_model_folder=${OUTPUT_DIR}"/pre_trained_nontrs_model_folder"
pre_trained_automl_trs_model = AutoML(results_path=pre_trained_trs_model_folder)
pre_trained_automl_nontrs_model = AutoML(results_path=pre_trained_automl_nontrs_model_folder)

# load test dataset
test_trs_sv_data_file=${OUTPUT_DIR}"/test_trs_sv_data_file"
test_nontrs_sv_data_file=${OUTPUT_DIR}"/test_nontrs_sv_data_file"
test_trs_sv_data = pd.read_csv(test_trs_sv_data_file,sep="\t", header=0, index_col=None, keep_default_na=False)
test_nontrs_sv_data = pd.read_csv(test_nontrs_sv_data_file,sep="\t", header=0, index_col=None, keep_default_na=False)

# process data transformation
trs_transform_fit_pickle_file=${OUTPUT_DIR}"/trs_transform_fit_pickle_file"
nontrs_transform_fit_pickle_file=${OUTPUT_DIR}"/nontrs_transform_fit_pickle_file"
trs_transform_fit=pickle.load(open(trs_transform_fit_pickle_file, "rb"))
nontrs_transform_fit=pickle.load(open(trs_transform_fit_pickle_file,"rb"))
tf_test_trs_sv_data=trs_transform_fit.fit(test_trs_sv_data)
tf_test_nontrs_sv_data=nontrs_transform_fit.fit(test_nontrs_sv_data)

# run sv classification prediction
test_trs_sv_predictions = pre_trained_automl_trs_model.predict_all(tf_test_trs_sv_data)
test_nontrs_sv_predictions = pre_trained_automl_nontrs_model.predict_all(tf_test_nontrs_sv_data)
```
Note: The scripts in step 2 and 3 have been tested with python version 3.9.6
