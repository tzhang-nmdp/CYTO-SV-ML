# CYTO-SV-ML
## Introduction
CYTO-SV-ML is a snakemake based workflow for calling, quality control, annotations, and classification (germline/somatic) of strucutural variations from whole genome sequencing data. It's integrated multiple SV callers (Delly/Breakdancer/CNVnator/Manta/Lumpy/ichorCNV) from Parliament2 and Chromoseq pipelines. Empowered by a collection of Public SV databases (DGV (http://dgv.tcag.ca), gnomAD (gnomad.broadinstitute.org/help/sv-overview), 1000 Genomes (www.internationalgenome.org/phase-3-structural-variant-dataset), CytoAtlas (atlasgeneticsoncology.org), and COSMIC (cancer.sanger.ac.uk/cosmic)), benchmarking SV data have been produced via comprehensive annotation and filtering algorithms. IN two MDS patient cohorts (494/94), it's achieved > 0.90 AUC and retrieved cytogenetic informations in 11% patients whose clinical records are missing. 

#### CYTO-SV-ML pipeline design:
![CYTO-SV-ML PIPELINE WORKFLOW](workflow.png)

### CYTO-SV-ML snakemake workflow
![CYTO-SV-ML Snakemake Snakemake workflow](cyto-sv-ml_snakemake_workflow.png)

<!--ts-->
   * [Installation](#Installation)
   * [Pipeline Usage](#Pipeline-usage)
   * [Download SV related Resource](#Download-SV-related-Resource)

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Installation
============
### Install Conda/Docker environment
```
Conda and Docker are prerequisites and please install them according to official websites in your local environment. 
```

### Install parliament and chromoseq Docker
```
# download parliament and chromoseq Docker
sudo docker pull docker.io/dnanexus/parliament2:latest # Parliament: https://github.com/dnanexus/parliament2
sudo docker pull  docker.io/zatawada/docker-basespace_chromoseq_v2:master # ChromoSeq: https://github.com/genome/docker-basespace_chromoseq
```

### Install CYTO-SV-ML Snakemake pipeline
```
git clone https://github.com/tzhang-nmdp/CYTO-SV-ML.git
cd CYTO-SV-ML
mamba env create -n py27 -f py27.yaml
mamba env create -n py39 -f py39.yaml
mamba env create -n cyto-sv-ml -f cyto-sv-ml.yaml
conda activate cyto-sv-ml
pip install --upgrade snakemake # in case that snakemake version is old
# note: python-graphviz might be conflicted with the pre-existing packages in your environment. If it happened, please remove it from requirement file and install it separately.
```

### Download hg38 reference genome 
```
cd CYTO-SV-ML/reference
mkdir CYTO-SV-ML/reference/hg38
# hg38 from Broad Institute Google Cloud
wget https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
mv Homo_sapiens_assembly38.fasta reference/hg38/hs38.fasta
mv Homo_sapiens_assembly38.fasta.fai reference/hg38/hs38.fasta.fai
# Please donwload VEP GRCh38 from ensembl (chromoseq pipeline)
curl -O https://ftp.ensembl.org/pub/release-90/variation/VEP/homo_sapiens_vep_90_GRCh38.tar.gz
tar xzf homo_sapiens_vep_90_GRCh38.tar.gz
cd CYTO-SV-ML/reference/homo_sapiens/90_GRCh38
```

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Pipeline Usage
============
### 1. Run CYTO-SV-ML Snakemake preprocess pipeline for each sample
Please change the config.yaml according to your own environment settings:                                                    
{your_work_dir} --> input/output dir <br/>
{cyto_sv_ml_dir} --> cyto_sv_ml dir and software subdir and reference subdir <br/>
{Sample_1} --> sample id <br/>
Please create input folder for chromseq pipeline --> {your_work_dir}/in/${sample}.cram <br/>
Please create input folder for parliament2 pipeline --> {your_work_dir}/in/${sample}.bam <br/>
Please create output folder --> {your_work_dir}/out/ <br/>
```
conda activate cyto-sv-ml
snakemake --core ${number_of_cores} -s cyto-sv-ml-preprocess.smk --use-conda --config sample=${sample} gender=${gender}
snakemake --core ${number_of_cores} -s cyto-sv-ml-preprocess.smk --use-conda --config sample=${sample} gender=${gender} --report ${your_work_dir}/out/${sample}/${sample}_report.html
```
![CYTO-SV-ML snakemake report](cyto-sv-ml_snakemake_runtime_report.png)

### 2. Run CYTO-SV-ML Snakemake modeling pipeline for whole cohort
```
snakemake --core ${number_of_cores} -s cyto-sv-ml-automl.smk --config cohort_name=${cohort_name}
snakemake --core ${number_of_cores} -s cyto-sv-ml-automl.smk --config cohort_name=${cohort_name} --report ${your_work_dir}/out/${cohort_name}_report.html
# please kindly check sv_feature_variable.dict for sv feature variable dictionary
```
![CYTO-SV-ML snakemake report](cyto-sv-ml_snakemake_report.png)

### 3. Run CYTO-SV-ML Snakemake interface pipeline for Web-portal application
```
snakemake --core ${number_of_cores} -s cyto-sv-ml-interface.smk  --config cohort_name=${cohort_name} kfolds=${kfolds}
# to run the docker image in the local machine and open user interface with "http://localhost:8000/"
sudo docker run -d -p 8000:80 cyto-sv-ml-app:${sample_all}
```
The analysis summary of 494 MDS cohort using CYTO-SV-ML pipeline 
### Online Web-portal
```
https://cyto-sv-ml.b12x.org/
```
![CYTO-SV-ML Rshiny web-portal](cyto-sv-ml_web-portal.png)

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Download SV related Resource
=============================

## SV database Download (The official websites contain detailed information)
```
gnomAD: https://gnomad.broadinstitute.org/downloads
1000g: https://www.internationalgenome.org/phase-3-structural-variant-dataset
CytoAtlas: https://github.com/genome/docker-basespace_chromoseq/blob/master/workflow_files/chromoseq_translocations.bedpe
COSMIC: https://cancer.sanger.ac.uk/cosmic/download
```

## SV Tool Download
```
SURVIVOR: https://github.com/fritzsedlazeck/SURVIVOR
SVTyper: https://github.com/hall-lab/svtyper
SeqComplex: https://github.com/caballero/SeqComplex
Komplexity: https://github.com/eclarke/komplexity
```

## ML Pipeline Download
```
AUTOML: https://github.com/mljar/mljar-supervised
```
