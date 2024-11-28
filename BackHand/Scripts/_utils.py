import os
import pathlib
import yaml
import sys

try:
    BASE_PATH = sys.argv[1]
except IndexError:
    BASE_PATH = "/home/labs/nyosef/yoavnah/CellRangerIDE/"
REFERENCE_PATH = "/shareDB/CellRanger/"

# Config files declaration
CONST_PATH = os.path.join(BASE_PATH,"const_files")
CONFIG_GENERAL_PATH = os.path.join(CONST_PATH, "config_general.yml")
with open(CONFIG_GENERAL_PATH,'r') as f:
    config = yaml.safe_load(f)

CONFIG_SELECTED_PIPELINE_PATH = os.path.join(CONST_PATH, "config_" + config['pipeline'] + "_pipeline.yml")
with open(CONFIG_SELECTED_PIPELINE_PATH,'r') as g:
    config_selected_pipeline = yaml.safe_load(g) 

# Project building folders
project_name = config['project_name']
PROJECT_PATH = os.path.join(BASE_PATH, "Projects")
PROJECT_NAME_PATH = os.path.join(PROJECT_PATH, project_name)
FILE_PATH = os.path.join(PROJECT_NAME_PATH, "File_Path")
QC_PATH = os.path.join(FILE_PATH, "10x_QC")
BCL_PATH = os.path.join(FILE_PATH, "BCL")
CSV_PATH = os.path.join(FILE_PATH, "CSV")
FASTQ_PATH = os.path.join(FILE_PATH, "Fastqs")
CUSTOM_REF_PATH = os.path.join(FILE_PATH, "Custom_Ref")
FASTQC_PATH = os.path.join(FASTQ_PATH, "fastqc")
H5_PATH = os.path.join(FILE_PATH, "h5")
H5ADS_PATH = os.path.join(FILE_PATH, "h5ads")
DEMULTIPLEXED_H5ADS_PATH = os.path.join(H5ADS_PATH, "demultiplexed")
BEFORE_DEMULTI_H5ADS_PATH = os.path.join(H5ADS_PATH, "not_demultiplexed")
OUTPUT_PATH = os.path.join(PROJECT_NAME_PATH, "OutputFiles/")

# Shell files location

mkfastq_path = os.path.join(OUTPUT_PATH, "mkfastq.sh")
make_multi_path = os.path.join(OUTPUT_PATH, "make_multi.sh")
make_count_path = os.path.join(OUTPUT_PATH, "make_count.sh")
make_h5_path = os.path.join(OUTPUT_PATH, "make_h5.sh")
make_h5ads_path = os.path.join(OUTPUT_PATH, "make_h5ads.sh")
demulti_path = os.path.join(OUTPUT_PATH, "demulti.sh")
make_reference_path = os.path.join(OUTPUT_PATH, "make_ref.sh")
basic_run_path = os.path.join(OUTPUT_PATH, "basic_run.sh")
run_file_path = os.path.join(OUTPUT_PATH, "run.sh")

# Aligners

DEFAULT_ALIGNER_PATH = "/apps/easybd/easybuild/amd/software/CellRanger/8.0.1/bin/cellranger"

# csv format executing file

MULTI_TEMPLATE_CSV_PATH = os.path.join(CONST_PATH, "multi_template.csv")
MULTI_CSV_PATH = os.path.join(CSV_PATH, "multi.csv")
samplesheet_path = os.path.join(CSV_PATH,"samplesheet.csv")
feature_reference_path = os.path.join(CSV_PATH,"feature_reference.csv")

# References
gex_reference_path = os.path.join(REFERENCE_PATH, "refdata-gex-mm10-2020-A")
vdj_reference_path = os.path.join(REFERENCE_PATH, "refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0")
custom_gtf_reference_path = os.path.join(CUSTOM_REF_PATH, "custom_ref.gtf")
custom_fasta_reference_path = os.path.join(CUSTOM_REF_PATH, "custom_ref.fasta")

# keys
GENE_SYMBOLS_KEY = "gene_symbols"
GEX_KEY = "Gene Expression"
GENE_IDS_KEY = "gene_ids"
FEATURE_TYPE_KEY = "feature_types"
SAMPLE_ID_KEY = "sample_id"
CELL_GROUP_KEY = "cell_group"
MODEL_KEY = "model"
ANTIBODY_KEY = "Antibody Capture"
CLASSIFICATION_KEY = "Classification"

# program runner

R_RUNNER = "/home/labs/nyosef/yoavnah/miniconda3/envs/yoav_env/envs/r_env/bin/Rscript"
