import os
import pathlib
import yaml
import sys

try:
    BASE_PATH =                             sys.argv[1]
except IndexError:
    BASE_PATH =                             __file__
    for _ in range(3):
        BASE_PATH =                         os.path.dirname(BASE_PATH)
REFERENCE_PATH =                            "/shareDB/CellRanger/"

# Config files declaration
CONST_PATH =                                os.path.join(BASE_PATH,"const_files")
CONFIG_GENERAL_PATH =                       os.path.join(CONST_PATH, "config_general.yml")
with open(CONFIG_GENERAL_PATH,'r') as f:
    config =                                yaml.safe_load(f)

CONFIG_SELECTED_PIPELINE_PATH =             os.path.join(CONST_PATH, "config_" + config['pipeline'] + "_pipeline.yml")
with open(CONFIG_SELECTED_PIPELINE_PATH,'r') as g:
    config_selected_pipeline =              yaml.safe_load(g) 

# Project building folders
project_name =                              config['project_name']
PROJECT_PATH =                              os.path.join(BASE_PATH, "Projects")
PROJECT_NAME_PATH =                         os.path.join(PROJECT_PATH, project_name)
FILE_PATH =                                 os.path.join(PROJECT_NAME_PATH, "File_Path")
QC_PATH =                                   os.path.join(FILE_PATH, "10x_QC")
BCL_PATH =                                  os.path.join(FILE_PATH, "BCL")
CSV_PATH =                                  os.path.join(FILE_PATH, "CSV")
FASTQ_PATH =                                os.path.join(FILE_PATH, "Fastqs")
CUSTOM_REF_PATH =                           os.path.join(FILE_PATH, "Custom_Ref")
CELLBENDER_PATH =                           os.path.join(FILE_PATH, "Cellbender")
FASTQC_PATH =                               os.path.join(FASTQ_PATH, "fastqc")
H5_PATH =                                   os.path.join(FILE_PATH, "h5")
H5ADS_PATH =                                os.path.join(FILE_PATH, "h5ads")
DEMULTIPLEXED_H5ADS_PATH =                  os.path.join(H5ADS_PATH, "demultiplexed")
BEFORE_DEMULTI_H5ADS_PATH =                 os.path.join(H5ADS_PATH, "not_demultiplexed")
LOGS_PATH =                                 os.path.join(FILE_PATH, "Logs")
INCPM_PATH =                                os.path.join(FILE_PATH, "INCP")
FASTQ_INCPM_PATH =                          os.path.join(INCPM_PATH, "fastqs")
SAMPLE_SHEET_INCPM_PATH =                   os.path.join(INCPM_PATH, "sample_sheet")
CSV_INCPM_PATH =                            os.path.join(INCPM_PATH, "CSV")
OUTPUT_PATH =                               os.path.join(PROJECT_NAME_PATH, "OutputFiles/")

# Shell files location
mkfastq_path =                              os.path.join(OUTPUT_PATH, "mkfastq.sh")
make_multi_path =                           os.path.join(OUTPUT_PATH, "make_multi.sh")
make_count_path =                           os.path.join(OUTPUT_PATH, "make_count.sh")
make_qc_path =                              os.path.join(OUTPUT_PATH, "make_qc.sh")
make_h5_path =                              os.path.join(OUTPUT_PATH, "make_h5.sh")
make_h5ads_path =                           os.path.join(OUTPUT_PATH, "make_h5ads.sh")
demulti_path =                              os.path.join(OUTPUT_PATH, "demulti.sh")
make_reference_path =                       os.path.join(OUTPUT_PATH, "make_ref.sh")
basic_run_path =                            os.path.join(OUTPUT_PATH, "basic_run.sh")
run_file_path =                             os.path.join(OUTPUT_PATH, "run_on_wexac.sh")

# Aligners
DEFAULT_ALIGNER_PATH =                      "/apps/easybd/easybuild/amd/software/CellRanger/9.0.0/bin/cellranger"

# csv format executing file
MULTI_TEMPLATE_CSV_PATH =                   os.path.join(CONST_PATH, "multi_template.csv")
MULTI_CSV_PATH =                            os.path.join(CSV_PATH, "multi.csv")
MULTI_FLEX_CSV_PATH =                       os.path.join(CSV_PATH, "multi_flex.csv")
samplesheet_path =                          os.path.join(CSV_PATH,"samplesheet.csv")
feature_reference_path =                    os.path.join(CSV_PATH,"feature_reference.csv")

# References
gex_reference_path =                        os.path.join(REFERENCE_PATH, "refdata-gex-mm10-2020-A")
vdj_reference_path =                        os.path.join(REFERENCE_PATH, "refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0")
custom_gtf_reference_path =                 os.path.join(CUSTOM_REF_PATH, "custom_ref.gtf")
custom_fasta_reference_path =               os.path.join(CUSTOM_REF_PATH, "custom_ref.fasta")

# Logs
multi_log =                                 os.path.join(LOGS_PATH, "multi_" + config['id'] + ".log")
multi_error_log =                           os.path.join(LOGS_PATH, "multi_" + config['id'] + ".err.log")
count_log =                                 os.path.join(LOGS_PATH, "count_" + config['id'] + ".log")
count_error_log =                           os.path.join(LOGS_PATH, "count_" + config['id'] + ".err.log")
demultiplex_log =                           os.path.join(LOGS_PATH, "deMULTIplex_" + config['id'] + ".log")
demultiplex_error_log =                     os.path.join(LOGS_PATH, "deMULTIplex_" + config['id'] + ".err.log")
qc_log =                                    os.path.join(LOGS_PATH, "qc_" + config['id'] + ".log")
qc_error_log =                              os.path.join(LOGS_PATH, "qc_" + config['id'] + ".log")

# Keys
GENE_SYMBOLS_KEY =                          "gene_symbols"
GEX_KEY =                                   "Gene Expression"
GENE_IDS_KEY =                              "gene_ids"
FEATURE_TYPE_KEY =                          "feature_types"
SAMPLE_ID_KEY =                             "sample_id"
CELL_GROUP_KEY =                            "cell_group"
MODEL_KEY =                                 "model"
ANTIBODY_KEY =                              "Antibody Capture"
CLASSIFICATION_KEY =                        "Classification"

# Program runner

# Scripts
SCRIPTS_PATH =                              os.path.join(BASE_PATH, "BackHand/Scripts")
# Cellbender runner files path          
CELLBENDER_PYTHON_RUNNER =                  os.path.join(SCRIPTS_PATH, "cellbenderPipe.py")
CELLBENDER_SHELL_RUNNER =                   os.path.join(SCRIPTS_PATH, "run_cellbender.sh")