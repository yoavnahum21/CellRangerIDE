import os
import pathlib
import yaml
import sys
<<<<<<< HEAD

try:
    BASE_PATH =                             sys.argv[1]
except IndexError:
    BASE_PATH =                             __file__
    for _ in range(3):
        BASE_PATH =                         os.path.dirname(BASE_PATH)
=======
import argparse
import ast

parser = argparse.ArgumentParser()

parser.add_argument("--cellrangeride-path", "--cellrangeride_path", type=str, help="Project name")
parser.add_argument("--pipeline", "--pipeline", type=str, help="Pipeline type")
parser.add_argument("--project-name", "--project_name", type=str, help="Project name")
parser.add_argument("--id", "--id", type=str, help="Project ID")
parser.add_argument("--running-machine", "--running_machine", type=str, help="Running machine")
parser.add_argument("--aligner-software-path", "--aligner_software_path", type=str, help="Aligner software path")
parser.add_argument("--reference-genome", "--reference_genome", type=str, help="Reference genome")      
parser.add_argument("--cellbender-path", "--cellbender_path", type=str, help="Cellbender path")
parser.add_argument("--jobs-number", "--jobs_number", type=int, help="number of jobs")
parser.add_argument("--memory-size", "--memory_size", type=int, help="memory size") 
parser.add_argument("--cpu-cores", "--cpu_cores", type=int, help="CPU cores number")
parser.add_argument("--queue-name", "--queue_name", type=str, help="Queue name")
parser.add_argument("--fastq-folder-path", "--fastq_folder_path", type=lambda s: s.split(","), help="Fastq path")
parser.add_argument("--include-introns", "--include_introns", type=bool, help="Include introns in the analysis")
parser.add_argument("--create-bam", "--create_bam", type=bool, help="Create BAM file")
parser.add_argument("--incpm-link", "--incpm_link", type=lambda s: s.split(","), help="incpm link")
parser.add_argument("--expected_cells", "--expected_cells", type=int, help="Expected cells number")
parser.add_argument("--R1-length", "--R1_length", type=int, help="R1 length")
parser.add_argument("--R2-length", "--R2_length", type=int, help="R2 length")
parser.add_argument("--R1-vdj-length", "--R1_vdj_length", type=int, help="R1 VDJ length")
parser.add_argument("--R2-vdj-length", "--R2_vdj_length", type=int, help="R2 VDJ length")
parser.add_argument("--sample-names", "--sample_names", type=lambda s: s.split(","), help="Sample names")
parser.add_argument("--fastq-folders-names", "--fastq_folders_names", type=lambda s: s.split(","), help="Fastq folders names")
parser.add_argument("--multiplexing-method", "--multiplexing_method", type=str, help="Multiplexing method")
parser.add_argument("--feature-reference-csv", "--feature_reference_csv", type=str, help="Feature reference CSV")
parser.add_argument("--feature-types", "--feature_types", type=lambda s: s.split(","), help="Feature types")
parser.add_argument("--hto-ids", "--hto_ids", type=lambda s: s.split(","), help="HTO IDs")
parser.add_argument("--hto-names", "--hto_names", type=lambda s: s.split(","), help="HTO names")
parser.add_argument("--hto-reads", "--hto_reads", type=lambda s: s.split(","), help="HTO reads")
parser.add_argument("--hto-pattern", "--hto_pattern", type=lambda s: s.split(","), help="HTO pattern")
parser.add_argument("--hto-sequence", "--hto_sequence", type=lambda s: s.split(","), help="HTO sequence")
parser.add_argument("--hto-feature-type", "--hto_feature_type", type=lambda s: s.split(","), help="HTO feature type")
parser.add_argument("--sample-id-cmo", "--sample_id_cmo", type=lambda s: s.split(","), help="Sample ID CMO")
parser.add_argument("--cmo-barcode-csv", "--cmo_barcode_csv", type=str, help="CMO barcode CSV")
parser.add_argument("--cmo-id", "--cmo_id", type=lambda s: s.split(","), help="CMO ID")
parser.add_argument("--probe-set-path", "--probe_set_path", type=str, help="Probe set path")
parser.add_argument("--probe-barcode-csv", "--probe_barcode_csv", type=str, help="Probe barcode CSV")
parser.add_argument("--sample-id-probes", "--sample_id_probes", type=lambda s: s.split(","), help="Sample ID probes")
parser.add_argument("--probe-barcode-ids", "--probe_barcode_ids", type=lambda s: s.split(","), help="Probe barcode IDs")
parser.add_argument("--probe-description", "--probe_description", type=lambda s: s.split(","), help="Probe description")
parser.add_argument("--demultiplex-method", "--demultiplex_method", type=str, help="Demultiplex method")
parser.add_argument("--priors-probabilities", "--priors_probabilities", type=lambda s: s.split(","), help="Priors probabilities")
parser.add_argument("--adata-path", "--adata_path", type=str, help="adata path")
parser.add_argument("--chosen-pipeline", "--chosen_pipeline", type=str, help="Chosen pipeline")
args = parser.parse_args() 

# if args.incpm_link:
#     try:
#         args.incpm_link = ast.literal_eval(args.incpm_link)
#     except Exception:
#         print("Could not parse incpm_link argument")
#         print(Exception)


# Base path determination

BASE_PATH =                                args.cellrangeride_path if args.cellrangeride_path else str(pathlib.Path(__file__).parent.parent.parent.resolve())
>>>>>>> argparser
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
<<<<<<< HEAD
project_name =                              config['project_name']
=======
project_name =                              args.project_name if args.project_name else config['project_name']
id =                                        args.id if args.id else config['id']
>>>>>>> argparser
PROJECT_PATH =                              os.path.join(BASE_PATH, "Projects")
PROJECT_NAME_PATH =                         os.path.join(PROJECT_PATH, project_name)
FILE_PATH =                                 os.path.join(PROJECT_NAME_PATH, "File_Path")
QC_PATH =                                   os.path.join(FILE_PATH, "10x_QC")
BCL_PATH =                                  os.path.join(FILE_PATH, "BCL")
CSV_PATH =                                  os.path.join(FILE_PATH, "CSV")
FASTQ_PATH =                                os.path.join(FILE_PATH, "Fastqs")
CUSTOM_REF_PATH =                           os.path.join(FILE_PATH, "Custom_Ref")
<<<<<<< HEAD
CELLBENDER_PATH =                           os.path.join(FILE_PATH, "Cellbender")
=======
>>>>>>> argparser
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
<<<<<<< HEAD
multi_log =                                 os.path.join(LOGS_PATH, "multi_" + config['id'] + ".log")
multi_error_log =                           os.path.join(LOGS_PATH, "multi_" + config['id'] + ".err.log")
count_log =                                 os.path.join(LOGS_PATH, "count_" + config['id'] + ".log")
count_error_log =                           os.path.join(LOGS_PATH, "count_" + config['id'] + ".err.log")
demultiplex_log =                           os.path.join(LOGS_PATH, "deMULTIplex_" + config['id'] + ".log")
demultiplex_error_log =                     os.path.join(LOGS_PATH, "deMULTIplex_" + config['id'] + ".err.log")
qc_log =                                    os.path.join(LOGS_PATH, "qc_" + config['id'] + ".log")
qc_error_log =                              os.path.join(LOGS_PATH, "qc_" + config['id'] + ".log")
=======
multi_log =                                 os.path.join(LOGS_PATH, "multi_" + id + ".log")
multi_error_log =                           os.path.join(LOGS_PATH, "multi_" + id + ".err.log")
demultiplex_log =                           os.path.join(LOGS_PATH, "deMULTIplex_" + id + ".log")
demultiplex_error_log =                     os.path.join(LOGS_PATH, "deMULTIplex_" + id + ".err.log")
qc_log =                                    os.path.join(LOGS_PATH, "qc_" + id + ".log")
qc_error_log =                              os.path.join(LOGS_PATH, "qc_" + id + ".log")
>>>>>>> argparser

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

<<<<<<< HEAD
# Scripts
SCRIPTS_PATH =                              os.path.join(BASE_PATH, "BackHand/Scripts")
# Cellbender runner files path          
CELLBENDER_PYTHON_RUNNER =                  os.path.join(SCRIPTS_PATH, "cellbenderPipe.py")
CELLBENDER_SHELL_RUNNER =                   os.path.join(SCRIPTS_PATH, "run_cellbender.sh")
=======

# Cellbender runner files path
CELLBENDER_PYTHON_RUNNER =                  "/home/projects/nyosef/yoavnah/CellRangerIDE/BackHand/Scripts/cellbender.py"
CELLBENDER_SHELL_RUNNER  =                  "/home/projects/nyosef/yoavnah/CellRangerIDE/BackHand/Scripts/run_cellbender.sh"

# Supported pipelines
SUPPORTED_PIPELINES =                       ['mkfastq', 'multi', 'count', 'demulti', 'mkref', 'cellbender', 'flex', 'velocyto', 'QC']
>>>>>>> argparser
