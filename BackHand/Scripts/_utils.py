import os
import pathlib
BASE_PATH = "/home/labs/nyosef/yoavnah/CellRangerIDE/"
import yaml

# Project building folders
CONST_PATH = os.path.join(BASE_PATH,"const_files")
CONFIG_PATH = os.path.join(CONST_PATH, "config.yml")
with open(CONFIG_PATH,'r') as f:
    config = yaml.safe_load(f)

project_name = config['project_name']
PROJECT_PATH = os.path.join(BASE_PATH, "Projects")
PROJECT_NAME_PATH = os.path.join(PROJECT_PATH, project_name)
FILE_PATH = os.path.join(PROJECT_NAME_PATH, "File_Path")
QC_PATH = os.path.join(FILE_PATH, "10x_QC")
BCL_PATH = os.path.join(FILE_PATH, "BCL")
CSV_PATH = os.path.join(FILE_PATH, "CSV")
FASTQ_PATH = os.path.join(FILE_PATH, "Fastqs")
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
basic_run_path = os.path.join(OUTPUT_PATH, "basic_run.sh")

# csv format executing file

PROGRAM_PATH = os.path.join(BASE_PATH, "program.csv")
MULTI_TEMPLATE_CSV_PATH = os.path.join(CONST_PATH, "multi_template.csv")
MULTI_CSV_PATH = os.path.join(CSV_PATH, "multi.csv")

# References
REFERENCE_PATH = os.path.join(CONST_PATH, "transcriptome")
gex_reference_path = os.path.join(REFERENCE_PATH, "refdata-gex-GRCh38-2020-A")
vdj_reference_path = os.path.join(REFERENCE_PATH, "refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0")

# keys
GENE_SYMBOLS_KEY = "gene_symbols"
GENE_IDS_KEY = "gene_ids"
FEATURE_TYPE_KEY = "feature_types"
SAMPLE_ID_KEY = "sample_id"
CELL_GROUP_KEY = "cell_group"
MODEL_KEY = "model"
ANTIBODY_KEY = "Antibody Capture"
CLASSIFICATION_KEY = "Classification"

# experimental design
MULTIPLEX_BARCODE_TO_SAMPLE_ID = {
    "AT_Epi_Naive2": {
        "AT_Epi": "AT_Epi_Naive2",
    },
    "AT_Epi_Naive3": {
        "AT_Epi": "AT_Epi_Naive3",
    },
    "AT_Epi_Act147": {
        "AT_Epi": "AT_Epi_Act147",
    },
    "AT_Epi_Act150": {
        "AT_Epi": "AT_Epi_Act150",
    }
}
SAMPLE_ID_TO_EXPECTED_BARCODES = {}
for barcode, barcode_dict in MULTIPLEX_BARCODE_TO_SAMPLE_ID.items():
    for sample_id, value in barcode_dict.items():
        if sample_id not in SAMPLE_ID_TO_EXPECTED_BARCODES:
            barcodes = []
        else:
            barcodes = SAMPLE_ID_TO_EXPECTED_BARCODES[sample_id]
        if value is not None and barcode not in barcodes:
            barcodes.append(barcode)
        SAMPLE_ID_TO_EXPECTED_BARCODES[sample_id] = barcodes
