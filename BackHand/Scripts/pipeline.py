from _utils import * 
import os
import sys
import shutil
import subprocess
import demultiplex as demultiplex
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import requests
import subprocess

MEM_PER_TASK = 70
CPUS_PER_TASK = 4
MAX_JOBS = 400
EXPECTED_CELLS = 10000

class pipeline:

    def __init__(self) -> None:

        if not os.path.exists(PROJECT_NAME_PATH):
            os.makedirs(PROJECT_NAME_PATH)
            os.makedirs(FILE_PATH)
            os.makedirs(QC_PATH)
            os.makedirs(BCL_PATH)
            os.makedirs(CSV_PATH)
            os.makedirs(FASTQ_PATH)
            os.makedirs(CUSTOM_REF_PATH)
            os.makedirs(FASTQC_PATH)
            os.makedirs(H5_PATH)
            os.makedirs(H5ADS_PATH)
            os.makedirs(BEFORE_DEMULTI_H5ADS_PATH)
            os.makedirs(DEMULTIPLEXED_H5ADS_PATH)
            os.makedirs(INCPM_PATH)
            os.makedirs(FASTQ_INCPM_PATH)
            os.makedirs(SAMPLE_SHEET_INCPM_PATH)
            os.makedirs(CSV_INCPM_PATH)
            os.makedirs(OUTPUT_PATH)

        else:
            print("Welcome to {}".format(PROJECT_NAME_PATH))
        
        self.pipeline                   = config['pipeline']

        self.donor                      = config['donor']

        self.seq_run                    = config['seq_run']

        self.id                         = config['id']

        self.output_destination         = config['output_destination']
        if self.output_destination      == None:
            self.output_destinatxion    = OUTPUT_PATH
        
        self.running_machine            = config['running_machine']
        if self.running_machine         == 'Default':
            self.running_machine        = "Wexac"

        self.aligner_software_path      = config['aligner_software_path']
        if self.aligner_software_path   == None:
            self.aligner_software_path  = DEFAULT_ALIGNER_PATH

        if self.pipeline                 == "mkfastq":
            self.BCL_path                = config_selected_pipeline['BCL_path']
            self.Lanes                   = config_selected_pipeline['Lanes']
            self.sample_name             = config_selected_pipeline['sample_name']
            self.index                   = config_selected_pipeline['index']
            # one day the rest of the options will be implemented as well

        elif self.pipeline                   == "count":
            self.aligner_ref_genome_path     = config_selected_pipeline['alignment_ref_genome_file']
            if self.aligner_ref_genome_path  == None:
                self.aligner_ref_genome_path = gex_reference_path
            self.aligner_ref_vdj_path        = config_selected_pipeline['alignment_ref_vdj_file']
            self.fastq_path                  = config_selected_pipeline['BCL_path']
            self.sample_name                 = config_selected_pipeline['sample_name']
            # one day the rest of the options will be implemented as well

        elif self.pipeline                   == "multi":
            self.aligner_ref_genome_path     = config_selected_pipeline['alignment_ref_genome_file']
            self.create_bam                  = config_selected_pipeline['create_bam']
            self.aligner_ref_genome_path     = config_selected_pipeline['alignment_ref_genome_file']
            self.expected_cells               = config_selected_pipeline['expected_cell']
            self.include_introns             = config_selected_pipeline['include_introns']
            if self.aligner_ref_genome_path  == None:
                self.aligner_ref_genome_path = gex_reference_path
            self.aligner_ref_vdj_path        = config_selected_pipeline['alignment_ref_vdj_file']
            self.sample_sheet                = config_selected_pipeline['INCPM_link']
            if self.sample_sheet             == None:
                self.fastq_path                     = config_selected_pipeline['fastq_path']
                self.multiplexing_method            = config_selected_pipeline['multiplexing_method']
                self.fastq_folders_name             = config_selected_pipeline['fastq_folders_name']
                self.lanes_used                     = config_selected_pipeline['lanes_used']
                if self.lanes_used[0]               == 'Default':
                    self.lanes_used                 =  ['any'] * len(self.lanes_used)
                self.feature_types                  = config_selected_pipeline['feature_types']
                if self.multiplexing_method         == "feature_barcode":
                    self.feature_reference_csv      = config_selected_pipeline['feature_reference_csv']
                    if self.feature_reference_csv   == 'Default':
                        self.feature_reference_csv  =   feature_reference_path
                        self.hto_id                 = config_selected_pipeline['hto_id']
                        self.hto_names              = config_selected_pipeline['hto_names']
                        if self.hto_names[0]        == 'Default':
                            self.hto_names          = self.hto_id
                        self.hto_pattern            = config_selected_pipeline['hto_pattern']
                        if self.hto_pattern[0]      == 'Default':
                            self.hto_pattern        = ['5PNNNNNNNNNN(BC)'] * len(self.hto_pattern)
                        self.hto_read               = config_selected_pipeline['hto_read']
                        if self.hto_read[0]         == 'Default':
                            self.hto_read           = ['R2'] * len(self.hto_read)
                        self.hto_sequence           = config_selected_pipeline['hto_sequence']
                        self.HTO_feature_type       = config_selected_pipeline['HTO_feature_type']
                        if self.HTO_feature_type[0] == 'Default':
                            self.HTO_feature_type   = ['Antibody Capture'] * len(self.HTO_feature_type)

                if self.multiplexing_method         == "cmo_barcode":
                    self.cmo_barcode_csv            = config_selected_pipeline['cmo_barcode_csv']
                    if self.cmo_barcode_csv         == 'Default':
                        self.sample_id_cmo              = config_selected_pipeline['sample_id_cmo']
                        self.cmo_id                 = config_selected_pipeline['cmo_id']

            else:
                incpm_path                          = config_selected_pipeline['INCPM_directory']
                if incpm_path                       == 'Default':
                    incpm_path                      = INCPM_PATH
                incpm_link                          = config_selected_pipeline['INCPM_link']
                self.sample_names                   = []
                self.lanes_used                     = []
                self.fastq_folders_name             = []
                self.R1_gex_len                     = config_selected_pipeline['R1_gex_length']
                self.R2_gex_len                     = config_selected_pipeline['R2_gex_length']
                self.R1_vdj_len                     = config_selected_pipeline['R1_vdj_length']
                self.R2_vdj_len                     = config_selected_pipeline['R2_vdj_length']
                self.multiplexing_method            = config_selected_pipeline['multiplexing_method']
                self.feature_reference_csv          = config_selected_pipeline['feature_reference_csv']
                self.hto_id                         = config_selected_pipeline['hto_id']
                self.hto_names                      = config_selected_pipeline['hto_names']
                if self.hto_names[0]                == 'Default':
                    self.hto_names                  = self.hto_id
                self.hto_pattern                    = config_selected_pipeline['hto_pattern']
                if self.hto_pattern[0]              == 'Default':
                    self.hto_pattern                = ['5PNNNNNNNNNN(BC)'] * len(self.hto_pattern)
                self.hto_read                       = config_selected_pipeline['hto_read']
                if self.hto_read[0]                 == 'Default':
                    self.hto_read                   = ['R2'] * len(self.hto_read)
                self.hto_sequence                   = config_selected_pipeline['hto_sequence']
                self.HTO_feature_type               = config_selected_pipeline['HTO_feature_type']
                if self.HTO_feature_type[0]         == 'Default':
                    self.HTO_feature_type           = ['Antibody Capture'] * len(self.HTO_feature_type)
                self.feature_types                  = []
                self.fastq_path                     = []
                feautre_type_list                   = ["gene expression","Antibody Capture","CRISPR Guide Capture","Multiplexing Capture","VDJ-B","VDJ-T","VDJ-T-GD","Antigen Capture #5' Antigen Capture only"]    
                for i in range(len(incpm_link)):
                    self.sample_sheet_path              = os.path.join(SAMPLE_SHEET_INCPM_PATH, "SampleSheet" + str(i+1) + ".csv")
                    copy_process_samplesheet            = os.path.join(SAMPLE_SHEET_INCPM_PATH, "copy_samplesheet" + str(i+1) + ".sh")
                    with open(copy_process_samplesheet, "w+") as f:
                        f.write(
f"wget --no-check-certificate {incpm_link[i]}SampleSheet.csv --output-document={self.sample_sheet_path}"
                    )
                    os.system(copy_process_samplesheet)

                    #### patch if csv is not valid ###
                    with open(self.sample_sheet_path, "r") as f:
                        tmp_csv_file = f.readlines()
                    
                    comma_max_count = 0
                    for row in tmp_csv_file:
                        comma_count = row.count(",")
                        if comma_count > comma_max_count:
                            comma_max_count = comma_count

                    phase = 0
                    for index, row in enumerate(tmp_csv_file):
                        row =  row[:-1] + ","*(comma_max_count-row.count(",")) + "\n"
                        if row.split(",")[1]:
                            phase += 1
                            continue
                        tmp_csv_file[index - phase] = row
                    
                    with open(self.sample_sheet_path, "w") as f:
                        f.writelines(tmp_csv_file)

                    #### end of patch ###
                    self.sample_sheet                   = pd.read_csv(self.sample_sheet_path)
                    samplesheet_relevant_index          = self.sample_sheet[self.sample_sheet.isin(['Sample_ID']).any(axis=1)].index[0]
                    df                                  = self.sample_sheet.iloc[samplesheet_relevant_index+1:].reset_index(drop=True)
                    self.lanes_used                     += df.iloc[:,0].tolist()
                    self.sample_names                   += df.iloc[:,1].tolist()
                    self.fastq_folders_name             += df.iloc[:,1].tolist()


                answer = str(input("Do you wish to use all libraries ??y/n\n"))
                while (answer == 'n' or answer == 'N'):
                    print("which library would you like to remove")               
                    for index, library in enumerate(self.fastq_folders_name):
                        print(str(index + 1)+". " + library)
                    rm_elem = int(input())
                    self.lanes_used.pop(rm_elem-1)
                    self.sample_names.pop(rm_elem-1)
                    self.fastq_folders_name.pop(rm_elem-1)
                    answer = str(input("what about now, do you wish to use all libraries??y/n\n"))
                print("keep going on with the next list")

                for index, library in enumerate(self.fastq_folders_name):
                    print(str(index + 1)+". " + library)
                
                for sample in self.sample_names:
                    self.fastq_path.append(f"{FASTQ_INCPM_PATH}/{sample}/")
                    print(f"which feautre type would you like to use for the sample: {sample}")
                    type = int(input("1.Gene Expression\n2.Antibody Capture\n3.CRISPR Guide Capture\n4.Multiplexing Capture\n5.VDJ-B\n6.VDJ-T\n7.VDJ-T-GD\n8.Antigen Capture #5' Antigen Capture only\n"))
                    self.feature_types.append(feautre_type_list[type-1])
                copy_process_fastqs                 = os.path.join(FASTQ_INCPM_PATH, "copy_fastq.sh")
                if not os.listdir(FASTQ_INCPM_PATH):
                    for link in incpm_link:
                        for sample in self.sample_names:
                            with open(copy_process_fastqs, "a") as f:
                                f.write(
    f'wget -r -np -nd -A "*.fastq.gz" --no-check-certificate {link}{sample}/ -P {FASTQ_INCPM_PATH}/{sample}/\n'
                                ) 
                    os.system(copy_process_fastqs)

            
            # one day the rest of the options will be implemented as well
        elif self.pipeline                   == "mkref":
            self.original_fasta_file_path    = config_selected_pipeline['original_fasta_file_path']
            self.customized_fasta_file_path  = config_selected_pipeline['customized_fasta_file_path']
            self.fasta_of_new_gene           = config_selected_pipeline['fasta_of_new_gene']
            self.original_gtf_file_path      = config_selected_pipeline['original_gtf_file_path']
            self.customized_gtf_file_path    = config_selected_pipeline['customized_gtf_file_path']
            self.gtf_of_new_gene             = config_selected_pipeline['gtf_of_new_gene']
            self.filter_gtf_attribute        = config_selected_pipeline['filter_gtf_attribute']
            self.genome_name                 = config_selected_pipeline['genome_name']


        elif self.pipeline                      == "demulti":
            self.demultiplex_method             = config_selected_pipeline['demultiplex_method']
            self.adata_path                     = config_selected_pipeline['adata_path']
            self.sample_names                   = config_selected_pipeline['sample_names']
            self.hto_list_per_sample            = config_selected_pipeline['hto_list_per_sample']
            self.sample_id_to_expected_barcodes = {}
            for index, sample in enumerate(self.sample_names):
                self.sample_id_to_expected_barcodes.update({sample: self.hto_list_per_sample[index]})


        elif self.pipeline                   == "vdj":
            pass


        elif self.pipeline                   == "aggr":
            pass

        # Runtime parameters

        self.cpus_number                = config['cpus_number']
        self.jobs_number                = config['jobs_number']
        self.total_memory_used          = config['total_memory_used']
        self.queue                      = config['queue']
    
    def __str__(self) -> str:
        pass


    def mkfastq(self) -> None:

        '''
Required:
    --run=PATH          Path of Illumina BCL run folder.

Optional:
# Sample Sheet
    --id=NAME           Name of the folder created by mkfastq. If not supplied,
                            will default to the name of the flowcell referred to
                            by the --run argument.
    --csv=PATH
    --samplesheet=PATH
    --sample-sheet=PATH
                        Path to the sample sheet. The sample sheet can either be
                            a simple CSV with lane, sample and index columns, or
                            an Illumina Experiment Manager-compatible sample
                            sheet. Sample sheet indexes can refer to 10x sample
                            index set names (e.g., SI-GA-A12).
    --simple-csv=PATH   Deprecated. Same meaning as --csv.
    --force-single-index
                        If 10x-supplied i7/i5 paired indices are specified,
                            but the flowcell was run with only one sample
                            index, allow the demultiplex to proceed using
                                the i7 half of the sample index pair.
    --filter-single-index
                        Only demultiplex samples identified
                            by an i7-only sample index, ignoring dual-indexed
                            samples.  Dual-indexed samples will not be
                            demultiplexed.
    --filter-dual-index
                        Only demultiplex samples identified
                          by i7/i5 dual-indices (e.g., SI-TT-A6), ignoring single-
                          index samples.  Single-index samples will not be 
                          demultiplexed.
    --rc-i2-override=BOOL
                        Indicates if the bases in the I2 read are emitted as 
                          reverse complement by the sequencing workflow.
                          Set to 'true' for the Reverse Complement Workflow
                          (Workflow B)/ NovaSeq Reagent Kit v1.5 or greater.
                          Set to 'false' for the Forward Strand Workflow
                          (Workflow A) / older NovaSeq Reagent Kits.
                          NOTE: this parameter is autodetected 
                          and should only be passed in special circumstances.
    
# bcl2fastq Pass-Through
    --lanes=NUMS        Comma-delimited series of lanes to demultiplex. Shortcut
                            for the --tiles argument.
    --use-bases-mask=MASK
                        Same as bcl2fastq; override the read lengths as
                            specified in RunInfo.xml. See Illumina bcl2fastq
                            documentation for more information.
    --delete-undetermined
                        Delete the Undetermined FASTQ files left by bcl2fastq
                            Useful if your sample sheet is only expected to
                            match a subset of the flowcell.
    --output-dir=PATH   Same as in bcl2fastq. Folder where FASTQs, reports and
                            stats will be generated.
    --project=NAME      Custom project name, to override the samplesheet or to
                            use in conjunction with the --csv argument.

# Martian Runtime
    --jobmode=MODE      Job manager to use. Valid options: local (default), sge,
                            lsf, or a .template file
    --localcores=NUM    Set max cores the pipeline may request at one time. Only
                            applies to local jobs.
    --localmem=NUM      Set max GB the pipeline may request at one time. Only
                            applies to local jobs.
    --localvmem=NUM     Set max virtual address space in GB for the pipeline.
                            Only applies to local jobs.
    --mempercore=NUM    Reserve enough threads for each job to ensure enough
                        memory will be available, assuming each core on your
                        cluster has at least this much memory available. Only
                            applies in cluster jobmodes.
    --maxjobs=NUM       Set max jobs submitted to cluster at one time. Only
                            applies in cluster jobmodes.
    --jobinterval=NUM   Set delay between submitting jobs to cluster, in ms.
                            Only applies in cluster jobmodes.
    --overrides=PATH    The path to a JSON file that specifies stage-level
                            overrides for cores and memory. Finer-grained
                            than --localcores, --mempercore and --localmem.
                            Consult the 10x support website for an example
                            override file.

    --uiport=PORT       Serve web UI at http://localhost:PORT
    --disable-ui        Do not serve the UI.
    --noexit            Keep web UI running after pipestance completes or fails.
    --nopreflight       Skip preflight checks.
        '''

        with open(mkfastq_path, "w") as f:

    
            f.write(
            f"""#!/usr/bin/bash

    #SBATCH --job-name={''}
    #SBATCH --output={''}
    #SBATCH --partition={''}
    #SBATCH --mem={''}G
    #SBATCH --cpus-per-task={''}

    cellranger mkfastq --run={''} --csv={''} --output-dir={''} \
        --maxjobs={MAX_JOBS} \
        --localcores={CPUS_PER_TASK} \
        --localmem={MEM_PER_TASK}
            """
            )



    def count(self):

        '''
        --id <ID>                                     A unique run id and output folder name [a-zA-Z0-9_-]+
        --description <TEXT>                          Sample description to embed in output files [default: ]
        --transcriptome <PATH>                        Path of folder containing 10x-compatible transcriptome reference
        --fastqs <PATH>                               Path to input FASTQ data
        --project <TEXT>                              Name of the project folder within a mkfastq or bcl2fastq-generated folder from which to pick FASTQs
        --sample <PREFIX>                             Prefix of the filenames of FASTQs to select
        --lanes <NUMS>                                Only use FASTQs from selected lanes
        --libraries <CSV>                             CSV file declaring input library data sources
        --feature-ref <CSV>                           Feature reference CSV file, declaring Feature Barcode constructs and associated barcodes
        --target-panel <CSV>                          The target panel CSV file declaring the target panel used, if any. Default analysis will exclude intronic mapped reads, which is the recommended mode for targeted assay. Use include-introns=true to
                                                      include intronic mapped reads in analysis
        --expect-cells <NUM>                          Expected number of recovered cells, used as input to cell calling algorithm
        --force-cells <NUM>                           Force pipeline to use this number of cells, bypassing cell calling algorithm. [MINIMUM: 10]
        --no-bam                                      Set --no-bam to not generate the BAM file. This will reduce the total computation time for the pipestance and the size of the output directory. If unsure, we recommend not to use this option. BAM file
                                                      could be useful for troubleshooting and downstream analysis
        --nosecondary                                 Disable secondary analysis, e.g. clustering. Optional
        --r1-length <NUM>                             Hard trim the input Read 1 to this length before analysis
        --r2-length <NUM>                             Hard trim the input Read 2 to this length before analysis
        --include-introns <true|false>                Include intronic reads in count (default=true unless --target-panel is specified in which case default=false)
        --chemistry <CHEM>                            Assay configuration. NOTE: by default the assay configuration is detected automatically, which is the recommended mode. You usually will not need to specify a chemistry. Options are: 'auto' for
                                                      autodetection, 'threeprime' for Single Cell 3', 'fiveprime' for  Single Cell 5', 'SC3Pv1' or 'SC3Pv2' or 'SC3Pv3' for Single Cell 3' v1/v2/v3, 'SC3Pv3LT' for Single Cell 3' v3 LT, 'SC3Pv3HT' for Single
                                                      Cell 3' v3 HT, 'SC5P-PE' or 'SC5P-R2' for Single Cell 5', paired-end/R2-only, 'SC-FB' for Single Cell Antibody-only 3' v2 or 5'. To analyze the GEX portion of multiome data, chemistry must be set to
                                                      'ARC-v1'; 'ARC-v1' chemistry cannot be autodetected [default: auto]
        --no-libraries                                Proceed with processing using a --feature-ref but no Feature Barcode libraries specified with the 'libraries' flag
        --check-library-compatibility <true|false>    Whether to check for barcode compatibility between libraries. [default: true]
        --no-target-umi-filter                        Turn off the target UMI filtering subpipeline. Only applies when --target-panel is used
        --dry                                         Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
        --jobmode <MODE>                              Job manager to use. Valid options: local (default), sge, lsf, slurm or path to a .template file. Search for help on "Cluster Mode" at support.10xgenomics.com for more details on configuring the
                                                      pipeline to use a compute cluster [default: local]
        --localcores <NUM>                            Set max cores the pipeline may request at one time. Only applies to local jobs
        --localmem <NUM>                              Set max GB the pipeline may request at one time. Only applies to local jobs
        --localvmem <NUM>                             Set max virtual address space in GB for the pipeline. Only applies to local jobs
        --mempercore <NUM>                            Reserve enough threads for each job to ensure enough memory will be available, assuming each core on your cluster has at least this much memory available. Only applies to cluster jobmodes
        --maxjobs <NUM>                               Set max jobs submitted to cluster at one time. Only applies to cluster jobmodes
        --jobinterval <NUM>                           Set delay between submitting jobs to cluster, in ms. Only applies to cluster jobmodes
        --overrides <PATH>                            The path to a JSON file that specifies stage-level overrides for cores and memory. Finer-grained than --localcores, --mempercore and --localmem. Consult https://support.10xgenomics.com/ for an example
                                                      override file
        --uiport <PORT>                               Serve web UI at http://localhost:PORT
        --disable-ui                                  Do not serve the web UI
        --noexit                                      Keep web UI running after pipestance completes or fails
        --nopreflight                                 Skip preflight checks
    -h, --help                                        Print help information
        '''

        with open(make_count_path, "w") as f:
            if(len(self.fastq_path)!=1):
                print("For the count pipeline, you should have one fastq directory!")
                return # should be exit probably
            f.write(
            f"""#!/usr/bin/bash

    #SBATCH --job-name={''}
    #SBATCH --output={''}
    #SBATCH --partition={''}
    #SBATCH --mem={''}G
    #SBATCH --cpus-per-task={''}

    {self.aligner_software_path} count \
    --id={self.id} \
    --fastqs={self.fastq_folders_name[0]} \
    --sample={self.sample_name} \
    --transcriptime={self.aligner_ref_genome_path} \
    --maxjobs={MAX_JOBS} \
    --localcores={CPUS_PER_TASK} \
    --localmem={MEM_PER_TASK}

    mv {self.id} {H5_PATH}
    mv {make_count_path} {H5_PATH}

            """
            )



    def multiplex(self):

        '''
cellranger-multi 
Analyze multiplexed data or combined gene expression/immune profiling/feature barcode data

USAGE:
    cellranger multi [OPTIONS] --id <ID> --csv <CSV>
OPTIONS:
        --id <ID>               A unique run id and output folder name [a-zA-Z0-9_-]+
        --description <TEXT>    Sample description to embed in output files [default: ]
        --csv <CSV>             Path of CSV file enumerating input libraries and analysis parameters
        --dry                   Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop

        --jobmode <MODE>        Job manager to use. Valid options: local (default), sge, lsf, slurm or path to a .template file. Search for help on "Cluster Mode" at support.10xgenomics.com for more details on configuring the pipeline to use a compute
                                cluster [default: local]
        --localcores <NUM>      Set max cores the pipeline may request at one time. Only applies to local jobs
        --localmem <NUM>        Set max GB the pipeline may request at one time. Only applies to local jobs
        --localvmem <NUM>       Set max virtual address space in GB for the pipeline. Only applies to local jobs
        --mempercore <NUM>      Reserve enough threads for each job to ensure enough memory will be available, assuming each core on your cluster has at least this much memory available. Only applies to cluster jobmodes
        --maxjobs <NUM>         Set max jobs submitted to cluster at one time. Only applies to cluster jobmodes
        --jobinterval <NUM>     Set delay between submitting jobs to cluster, in ms. Only applies to cluster jobmodes
        --overrides <PATH>      The path to a JSON file that specifies stage-level overrides for cores and memory. Finer-grained than --localcores, --mempercore and --localmem. Consult https://support.10xgenomics.com/ for an example override file
        --uiport <PORT>         Serve web UI at http://localhost:PORT
        --disable-ui            Do not serve the web UI
        --noexit                Keep web UI running after pipestance completes or fails
        --nopreflight           Skip preflight checks

        '''

        # Here we create a compatible csv file for the multi pipeline.

#        if len(config_selected_pipeline['fastq_path']) == 0 or len(config_selected_pipeline['fastq_path']) != len(config_selected_pipeline['fastq_folders_name']) or len(config_selected_pipeline['feature_type'])  != len(config['fastq_folders_name']): # None parameters should be entered as condition as well
#            print("Your config is not valid!, plz try again later")
#            return

        self.multi_csv = pd.DataFrame(columns=['[gene-expression]',None, None, None])
        col = ['[gene-expression]',None, None, None]
        self.multi_csv.loc[0] = ['reference',None, None, None]
        row_reference_index = self.multi_csv[self.multi_csv.iloc[:, 0] == 'reference'].index[0]
        self.multi_csv.iloc[row_reference_index, 1] = self.aligner_ref_genome_path
        if self.create_bam is True:
            self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["create-bam", "true", None, None]], columns=col)], ignore_index=True)
        else:
            self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["create-bam", "false", None, None]], columns=col)], ignore_index=True)
        if self.R1_gex_len is not None and self.R2_gex_len is not None:
            self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["r1-length", self.R1_gex_len, None, None]], columns=col)], ignore_index=True)
            self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["r2-length", self.R2_gex_len, None, None]], columns=col)], ignore_index=True)
        if self.expected_cells is not None:
            self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["expected-cells", self.expected_cells, None, None]], columns=col)], ignore_index=True)
        if self.include_introns is not None:       
            self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["include-introns", self.include_introns, None, None]], columns=col)], ignore_index=True)
           
        if self.aligner_ref_vdj_path != None:
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[vdj]', None, None, None]], columns=col)], ignore_index=True)
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['reference', self.aligner_ref_vdj_path, None, None]], columns=col)], ignore_index=True)
            if self.R1_vdj_len is not None and self.R2_vdj_len is not None:
                self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["r1-length", self.R1_vdj_len, None, None]], columns=col)], ignore_index=True)
                self.multi_csv = self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([["r2-length", self.R2_vdj_len, None, None]], columns=col)], ignore_index=True)

        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[libraries]', None, None, None]], columns=col)], ignore_index=True)
        variable_libraries_names = ['fastq_id', 'fastqs', 'lanes', 'feature_types']
        variable_libraries_lists = list(zip(self.sample_names, self.fastq_path, self.lanes_used, self.feature_types))
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([variable_libraries_names], columns=col)], ignore_index=True)
        library_data_df = pd.DataFrame(variable_libraries_lists, columns=col)
        self.multi_csv = pd.concat([self.multi_csv, library_data_df], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        if (self.multiplexing_method == 'cmo_barcode'):
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[samples]', None, None, None]], columns=col)], ignore_index=True)
            variable_samples_names = ['sample_id', 'cmo_ids', 'description', None]
            variable_samples_list  = list(zip(self.sample_id_cmo,self.cmo_id, [None] * len(self.sample_id_cmo),[None] * len(self.sample_id_cmo)))
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([variable_samples_names], columns=col)], ignore_index=True)
            samples_data_df = pd.DataFrame(variable_samples_list, columns=col)
            self.multi_csv = pd.concat([self.multi_csv, samples_data_df], ignore_index=True)
            self.multi_csv.to_csv(MULTI_CSV_PATH, header=False, index=False)
        elif (self.multiplexing_method == 'feature_barcode'):
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[feature]', None, None, None]], columns=col)], ignore_index=True)
            if(self.feature_reference_csv != 'Default'):
                self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['reference', self.feature_reference_csv, None, None]], columns=col)], ignore_index=True)
            else:
                sample_sheet = pd.DataFrame([self.hto_id,self.hto_names,self.hto_read,self.hto_pattern,self.hto_sequence,self.HTO_feature_type]).transpose()
                sample_sheet.columns = ['id','name','read','pattern','sequence','feature_type']
                sample_sheet.to_csv(feature_reference_path, index=False)
                self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['reference', feature_reference_path, None, None]], columns=col)], ignore_index=True)
        else:
            print("Proceeding without hashing methods")
        
        self.multi_csv.to_csv(MULTI_CSV_PATH, index=False)

        with open(make_multi_path, "w+") as f:
            f.write(
            f"""#!/usr/bin/bash
    
{self.aligner_software_path} multi\
--id={self.id}\
--csv={MULTI_CSV_PATH}\
--maxjobs={self.jobs_number}\
--localcores={self.cpus_number}\
--localmem={self.total_memory_used}

            """ 
            )
        if self.running_machine == "Wexac":
            with open(run_file_path, "w+") as g:
                g.write(
            f"""
#BSUB -J make_multi
#BSUB -q {self.queue}
#BSUB -R rusage[mem={self.total_memory_used}G]
#BSUB -R affinity[thread*{jobs_number}]
#BSUB -R "span[hosts=1]"
#BSUB -oo {multi_log}
#BSUB -eo {multi_error_log} 

sh {make_multi_path}

            """

            )
    
    def demultiplex(self):    
        
        if self.demultiplex_method == 'hashsolo':
            for index in range(len(self.sample_names)):
                adata = sc.read_10x_h5(self.adata_path[index], gex_only=False)
                sample, hashtags = list(self.sample_id_to_expected_barcodes.items())[index]
                number_of_noise_barcodes = len(hashtags) - 1
                adata = demultiplex.hashsolo({sample: hashtags}, sample=sample, adata=adata, number_of_noise_barcodes=number_of_noise_barcodes, plot=False)
                
                adata.var[GENE_SYMBOLS_KEY] = adata.var.index
                adata.var.set_index(GENE_IDS_KEY, inplace=True)

                adata.obs[SAMPLE_ID_KEY] = sample
                adata.obs.index = adata.obs.index + "-" + adata.obs[SAMPLE_ID_KEY]

                adata_path = os.path.join(DEMULTIPLEXED_H5ADS_PATH, f"{sample}.h5ad")
                adata.write(adata_path)

        elif self.demultiplex_method == 'demultiplex2':
            log_file = open(demultiplex_log, "w+")
            os.dup2(log_file.fileno(), sys.stdout.fileno())
            for index in range(len(self.sample_names)):
                adata = sc.read_10x_h5(self.adata_path[index], gex_only=False)
                sample, hashtags = list(self.sample_id_to_expected_barcodes.items())[index]
                adata = demultiplex.demultiplex2({sample: hashtags}, sample=sample, adata=adata, plot=False)

                adata.var[GENE_SYMBOLS_KEY] = adata.var.index
                adata.var.set_index(GENE_IDS_KEY, inplace=True)

                adata.obs[SAMPLE_ID_KEY] = sample
                adata.obs.index = adata.obs.index + "-" + adata.obs[SAMPLE_ID_KEY]

                adata_path = os.path.join(DEMULTIPLEXED_H5ADS_PATH, f"{sample}.h5ad")
                adata.write(adata_path)
            
            sys.stdout = sys.__stdout__
            log_file.close()

        else:
            print("suupppp, check out you demultiplex method, it's not valid!!")

    def aggregate(self):
        '''
        cellranger-aggr 
Aggregate data from multiple Cell Ranger runs

USAGE:
    cellranger aggr [OPTIONS] --id <ID> --csv <CSV>

OPTIONS:    
        --id <ID>               A unique run id and output folder name [a-zA-Z0-9_-]+
        --description <TEXT>    Sample description to embed in output files [default: ]
        --csv <CSV>             Path of CSV file enumerating 'cellranger count/vdj/multi' outputs
        --normalize <MODE>      Library depth normalization mode [default: mapped] [possible values: mapped, none]
        --dry                   Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
        --nosecondary           Disable secondary analysis, e.g. clustering
        --jobmode <MODE>        Job manager to use. Valid options: local (default), sge, lsf, slurm or path to a .template file. Search for help on "Cluster Mode" at support.10xgenomics.com for more
                                details on configuring the pipeline to use a compute cluster [default: local]
        --localcores <NUM>      Set max cores the pipeline may request at one time. Only applies to local jobs
        --localmem <NUM>        Set max GB the pipeline may request at one time. Only applies to local jobs
        --localvmem <NUM>       Set max virtual address space in GB for the pipeline. Only applies to local jobs
        --mempercore <NUM>      Reserve enough threads for each job to ensure enough memory will be available, assuming each core on your cluster has at least this much memory available. Only applies to
                                cluster jobmodes
        --maxjobs <NUM>         Set max jobs submitted to cluster at one time. Only applies to cluster jobmodes
        --jobinterval <NUM>     Set delay between submitting jobs to cluster, in ms. Only applies to cluster jobmodes
        --overrides <PATH>      The path to a JSON file that specifies stage-level overrides for cores and memory. Finer-grained than --localcores, --mempercore and --localmem. Consult
                                https://support.10xgenomics.com/ for an example override file
        --uiport <PORT>         Serve web UI at http://localhost:PORT
        --disable-ui            Do not serve the web UI
        --noexit                Keep web UI running after pipestance completes or fails
        --nopreflight           Skip preflight checks
    -h, --help                  Print help information
        '''
        pass

    def make_anndatas(self):
        pass

    def make_custom_reference(self):
        with open(make_reference_path, "w+") as f:
            f.write(
            f"""#!/usr/bin/bash

    # SBATCH --job-name={''}
    # SBATCH --output={''}
    # SBATCH --partition={''}
    # SBATCH --mem={''}G
    # SBATCH --cpus-per-task={''}
    
    bsub \ 
    -q gsla_high_gpu \ 
    -gpu num=1:j_exclusive=yes \ 
    -R rusage[mem={self.total_memory_used}G] \ 
    -R affinity[thread*{self.cpus_number}] \ 
    -Is /bin/bash
    
    {self.aligner_software_path} mkgtf \ 
        {self.original_gtf_file_path} \ 
        {self.customized_gtf_file_path}"""
            )
            for attribute in self.filter_gtf_attribute:
                if attribute is not None:
                    f.write(f" \\\n        --attribute=gene_biotype:{attribute}")
                else:
                    f.write(f"\\\n")
                    break
        
            f.write(
                f"""

     {self.aligner_software_path} mkref \ 
        --genome={self.genome_name} \ 
        --fasta={self.original_fasta_file_path} \ 
        --genes={self.customized_gtf_file_path} \ 
        """
            )

    def cellbender(self):
        pass
