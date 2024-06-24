from _utils import * 
import os
import shutil
import subprocess
import hashsolo
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

MEM_PER_TASK = 70
CPUS_PER_TASK = 4
MAX_JOBS = 400
SAMPLE_NAME = "EXAMPLE01" 
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
            os.makedirs(FASTQC_PATH)
            os.makedirs(H5_PATH)
            os.makedirs(H5ADS_PATH)
            os.makedirs(BEFORE_DEMULTI_H5ADS_PATH)
            os.makedirs(DEMULTIPLEXED_H5ADS_PATH)
            os.makedirs(OUTPUT_PATH)

        else:
            print("This is an exiting project!")

        self.donor = config['donor']

        self.seq_run = config['seq_run']

        self.id = config['id']

        self.output_destination = config['output_destination']
        if self.output_destination == None:
            self.output_destinatxion = OUTPUT_PATH
        
        
        self.aligner_software_path = config['aligner_software_path']
        if self.aligner_software_path == None:
            self.aligner_software_path = "/apps/RH7U2/general/cellranger/7.1.0/cellranger"

        self.aligner_ref_genome_path = config['alignment_ref_genome_file']
        if self.aligner_ref_genome_path == None:
            self.aligner_ref_genome_path = gex_reference_path

        self.aligner_ref_vdj_path = config['alignment_ref_vdj_file']
        if self.aligner_ref_vdj_path == None:
            self.aligner_ref_vdj_path = vdj_reference_path

        self.user = config['weizmann_user']

        self.shell_file = config['shell_file']
        if self.shell_file == None:
            self.shell_file = OUTPUT_PATH + str(self.id) + "_cellrangerpipe.sh"
        
        self.Sample_sheet_address = config['Sample_sheet_address']
        
        self.fastq_folders_name = config['fastq_folders_name']

        self.fastq_path = config['fastq_path']

        self.lanes_used = config['lanes_used']

        self.feature_type = config['feature_type']

        self.multi_csv = None
        
        self.multiplexing_method = config['multiplexing_method']
        if ((self.multiplexing_method != 'feature_barcode') and (self.multiplexing_method != 'cmo_barcode') and (self.multiplexing_method != None)):
            print("Your multiplexing method isn't valid")
            exit(1)

        self.feature_ref_csv = config["feature_reference_csv"]
        if self.feature_ref_csv == None:
            self.hto_id = config["hto_id"]
            self.hto_names = config["hto_names"]
            if self.hto_names[0] == None:
                self.hto_names = self.hto_id
            self.hto_read = config["hto_read"]
            if self.hto_read[0] == None:
                self.hto_read = ['R2'] * len(self.hto_read)
            self.hto_pattern = config["hto_pattern"]
            if self.hto_pattern[0] == None:
                self.hto_pattern = ['5PNNNNNNNNNN(BC)'] * len(self.hto_pattern)
            self.hto_sequence = config["hto_sequence"]
            self.HTO_feature_type = config["HTO_feature_type"]
            if self.HTO_feature_type[0] == None:
                self.HTO_feature_type = ['Antibody Capture'] * len(self.HTO_feature_type)

        
    
    def __str__(self) -> str:
        pass


    def read_program_from_csv(self): 

        #TBD - decide about a format that should read the next files names (Use: self.program)

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

    cellranger mkfastq --run={''} --csv={''} --output-dir={''} --maxjobs={''} --localcores={''} --localmem={''}
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
            f.write(
            f"""#!/usr/bin/bash

    #SBATCH --job-name={''}
    #SBATCH --output={''}
    #SBATCH --partition={''}
    #SBATCH --mem={''}G
    #SBATCH --cpus-per-task={''}

    cellranger count --run={''} --csv={''} --output-dir={''} --maxjobs={''} --localcores={''} --localmem={''}
            """
            )


    def make_samplesheet(self):
        pass


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
        if len(config['fastq_path']) == 0 or len(config['fastq_path']) != len(config['fastq_folders_name']) or len(config['feature_type'])  != len(config['fastq_folders_name']): # None parameters should be entered as condition as well
            print("Your config is not valid!, plz try again later")
            return
        self.multi_csv = pd.DataFrame(columns=['[gene-expression]',None, None, None])
        col = ['[gene-expression]',None, None, None]
        self.multi_csv.loc[0] = ['reference',None, None, None]
        row_reference_index = self.multi_csv[self.multi_csv.iloc[:, 0] == 'reference'].index[0]
        self.multi_csv.iloc[row_reference_index, 1] = self.aligner_ref_genome_path
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[libraries]', None, None, None]], columns=col)], ignore_index=True)
        variable_libraries_names = ['fastq_id', 'fastqs', 'lanes', 'feature_types']
        variable_libraries_lists = list(zip(self.fastq_folders_name, self.fastq_path, self.lanes_used, self.feature_type))
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([variable_libraries_names], columns=col)], ignore_index=True)
        library_data_df = pd.DataFrame(variable_libraries_lists, columns=col)
        self.multi_csv = pd.concat([self.multi_csv, library_data_df], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([[None, None, None, None]], columns=col)], ignore_index=True)
        if (self.multiplexing_method == 'cmo_barcode'):
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[samples]', None, None, None]], columns=col)], ignore_index=True)
            variable_samples_names = ['sample_id', 'cmo_ids', 'description', None]
            variable_samples_list  = list(zip(self.fastq_folders_name,[None, None, None, None], [None, None, None, None], [None, None, None, None]))
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([variable_samples_names], columns=col)], ignore_index=True)
            samples_data_df = pd.DataFrame(variable_samples_list, columns=col)
            self.multi_csv = pd.concat([self.multi_csv, samples_data_df], ignore_index=True)
            self.multi_csv.to_csv(MULTI_CSV_PATH, header=False, index=False)
        elif (self.multiplexing_method == 'feature_barcode'):
            self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['[feature]', None, None, None]], columns=col)], ignore_index=True)
            if(self.feature_ref_csv != "Default"):
                self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['reference', self.feature_ref_csv, None, None]], columns=col)], ignore_index=True)
            else:
                sample_sheet = pd.DataFrame([self.hto_id,self.hto_names,self.hto_read,self.hto_pattern,self.hto_sequence,self.HTO_feature_type]).transpose()
                sample_sheet.columns = ['id','name','read','pattern','sequence','feature_type']
                sample_sheet.to_csv(feature_reference_path, index=False)
                self.multi_csv = pd.concat([self.multi_csv, pd.DataFrame([['reference', feature_reference_path, None, None]], columns=col)], ignore_index=True)
        else:
            print("you chose wrong multiplexing method")
            return
        self.multi_csv.to_csv(MULTI_CSV_PATH, index=False)

        with open(make_multi_path, "w+") as f:
            f.write(
            f"""#!/usr/bin/bash

    # SBATCH --job-name={''}
    # SBATCH --output={''}
    # SBATCH --partition={''}
    # SBATCH --mem={''}G
    # SBATCH --cpus-per-task={''}

    {self.aligner_software_path} multi --id={self.id} --csv={MULTI_CSV_PATH} --maxjobs={MAX_JOBS} --localcores={CPUS_PER_TASK} --localmem={MEM_PER_TASK}
    mv {self.id} {H5_PATH}
    mv {make_multi_path} {H5_PATH}
            """
            )
    
    # def demultiplex(self, sample: str, adata: ad.AnnData, plot: bool = True):
    def demultiplex(self):    
        adata = sc.read_10x_h5(os.path.join(H5_PATH, "sample_filtered_feature_bc_matrix.h5"), gex_only=False)
        sample = "B16"
        adata = hashsolo.demultiplex(sample,adata,plot=False)
        adata.var[GENE_SYMBOLS_KEY] = adata.var.index
        adata.var.set_index(GENE_IDS_KEY, inplace=True)

        adata.obs[SAMPLE_ID_KEY] = sample
        adata.obs.index = adata.obs.index + "-" + adata.obs[SAMPLE_ID_KEY]

        adata_path = os.path.join(H5ADS_PATH, f"{sample}.h5ad")
        adata.write(adata_path)

        

    def reanalayze(self):
        '''
        cellranger-reanalyze 
Re-run secondary analysis (dimensionality reduction, clustering, etc)

USAGE:
    cellranger reanalyze [OPTIONS] --id <ID> --matrix <MATRIX_H5>

OPTIONS:
        --id <ID>                      A unique run id and output folder name [a-zA-Z0-9_-]+
        --description <TEXT>           Sample description to embed in output files [default: ]
        --matrix <MATRIX_H5>           A feature-barcode matrix containing data for one genome. Should be the filtered version, unless using --force-cells
        --params <PARAMS_CSV>          A CSV file specifying analysis parameters. Optional
        --barcodes <BARCODES_CSV>      A CSV file containing a list of cell barcodes to use for reanalysis, e.g. barcodes exported from Loupe Browser. Optional
        --genes <GENES_CSV>            A CSV file containing a list of feature IDs to use for reanalysis. For gene expression, this should correspond to the gene_id field in the reference GTF should be \(e.g. ENSG... for ENSEMBL-based references\).
                                       Optional
        --exclude-genes <GENES_CSV>    A CSV file containing a list of feature IDs to exclude from reanalysis. For gene expression, this should correspond to the gene_id field in the reference GTF \(e.g., ENSG... for ENSEMBL-based references\). The
                                       exclusion is applied after --genes. Optional
        --agg <AGGREGATION_CSV>        If the input matrix was produced by 'aggr', you may pass the same aggregation CSV in order to retain per-library tag information in the resulting .cloupe file.  This argument is required to enable chemistry batch
                                       correction. Optional
        --force-cells <NUM>            Force pipeline to use this number of cells, bypassing cell calling algorithm. [MINIMUM: 10]
        --dry                          Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
        --jobmode <MODE>               Job manager to use. Valid options: local (default), sge, lsf, slurm or path to a .template file. Search for help on "Cluster Mode" at support.10xgenomics.com for more details on configuring the pipeline to use a
                                       compute cluster [default: local]
        --localcores <NUM>             Set max cores the pipeline may request at one time. Only applies to local jobs
        --localmem <NUM>               Set max GB the pipeline may request at one time. Only applies to local jobs
        --localvmem <NUM>              Set max virtual address space in GB for the pipeline. Only applies to local jobs
        --mempercore <NUM>             Reserve enough threads for each job to ensure enough memory will be available, assuming each core on your cluster has at least this much memory available. Only applies to cluster jobmodes
        --maxjobs <NUM>                Set max jobs submitted to cluster at one time. Only applies to cluster jobmodes
        --jobinterval <NUM>            Set delay between submitting jobs to cluster, in ms. Only applies to cluster jobmodes
        --overrides <PATH>             The path to a JSON file that specifies stage-level overrides for cores and memory. Finer-grained than --localcores, --mempercore and --localmem. Consult https://support.10xgenomics.com/ for an example override file
        --uiport <PORT>                Serve web UI at http://localhost:PORT
        --disable-ui                   Do not serve the web UI
        --noexit                       Keep web UI running after pipestance completes or fails
        --nopreflight                  Skip preflight checks
        '''
        with open(make_multi_path, "w") as f:
            f.write(
            f"""#!/usr/bin/bash

    # SBATCH --job-name={''}
    # SBATCH --output={''}
    # SBATCH --partition={''}
    # SBATCH --mem={''}G
    # SBATCH --cpus-per-task={''}

    cellranger multi --id={self.id} --csv={MULTI_CSV_PATH} --maxjobs={MAX_JOBS} --localcores={CPUS_PER_TASK} --localmem={MEM_PER_TASK}
            """
            )

    def multi_with_5prime(self):
        pass

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
        pass

    
    def run_basic_pipeline(self):
        with open(self.shell_file, "w") as f:
            f.write(
            f"""#!/usr/bin/bash

            # {self.aligner_software_path} mkfastq --run={BCL_PATH}  --output-dir={FASTQ_PATH} --maxjobs={MAX_JOBS} --localcores={CPUS_PER_TASK} --localmem={MEM_PER_TASK}
            {self.aligner_software_path} count --id={self.id} --transcriptome={self.aligner_ref_genome_path} --fastqs={FASTQ_PATH} --sample={SAMPLE_NAME} --maxjobs={MAX_JOBS} --localcores={CPUS_PER_TASK} --localmem={MEM_PER_TASK}
            """
            )
        