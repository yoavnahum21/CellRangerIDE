from _utils import * 
import os
import shutil
import subprocess
import pandas as pd

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
        if self.output_destination == "Default":
            self.output_destinatxion = OUTPUT_PATH
        
        
        self.aligner_software_path = config['aligner_software_path']
        if self.aligner_software_path == "Default":
            self.aligner_software_path = "cellranger"

        self.aligner_ref_genome_path = config['alignment_ref_genome_file']
        if self.aligner_ref_genome_path == "Default":
            self.aligner_ref_genome_path = gex_reference_path

        self.aligner_ref_vdj_path = config['alignment_ref_vdj_file']
        if self.aligner_ref_vdj_path == "Default":
            self.aligner_ref_vdj_path = vdj_reference_path

        self.user = ['weizmann_user']

        self.shell_file = config['shell_file']
        if self.shell_file == "Default":
            self.shell_file = OUTPUT_PATH + self.id + "_cellrangerpipe.sh"
        
        self.Sample_sheet_address = config['Sample_sheet_address']

        self.multi_csv = pd.read_csv(MULTI_TEMPLATE_CSV_PATH)
        
    
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
        with open(make_multi_path, "w") as f:
            f.write(
            f"""#!/usr/bin/bash

    # SBATCH --job-name={''}
    # SBATCH --output={''}
    # SBATCH --partition={''}
    # SBATCH --mem={''}G
    # SBATCH --cpus-per-task={''}

    cellranger multi --id={self.id} --csv={} --maxjobs={''} --localcores={''} --localmem={''}
            """
            )
    
    def demultiplex(self):

        '''
        --id <ID>               A unique run id and output folder name [a-zA-Z0-9_-]+
        --description <TEXT>    Sample description to embed in output files [default: ]
        --csv <CSV>             Path of CSV file enumerating 'cellranger count/vdj/multi' outputs
        --normalize <MODE>      Library depth normalization mode [default: mapped] [possible values: mapped, none]
        --nosecondary           Disable secondary analysis, e.g. clustering
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
        '''
        pass

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
        pass



    def make_anndatas(self):
        pass
    
    def run_basic_pipeline(self):
        with open(self.shell_file, "w") as f:
            f.write(
            f"""#!/usr/bin/bash

            # {self.aligner_software_path} mkfastq --run={BCL_PATH}  --output-dir={FASTQ_PATH} --maxjobs={MAX_JOBS} --localcores={CPUS_PER_TASK} --localmem={MEM_PER_TASK}
            {self.aligner_software_path} count --id={self.id} --transcriptome={self.aligner_ref_genome_path} --fastqs={FASTQ_PATH} --sample={SAMPLE_NAME} --maxjobs={MAX_JOBS} --localcores={CPUS_PER_TASK} --localmem={MEM_PER_TASK}
            """
            )
        