print('''        cellranger-reanalyze 
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
        --nopreflight                  Skip preflight checks''')