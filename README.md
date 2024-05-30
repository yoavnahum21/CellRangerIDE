# CellRangerIDE
This repo includes the future program to come for running various cellranger pipelines and additional algorithms. Have a look in the this README file!

Hello CellRanger users!

Here you can use various cellranger pipelines of 10x genomics.

Hopefully we will set an infrastracture that will make those pipelines more accessable 
and also new algorithems which are not included in 10x genomics.

The API should contain:

1) Converting BCL to Fastqs
2) Count pipelines
3) Multiplexing
4) Demultiplexing 
5) Quality checks
6) Dealing with various samples: GEX, VDJ, Antibodies etc.
7) PCA & TSNE & UMAP presentations
8) Optional: Preprocessing algorithems


All you need to do:

1) Fill out the config file in the const_files folder
2) If you wish to create .h5 you should add as well a reference transcriptome.
   You can find the Reference transcriptome on wexac in the next directory:
   "/home/labs/nyosef/yoavnah/CellRangerIDE/const_files/transcriptome"
3) If you use multiple samples, you should use "cellranger multi", therefore you have to fill out the multi_template.csv as well.
