import anndata as ad
import os
import scanpy as sc
import pandas as pd
import subprocess
from _utils import CSV_PATH, R_RUNNER, OUTPUT_PATH, DEMULTIPLEXED_H5ADS_PATH
from ridgeplot import ridgeplot
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from _utils import (
    FEATURE_TYPE_KEY,
    ANTIBODY_KEY,
    CLASSIFICATION_KEY,
    H5_PATH
)

def hashsolo(sample_id_to_expected_barcodes: dict, sample: str, adata: ad.AnnData, number_of_noise_barcodes: int, plot: bool = True) -> ad.AnnData:

    bdata = adata.copy()
    all_barcodes = list(adata.var[adata.var[FEATURE_TYPE_KEY] == ANTIBODY_KEY].index.unique())

    expected_barcodes = sample_id_to_expected_barcodes[sample]
    n_obs_prev = adata.n_obs

    for barcode in expected_barcodes:
        values = np.squeeze(bdata.X[:, bdata.var.index == barcode].toarray())
        bdata.obs[barcode] = values
        bdata.obs[f"{barcode}_log"] = np.log1p(values)

    for barcode in all_barcodes:
        adata = adata[:, adata.var.index != barcode]

    sc.external.pp.hashsolo(bdata, expected_barcodes, number_of_noise_barcodes)

    df = bdata.obs[[
        "most_likely_hypothesis", 
        "cluster_feature", 
        "negative_hypothesis_probability", 
        "singlet_hypothesis_probability",
        "doublet_hypothesis_probability",
        CLASSIFICATION_KEY,
    ]]
    df.to_csv(os.path.join(H5_PATH, f"{sample}_hashsolo.csv"))

    if plot:
        plot_qc(sample, bdata, expected_barcodes)    ##### TBD

    mask = bdata.obs[CLASSIFICATION_KEY].isin(expected_barcodes)
    adata = adata[mask]
    bdata = bdata[mask]
    adata.obs[CLASSIFICATION_KEY] = bdata.obs[CLASSIFICATION_KEY]

    n_obs_post = adata.n_obs
    print(f"{sample}: {n_obs_prev} -> {n_obs_post}, removed {(n_obs_prev - n_obs_post)/n_obs_prev}")

    return adata

def demultiplex2(sample_id_to_expected_barcodes: dict, sample: str, adata: ad.AnnData, plot: bool = True) -> ad.AnnData:
    
    bdata = adata.copy()
    all_barcodes = list(adata.var[adata.var[FEATURE_TYPE_KEY] == ANTIBODY_KEY].index.unique())
    
    expected_barcodes = sample_id_to_expected_barcodes[sample]
    n_obs_prev = adata.n_obs

    for barcode in expected_barcodes:
        values = np.squeeze(bdata.X[:, bdata.var.index == barcode].toarray())
        bdata.obs[barcode] = values
        bdata.obs[f"{barcode}_log"] = np.log1p(values)

    for barcode in all_barcodes:
        adata = adata[:, adata.var.index != barcode]

    df = bdata.obs[expected_barcodes] # The problem is here let's continue tomorrow !! 
    csv_barcodes_path = os.path.join(CSV_PATH, "cell_hashing.csv") 
    df.to_csv(csv_barcodes_path, index=True)
    
    r_script_file = os.path.join(OUTPUT_PATH, "deMULTIplex2.R")
    
    with open(r_script_file, "w+") as f:
        f.write(
f"""
library(deMULTIplex2)
library(dplyr)
library(ggplot2)
log_file <- '{os.path.join(DEMULTIPLEXED_H5ADS_PATH, 'deMULTIplex_log.txt')}'
sink(log_file, append = TRUE) # Redirect output to log file

cat('Starting demultiplexing process...\\n')
mat <- read.csv('{csv_barcodes_path}')
rownames(mat) <- mat$X
mat$X <- NULL
colnames(mat) <- gsub('X', '', colnames(mat))
res <- demultiplexTags(mat, plot.path = '{DEMULTIPLEXED_H5ADS_PATH}', plot.name = 'demultiplex_plot', plot.diagnostics = TRUE)
prob_mtx_df <- as.data.frame(res$prob_mtx)
res_mtx_df <- as.data.frame(res$res_mtx)
final_assign_df <- as.data.frame(res$final_assign)
coefs <- as.data.frame(res$coefs)
write.csv(prob_mtx_df, '{os.path.join(DEMULTIPLEXED_H5ADS_PATH, 'prob_mtx_df_deMULTIplex2.csv')}')
write.csv(res_mtx_df, '{os.path.join(DEMULTIPLEXED_H5ADS_PATH, 'res_mtx_df_deMULTIplex2.csv')}')
write.csv(final_assign_df, '{os.path.join(DEMULTIPLEXED_H5ADS_PATH, 'final_assign_deMULTIplex2.csv')}')
write.csv(coefs, '{os.path.join(DEMULTIPLEXED_H5ADS_PATH, 'coefs_deMULTIplex2.csv')}')

"""
        )
    with open("/home/labs/nyosef/yoavnah/CellRangerIDE/Projects/Oier/File_Path/h5ads/demultiplexed/deMULTIplex_log.log", "w+") as g:
        process = subprocess.Popen([f"{R_RUNNER}", f"{r_script_file}"],
                                   stdout= g,
                                   stderr= subprocess.STDOUT)
        process.wait() 
        
    mask = pd.read_csv(os.path.join(DEMULTIPLEXED_H5ADS_PATH, 'final_assign_deMULTIplex2.csv'))
    for i in range(len(adata.obs)):
        if mask['Unnamed: 0'][i] != adata.obs.iloc[i]._name:
            mask = pd.concat([mask.iloc[:i], pd.DataFrame({'Unnamed: 0': adata.obs.iloc[i]._name, 'res$final_assign': ['negative']}), mask.iloc[i:]]).reset_index(drop=True)
    mask1 = mask["res$final_assign"].isin(expected_barcodes)
    adata = adata[mask1]
    mask = mask[mask['res$final_assign'].isin(expected_barcodes)]
    mask.reset_index(drop=True,inplace=True)
    mask = mask.set_index('Unnamed: 0')
    adata.obs[CLASSIFICATION_KEY] = mask['res$final_assign']

    n_obs_post = adata.n_obs
    print(f"{sample}: {n_obs_prev} -> {n_obs_post}, {((n_obs_prev - n_obs_post) * 100 /n_obs_prev):.2f}% of the cells were removed!")
    return adata

def plot_qc(sample: str, adata: ad.AnnData, barcodes: List[str]):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111, projection="3d" if len(barcodes) == 3 else None)
    for classification in adata.obs[CLASSIFICATION_KEY].unique():
        bdata = adata[adata.obs[CLASSIFICATION_KEY] == classification]
        ax.scatter(*[bdata.obs[barcode] for barcode in barcodes], label=classification)
    plt.legend()
    plt.show()
    path = os.path.join(H5_PATH, f"{sample}_scatter.png")
    fig.savefig(path, dpi=400)

    log_keys = [f"{barcode}_log" for barcode in barcodes]
    df = adata.obs[log_keys]
    fig = ridgeplot(
        samples=df.values.T,
        bandwidth=0.5,
        labels=barcodes,
        spacing=5/9
    )
    fig.show()
    path = os.path.join(H5_PATH, f"{sample}_ridge.png")
    fig.write_image(path)
