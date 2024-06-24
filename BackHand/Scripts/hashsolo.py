import anndata as ad
import os
import scanpy as sc
from ridgeplot import ridgeplot
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from _utils import (
    FEATURE_TYPE_KEY,
    ANTIBODY_KEY,
    CLASSIFICATION_KEY,
    SAMPLE_ID_TO_EXPECTED_BARCODES,
    H5_PATH
)

def demultiplex(sample: str, adata: ad.AnnData, plot: bool = True) -> ad.AnnData:

    bdata = adata.copy()
    all_barcodes = list(adata.var[adata.var[FEATURE_TYPE_KEY] == ANTIBODY_KEY].index.unique())

    expected_barcodes = SAMPLE_ID_TO_EXPECTED_BARCODES[sample]
    n_obs_prev = adata.n_obs

    for barcode in expected_barcodes:
        values = np.squeeze(bdata.X[:, bdata.var.index == barcode].toarray())
        bdata.obs[barcode] = values
        bdata.obs[f"{barcode}_log"] = np.log1p(values)

    for barcode in all_barcodes:
        adata = adata[:, adata.var.index != barcode]

    sc.external.pp.hashsolo(bdata, expected_barcodes, number_of_noise_barcodes=1)

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
