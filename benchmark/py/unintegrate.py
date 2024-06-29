#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/11/30 17:33

# @Author : Fanjie
"""
import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import os

dir = 'inte_bench'
with open(os.path.join(dir, "raw_py.pickle"), "rb") as file:
    adata = pickle.load(file)

"""
'n_counts': number of counts per cell (count depth)
'log_counts': ``np.log`` of counts per cell
'n_genes': number of counts per gene
'mito_frac': percent of mitochondrial gene counts as an indicator of cell viability
"""

count_matrix = adata.X
adata.obs["n_counts"] = count_matrix.sum(1)
adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
adata.obs["n_genes"] = (count_matrix > 0).sum(1)


# for each cell compute fraction of counts in mito genes vs. all genes
mito_genes = adata.var_names.str.match("^MT-")
mt_sum = np.sum(adata[:, mito_genes].X, axis=1)
total_sum = np.sum(adata.X, axis=1)
adata.obs["mito_frac"] = mt_sum / total_sum

# save raw data object
# TypeError: Can't implicitly convert non-string objects to strings
# Error raised while writing key 'unilung_sample' of <class 'h5py._hl.group.Group'> to /
# ad.obs['column_name'] = ad.obs['column_name'].astype(str)
adata.write_h5ad(os.path.join(dir, 'rawcount.h5ad'))

# with open(os.path.join(dir, "raw_py.pickle"), "wb") as file:
#     pickle.dump(adata, file)


# unintegrated
import matplotlib.pyplot as plt
import scanpy as sc
import os

def run_raw(in_h5ad, out_h5ad, out_umap, out_unintegrate, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_unintegrate: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start unintegrate')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
    adata = sc.read_h5ad(in_h5ad)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=top_genes,
        layer="counts",
        batch_key='batch',
        subset=True
    )
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.embedding(
        adata, basis="umap", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['unintegrate_batch', 'unintegrate_ann_level1', 'unintegrate_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_unintegrate, dpi=300, format='pdf', bbox_inches='tight')
    sc.pl.umap(
        adata, color=['batch', 'map_ann_level1', 'map_ann_level2'],
        ncols=1, title=['unintegrate_batch', 'unintegrate_ann_level1', 'unintegrate_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, format='pdf', bbox_inches='tight')
    adata.obsm['X_umapraw'] = adata.obsm['X_umap']
    print("Save output")
    adata.write_h5ad(out_h5ad)
    print("Done unintegerated")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'unintegrate.h5ad')
out_umap = os.path.join(figdir, 'uninteg_umap.pdf')
out_integrate = os.path.join(figdir, 'uninteg_embed.pdf')

run_raw(in_h5ad, out_h5ad, out_umap, out_integrate, 2000)
