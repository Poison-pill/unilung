#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/11/30 17:38

# @Author : Fanjie
"""
import matplotlib.pyplot as plt
import scanpy as sc
import os
import scvi
from scvi.model.utils import mde


def run_scVI(in_h5ad, out_h5ad, out_umap, out_scvi, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start scVI integration')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
    adata = sc.read_h5ad(in_h5ad)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=top_genes,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        batch_key='batch',
        subset=True
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    print("Setup scVI model")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key='batch')
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, n_hidden=128, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    adata.obsm["X_mde_scVI"] = mde(adata.obsm["X_scVI"])
    sc.pl.embedding(
        adata, basis="X_mde_scVI", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['scVI_batch', 'scVI_ann_level1', 'scVI_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_scvi, dpi=300, format='pdf', bbox_inches='tight')
    sc.pp.neighbors(adata, use_rep="X_scVI", key_added='scVI', n_neighbors=30)
    sc.tl.umap(adata, neighbors_key='scVI', min_dist=0.3)
    sc.pl.umap(
        adata, neighbors_key='scVI', color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['scVI_batch', 'scVI_ann_level1', 'scVI_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, format='pdf', bbox_inches='tight')
    adata.obsm['X_umapscVI'] = adata.obsm['X_umap']
    print("Save output")
    adata.write_h5ad(out_h5ad)
    print("Done scVI")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'scvi.h5ad')
out_scvi = os.path.join(figdir, 'scvi_embed.pdf')
out_umap = os.path.join(figdir, 'scvi_umap.pdf')

run_scVI(in_h5ad, out_h5ad, out_umap, out_scvi, 2000)