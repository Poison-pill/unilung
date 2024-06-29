#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/11/29 18:45

# @Author : Fanjie
"""
## Scanorama
import matplotlib.pyplot as plt
import scanpy as sc
import scanorama
import numpy as np
import pickle
import os

def run_scanorama(in_h5ad, out_h5ad, out_umap, out_scanorama, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start scanorama integration')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
    adata = sc.read_h5ad(in_h5ad)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=top_genes,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        batch_key='batch',
        subset=True
    )
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.3)
    adata.obsm['X_umapraw'] = adata.obsm['X_umap']
    print("Scanorama")
    sc.external.pp.scanorama_integrate(adata, key='batch', basis='X_pca', adjusted_basis='X_scanorama')
    sc.pp.neighbors(adata, use_rep='X_scanorama', key_added='scanorama', n_neighbors=30)
    sc.tl.umap(adata, neighbors_key='scanorama', min_dist=0.3)
    sc.pl.embedding(
        adata, basis="X_scanorama", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['scanorama_batch', 'scanorama_ann_level1', 'scanorama_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_scanorama, dpi=300, format='pdf', bbox_inches='tight')
    sc.pl.umap(
        adata, neighbors_key='scanorama', color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['scanorama_batch', 'scanorama_ann_level1', 'scanorama_ann_level2']
    )
    plt.savefig(out_umap, dpi=300,  bbox_inches='tight')
    adata.obsm['X_umapscanorama'] = adata.obsm['X_umap']
    print("Save output")
    with open(os.path.join(outdir, "scanorama.pickle"), "wb") as file:
        pickle.dump(adata, file)
    adata.write(out_h5ad)
    print("Done scanorama")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'scanorama.h5ad')
out_umap = os.path.join(figdir, 'scanorama_umap.pdf')
out_scanorama = os.path.join(figdir, 'scanorama_embed.pdf')

run_scanorama(in_h5ad, out_h5ad, out_umap, out_scanorama, 2000)
