#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/11/29 18:43

# @Author : Fanjie
"""
## harmony
import matplotlib.pyplot as plt
import scanpy as sc
import pickle
import os

def run_harmony(in_h5ad, out_h5ad, out_umap, out_harmony, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start Harmony integration')
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
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.3)
    adata.obsm['X_umapraw'] = adata.obsm['X_umap']
    print("harmony")
    sc.external.pp.harmony_integrate(adata, key='batch', basis='X_pca')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', key_added='harmony', n_neighbors=30)
    sc.pl.embedding(
        adata, basis="X_pca_harmony", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['harmony_batch', 'harmony_ann_level1', 'harmony_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_harmony, dpi=300, format='pdf', bbox_inches='tight')
    sc.tl.umap(adata, neighbors_key='harmony', min_dist=0.3)
    sc.pl.umap(
        adata, color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['harmony_batch', 'harmony_ann_level1', 'harmony_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, bbox_inches='tight')
    print("Save output")
    with open(os.path.join(outdir, "merge.pickle"), "wb") as file:
        pickle.dump(in_ad, file)
    adata.write(out_h5ad)
    print("Done harmony")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'harmony.h5ad')
out_umap = os.path.join(figdir, 'harmony_umap.pdf')
out_harmony = os.path.join(figdir, 'harmony_embed.pdf')

run_harmony(in_h5ad, out_h5ad, out_umap, out_harmony, 2000)

