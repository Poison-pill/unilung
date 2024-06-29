#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/12/23 23:04

# @Author : Fanjie
"""
import matplotlib.pyplot as plt
import scanpy as sc
import cellhint
import numpy as np
import pickle
import os

def run_cellhint(in_h5ad, out_h5ad, out_umap, out_cellhint, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start cellhint integration')
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
    print("Cellhint")
    cellhint.integrate(adata, 'batch', 'map_celltype')
    sc.pp.neighbors(adata, use_rep='X_pca', key_added='cellhint', n_neighbors=15)
    sc.tl.umap(adata, neighbors_key='cellhint', min_dist=0.3)
    sc.pl.embedding(
        adata, basis="X_pca", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['cellhint_batch', 'cellhint_ann_level1', 'cellhint_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_cellhint, dpi=300, format='pdf', bbox_inches='tight')
    sc.pl.umap(
        adata, neighbors_key='cellhint', color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['cellhint_batch', 'cellhint_ann_level1', 'cellhint_ann_level2']
    )
    plt.savefig(out_umap, dpi=300,  bbox_inches='tight')
    adata.obsm['X_umapcellhint'] = adata.obsm['X_umap']
    print("Save output")
    with open(os.path.join(outdir, "cellhint.pickle"), "wb") as file:
        pickle.dump(adata, file)
    adata.write(out_h5ad)
    print("Done cellhint")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'cellhint.h5ad')
out_umap = os.path.join(figdir, 'cellhint_umap.pdf')
out_cellhint = os.path.join(figdir, 'cellhint_embed.pdf')

run_cellhint(in_h5ad, out_h5ad, out_umap, out_cellhint, 2000)

def run_cellhint_harmonize(in_h5ad, out_h5ad, out_umap, out_cellhint, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start cellhint integration')
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
    print("Cellhint")
    alignment = cellhint.harmonize(adata, 'batch', 'map_celltype')
    adata.obs[['reannotation', 'group']] = alignment.reannotation[['reannotation', 'group']].loc[adata.obs_names]
    cellhint.integrate(adata, 'batch', 'reannotation')
    sc.pp.neighbors(adata, use_rep='X_pca', key_added='cellhint_harm', n_neighbors=30)
    sc.tl.umap(adata, neighbors_key='cellhint_harm', min_dist=0.3)
    sc.pl.embedding(
        adata, basis="X_pca", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['cellhint_harm_batch', 'cellhint_harm_ann_level1', 'cellhint_harm_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_cellhint, dpi=300, format='pdf', bbox_inches='tight')
    sc.pl.umap(
        adata, neighbors_key='cellhint_harm', color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['cellhint_harm_batch', 'cellhint_harm_ann_level1', 'cellhint_harm_ann_level2']
    )
    plt.savefig(out_umap, dpi=300,  bbox_inches='tight')
    adata.obsm['X_umapcellhint_harm'] = adata.obsm['X_umap']
    print("Save output")
    with open(os.path.join(outdir, "cellhint_harm.pickle"), "wb") as file:
        pickle.dump(adata, file)
    adata.write(out_h5ad)
    print("Done cellhint")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'cellhint_harm.h5ad')
out_umap = os.path.join(figdir, 'cellhint_harm_umap.pdf')
out_cellhint = os.path.join(figdir, 'cellhint_harm_embed.pdf')

run_cellhint_harmonize(in_h5ad, out_h5ad, out_umap, out_cellhint, 2000)