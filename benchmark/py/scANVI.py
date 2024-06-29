#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/11/30 17:40

# @Author : Fanjie
"""
import matplotlib.pyplot as plt
from scvi.model.utils import mde
import scanpy as sc
import scvi
import os

def run_scANVI(in_h5ad, out_h5ad, out_umap, out_scanvi, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start scANVI integration')
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
    print("Setup scANVI model")
    scvi.model.SCANVI.setup_anndata(adata, layer="counts", labels_key='map_celltype', unlabeled_category="Unknown", batch_key='batch')
    vae = scvi.model.SCANVI(adata)
    vae.train(max_epochs=20, n_samples_per_label=100)
    adata.obsm["X_scANVI"] = vae.get_latent_representation()
    adata.obsm["X_mde_scANVI"] = mde(adata.obsm["X_scANVI"])
    sc.pl.embedding(
        adata, basis="X_mde_scANVI", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['scANVI_batch', 'scANVI_ann_level1', 'scANVI_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_scanvi, dpi=300, format='pdf', bbox_inches='tight')
    sc.pp.neighbors(adata, use_rep="X_scANVI", key_added='scANVI', n_neighbors=30)
    sc.tl.umap(adata, neighbors_key='scANVI', min_dist=0.3)
    sc.pl.umap(
        adata, neighbors_key='scANVI', color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['scANVI_batch', 'scANVI_ann_level1', 'scANVI_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, format='pdf', bbox_inches='tight')
    adata.obsm['X_umapscANVI'] = adata.obsm['X_umap']
    print("Save output")
    adata.write_h5ad(out_h5ad)
    print("Done scANVI")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'scanvi.h5ad')
out_umap = os.path.join(figdir, 'scanvi_umap.pdf')
out_scanvi = os.path.join(figdir, 'scanvi_embed.pdf')

run_scANVI(in_h5ad, out_h5ad, out_umap, out_scanvi, 2000)