#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/12/4 10:37

# @Author : Fanjie
"""
## ComBat
import scanpy as sc
from scvi.model.utils import mde
import matplotlib.pyplot as plt
import omicverse as ov
import pickle
import os

def run_combat(in_h5ad, out_h5ad, out_umap, out_combat, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_combat: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start ComBat integration')
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
    adata_combat = ov.single.batch_correction(adata, batch_key='batch',
                                              methods='combat', n_pcs=50)
    adata.obsm['X_combat'] = adata_combat.obsm['X_combat'].copy()
    adata.obsm["X_mde_combat"] = mde(adata.obsm["X_combat"])
    sc.pl.embedding(
        adata, basis="X_mde_combat", color=['batch', 'map_ann_level1', 'map_ann_level2'],
        ncols=1, title=['combat_batch', 'combat_ann_level1', 'combat_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_combat, dpi=300, format='pdf', bbox_inches='tight')
    sc.pp.neighbors(adata, use_rep="X_combat", key_added='combat', n_neighbors=30)
    sc.tl.umap(adata, neighbors_key='combat', min_dist=0.3)
    sc.pl.umap(
        adata, neighbors_key='combat', color=['batch', 'map_ann_level1', 'map_ann_level2'],
        ncols=1, title=['combat_batch', 'combat_ann_level1', 'combat_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, format='pdf', bbox_inches='tight')
    adata.obsm['X_umapcombat'] = adata.obsm['X_umap']
    print("Save output")
    with open(os.path.join(outdir, "combat.pickle"), "wb") as file:
        pickle.dump(adata, file)
    adata.write_h5ad(out_h5ad)
    print("Done combat")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'combat.h5ad')
out_umap = os.path.join(figdir, 'combat_umap.pdf')
out_combat = os.path.join(figdir, 'combat_embed.pdf')

run_combat(in_h5ad, out_h5ad, out_umap, out_combat, 2000)