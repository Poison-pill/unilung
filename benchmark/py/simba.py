#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/11/30 17:43

# @Author : Fanjie
"""
import scanpy as sc
from scvi.model.utils import mde
import matplotlib.pyplot as plt
import omicverse as ov
import pickle
import os


def run_simba(in_h5ad, out_h5ad, out_umap, out_simba, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_scvi: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start simba integration')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
    adata = sc.read_h5ad(in_h5ad)
    adata.var_names_make_unique()
    adata_simba = adata.copy()
    print("Setup simba model")
    workdir = '/public8/lilab/student/fjwu/btit1/inte_bench/models/result_simba'
    adata_simba = ov.single.pySIMBA(adata_simba, workdir)
    adata_simba.preprocess(batch_key='batch', min_n_cells=1, method='lib_size',
                           n_top_genes=top_genes, n_bins=5)
    adata_simba.gen_graph()
    adata_simba.train(num_workers=10)
    adata_new = adata_simba.batch_correction()
    adata.obsm['X_simba'] = adata_new.obsm['X_simba'].copy()
    adata.obsm["X_mde_simba"] = mde(adata.obsm["X_simba"])
    sc.pl.embedding(
        adata, basis="X_mde_simba", color=['batch', 'map_ann_level1', 'map_ann_level2'],
        ncols=1, title=['simba_batch', 'simba_ann_level1', 'simba_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_simba, dpi=300, format='pdf', bbox_inches='tight')
    sc.pp.neighbors(adata, use_rep="X_simba", key_added='simba', n_neighbors=30)
    sc.tl.umap(adata, neighbors_key='simba', min_dist=0.3)
    sc.pl.umap(
        adata, neighbors_key='simba', color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['simba_batch', 'simba_ann_level1', 'simba_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, format='pdf', bbox_inches='tight')
    adata.obsm['X_umapsimba'] = adata.obsm['X_umap']
    print("Save output")
    with open(os.path.join(outdir, "simba.pickle"), "wb") as file:
        pickle.dump(adata, file)
    adata.write_h5ad(out_h5ad)
    print("Done simba")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'raw.h5ad')
out_h5ad = os.path.join(outdir, 'simba.h5ad')
out_umap = os.path.join(figdir, 'simba_umap.pdf')
out_simba = os.path.join(figdir, 'simba_embed.pdf')

run_simba(in_h5ad, out_h5ad, out_umap, out_simba, 2000)