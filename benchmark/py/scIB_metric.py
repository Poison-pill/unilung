#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2023/12/7 10:37

# @Author : Fanjie
"""
import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import scib
import os


def get_resolutions(n=10, min=0.1, max=2):
    min = np.max([1, int(min * 10)])
    max = np.max([min, max * 5])
    frac = n / 5
    return [frac * x / n for x in range(min, max + 1)]


def caculate_lisi(adata, batch_key, label_key, cluster_key, embed, type_, scaled=True):
    try:
        from rpy2.robjects import r
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        LISI = importr("lisi")
    except ModuleNotFoundError as e:
        raise e

    X = pd.DataFrame(adata.obsm[embed])
    meta = adata.obs
    X.index = meta.index

    if type_ == 'ilisi':
        scores = LISI.compute_lisi(X, meta, [batch_key, cluster_key])
        scores = pandas2ri.rpy2py(scores)
        scores_med = np.nanmedian(scores[batch_key])
        if scaled:
            lisi = (scores_med - np.min(scores[batch_key])) / (np.max(scores[batch_key]) - np.min(scores[batch_key]))
        else:
            lisi = scores_med
    else:
        scores = LISI.compute_lisi(X, meta, [label_key, cluster_key])
        scores = pandas2ri.rpy2py(scores)
        scores_med = np.nanmedian(scores[label_key])
        if scaled:
            lisi = (np.min(scores[label_key]) - scores_med) / (np.max(scores[label_key]) - np.min(scores[label_key]))
        else:
            lisi = scores_med

    return lisi


def run_scib_metrics(in_h5ad, raw_h5ad, datadir, method, batch_key, label_key, cluster_key, out_csv, out_h5ad, out_sil):
    embedding_keys = {"Unintegrate": "X_pca", "Harmony": "X_pca_harmony", "ComBat": "X_combat", "SIMBA": "X_simba", "Scanorama": "X_scanorama",
                      "scVI": "X_scVI", "scANVI": "X_scANVI", "FastMNN": "X_mnn", "BBKNN": "X_pca", "LIGER": "X_iNMF", "SeuratV4_RPCA": "X_pca",
                      "SeuratV5_RPCA": "X_integrated.rpca", "SeuratV5_CCA": "X_integrated.cca", "SeuratV5_Harmony": "X_harmony",
                      "SeuratV5_Unintegrate": "X_pca", "SeuratV5_FastMNN": "X_integrated.mnn", "Cellhint": "X_pca", "Cellhint_harmonize": "X_pca"}
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 5))

    # compute neighbour graph from raw data
    orig_ad = sc.read_h5ad(raw_h5ad)
    sc.pp.normalize_total(orig_ad)
    sc.pp.log1p(orig_ad)
    sc.pp.highly_variable_genes(orig_ad, flavor="seurat_v3", n_top_genes=2000, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(orig_ad, max_value=10)
    sc.tl.pca(orig_ad, svd_solver='arpack')
    sc.pp.neighbors(orig_ad)

    in_ad = sc.read_h5ad(os.path.join(datadir, in_h5ad))
    embedding_key = embedding_keys[method]
    sc.pp.neighbors(in_ad, use_rep=embedding_key)
    scib.me.cluster_optimal_resolution(in_ad, resolutions=get_resolutions(n=10, max=1), cluster_key=cluster_key, label_key=label_key)

    print("Start computing biological conservation metrics for " + method)
    ari = scib.metrics.ari(in_ad, cluster_key=cluster_key, label_key=label_key)
    print("ARI DONE")
    nmi = scib.metrics.nmi(in_ad, cluster_key=cluster_key, label_key=label_key)
    print("NMI DONE")
    sil = scib.metrics.silhouette(adata=in_ad, group_key=cluster_key, embed=embedding_key, metric='euclidean', scale=True)
    print("sil DONE")
    clisi = caculate_lisi(in_ad, batch_key=batch_key, label_key=label_key, cluster_key=cluster_key, embed=embedding_key, type_='clisi', scaled=True)
    print("clisi DONE")

    print("Start computing batch correction metrics for " + method)
    silb = scib.metrics.silhouette_batch(in_ad, batch_key=batch_key, group_key=label_key, embed=embedding_key, metric='euclidean', return_all=True, scale=True, verbose=True)
    print("silb DONE")
    pcr = scib.metrics.pcr_comparison(adata_pre=orig_ad, adata_post=in_ad, covariate=batch_key, embed=embedding_key, n_comps=50, scale=True, verbose=True)
    print("pcr DONE")
    in_ad.obs[label_key] = in_ad.obs[label_key].astype('category')
    gra_knn = scib.metrics.graph_connectivity(in_ad, label_key='map_celltype')
    print("graph DONE")
    ilisi = caculate_lisi(in_ad, batch_key=batch_key, label_key=label_key, cluster_key=cluster_key, embed=embedding_key, type_='ilisi', scaled=True)
    print("ilisi DONE")


    ## save result
    output = pd.DataFrame()
    output.loc['NMI_cluster/label', 'value'] = nmi
    output.loc['ARI_cluster/label', 'value'] = ari
    output.loc['silhouette', 'value'] = sil
    output.loc['iLISI', 'value'] = ilisi
    output.loc['cLISI', 'value'] = clisi
    output.loc['graph_conn', 'value'] = gra_knn
    output.loc['pcr', 'value'] = pcr
    output.loc['silhouette_batch', 'value'] = silb[0]
    output.loc['batch_key', 'value'] = batch_key
    output.loc['cluster_key', 'value'] = cluster_key
    output.loc['integration_method', 'value'] = method
    output.T.to_csv(out_csv)
    # save silb result
    silb[1]['input_h5ad'] = in_h5ad
    silb[1]['unintegrated_h5ad'] = orig_ad
    silb[1]['batch_key'] = batch_key
    silb[1]['cluster_key'] = cluster_key
    silb[1]['integration_method'] = method
    silb[1].to_csv(out_sil)

    in_ad.write(out_h5ad, compression='gzip')
    print("finish batch metrics")







