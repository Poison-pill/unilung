library(bbknnR)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(anndata)
library(cowplot)
library(SeuratDisk)

indir <- "inte_bench"
outdir <- "inte_bench/inte_data"
figdir <- "inte_bench/fig"
in_rds <- file.path(indir,'raw.rds')
out_h5seurat <- file.path(outdir,'bbknn.h5seurat')
out_h5ad <- file.path(outdir,'bbknn.h5ad')

bbknn_obj <- readRDS(in_rds)
bbknn_obj <- NormalizeData(bbknn_obj, normalization.method = "LogNormalize")
bbknn_obj <- FindVariableFeatures(bbknn_obj, selection.method = "vst", nfeatures = 2000)
bbknn_obj <- ScaleData(bbknn_obj, verbose = F)
bbknn_obj <- RunPCA(bbknn_obj, npcs = 30,verbose = F)

# bbknn
bbknn_obj <- RunBBKNN(bbknn_obj, reduction = "pca", batch_key = "batch",
                      run_TSNE = F, run_UMAP = F)
bbknn_obj <- FindNeighbors(bbknn_obj, reduction = "pca", dims = 1:30)
bbknn_obj <- RunUMAP(bbknn_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.bbknn")

saveRDS(bbknn_obj,file=file.path(outdir,'bbknn.rds'))

SaveH5Seurat(bbknn_obj,
             filename = out_h5seurat,
             assay = "RNA",
             overwrite = TRUE
)
Convert(out_h5seurat, out_h5ad, overwrite = F)

plot_grid(ncol = 1,
          DimPlot(bbknn_obj, reduction = "umap.bbknn", group.by = 'batch', label = F)+NoAxes(),
          DimPlot(bbknn_obj, reduction = "umap.bbknn", group.by = 'map_ann_level1', label = F)+NoAxes(),
          DimPlot(bbknn_obj, reduction = "umap.bbknn", group.by = 'map_ann_level2', label = F)+NoAxes(),
          DimPlot(bbknn_obj, reduction = "umap.bbknn", label = F)+NoAxes()
)
ggplot2::ggsave(filename = file.path(figdir,'bbknn_umap.pdf'),width = 10,height = 32)

