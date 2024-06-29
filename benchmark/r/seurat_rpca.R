library(Seurat)
library(patchwork)
library(future)
library(anndata)
library(tidyverse)
library(cowplot)
library(SeuratDisk)
plan("multisession", workers = 10)
options(future.globals.maxSize = 376 * 1024^3) 

indir <- "inte_bench"
outdir <- "inte_bench/inte_data"
figdir <- "inte_bench/fig"
in_rds <- file.path(indir,'raw.rds')
out_h5seurat <- file.path(outdir,'seurat_rpca.h5seurat')
out_h5ad <- file.path(outdir,'seurat_rpca.h5ad')

obj <- readRDS(in_rds)
seurat_list <- SplitObject(obj, split.by = "batch")

# normalize and identify HVG for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = seurat_list)

# run rpca

seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, reduction = "rpca", dims = 1:50)
seu_comb <- IntegrateData(anchors = seurat_anchors, k.weight = 50)

DefaultAssay(seu_comb) <- "integrated"
seu_comb <- ScaleData(seu_comb)
seu_comb <- RunPCA(seu_comb, npcs = 50,verbose = F)
seu_comb <- RunUMAP(seu_comb, reduction = "pca", dims = 1:30)
seu_comb <- FindNeighbors(seu_comb, reduction = "pca", dims = 1:30)
#seu_comb <- FindClusters(resolution = 1, algorithm = 4)# 4 means Leiden
saveRDS(seu_comb,file=file.path(outdir,'seurat_rpca.rds'))

SaveH5Seurat(seu_comb,
             filename = out_h5seurat,
             assay = "integrated",
             overwrite = TRUE
)
Convert(out_h5seurat, out_h5ad, overwrite = F)

plot_grid(ncol = 1,
          DimPlot(seu_comb, reduction = "umap", group.by = 'batch', label = F)+NoAxes(),
          DimPlot(seu_comb, reduction = "umap", group.by = 'map_ann_level1', label = F)+NoAxes(),
          DimPlot(seu_comb, reduction = "umap", group.by = 'map_ann_level2', label = F)+NoAxes(),
          DimPlot(seu_comb, reduction = "umap", label = F)+NoAxes()
)
ggplot2::ggsave(filename = file.path(figdir,'seuratrpca_umap.pdf'),width = 10,height = 32)