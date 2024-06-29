library(anndata)
library(tidyverse)
library(rliger)
library(cowplot)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)


indir <- "inte_bench"
outdir <- "inte_bench/inte_data"
figdir <- "inte_bench/fig"
in_rds <- file.path(indir,'raw.rds')
out_h5seurat <- file.path(outdir,'liger.h5seurat')
out_h5ad <- file.path(outdir,'liger.h5ad')

liger_obj <- readRDS(in_rds)
liger_obj <- NormalizeData(liger_obj, normalization.method = "LogNormalize")
liger_obj <- FindVariableFeatures(liger_obj, selection.method = "vst", nfeatures = 2000)
liger_obj <- ScaleData(liger_obj, split.by = "batch", do.center = FALSE)

# LIGER
liger_obj <- RunOptimizeALS(liger_obj, k = 30, split.by = "batch")
liger_obj <- RunQuantileNorm(liger_obj, split.by = "batch")

liger_obj <- RunUMAP(liger_obj, dims = 1:ncol(liger_obj[["iNMF"]]), reduction = "iNMF", n_neighbors = 15L,  min_dist = 0.3)
liger_obj <- FindNeighbors(liger_obj, reduction = "iNMF")
# liger_obj <- FindClusters(liger_obj, resolution = 1, algorithm = 4)# 4 means Leiden
saveRDS(liger_obj,file=file.path(outdir,'liger.rds'))


# have to convert all factor to character, or when later converting to h5ad, the factors will be numbers
i <- sapply(liger_obj@meta.data, is.factor)
liger_obj@meta.data[i] <- lapply(liger_obj@meta.data[i], as.character)


SaveH5Seurat(liger_obj,
             filename = out_h5seurat,
             overwrite = TRUE
)
# iNMF embedding will be in .obsm['iNMF']
Convert(out_h5seurat, "h5ad", overwrite = TRUE)


plot_grid(ncol = 1,
          DimPlot(liger_obj, reduction = "umap", group.by = 'batch', label = F)+NoAxes(),
          DimPlot(liger_obj, reduction = "umap", group.by = 'map_ann_level1', label = F)+NoAxes(),
          DimPlot(liger_obj, reduction = "umap", group.by = 'map_ann_level2', label = F)+NoAxes(),
          DimPlot(liger_obj, reduction = "umap", label = F)+NoAxes()
)
ggplot2::ggsave(filename = file.path(figdir,'liger_umap.pdf'),width = 10,height = 32)