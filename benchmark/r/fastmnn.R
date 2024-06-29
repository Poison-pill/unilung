library(anndata)
library(tidyverse)
library(cowplot)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)


indir <- "inte_bench"
outdir <- "inte_bench/inte_data"
figdir <- "inte_bench/fig"
in_rds <- file.path(indir,'raw.rds')
out_h5seurat <- file.path(outdir,'fastmnn.h5seurat')
out_h5ad <- file.path(outdir,'fastmnn.h5ad')

mnn_obj <- readRDS(in_rds)
mnn_obj <- NormalizeData(mnn_obj, normalization.method = "LogNormalize")
mnn_obj <- FindVariableFeatures(mnn_obj, selection.method = "vst", nfeatures = 2000)


## fastMNN
mnn_obj <- RunFastMNN(object.list = SplitObject(mnn_obj, split.by = "batch"))
mnn_obj <- RunUMAP(mnn_obj, reduction = "mnn", dims = 1:30, n_neighbors = 15L,  min_dist = 0.3)
mnn_obj <- FindNeighbors(mnn_obj, reduction = "mnn", dims = 1:30)
#mnn_obj <- FindClusters(mnn_obj, resolution = 1, algorithm = 4)# 4 means Leiden
saveRDS(mnn_obj,file=file.path(outdir,'fastmnn.rds'))

# have to convert all factor to character, or when later converting to h5ad, the factors will be numbers
i <- sapply(mnn_obj@meta.data, is.factor)
mnn_obj@meta.data[i] <- lapply(mnn_obj@meta.data[i], as.character)

SaveH5Seurat(mnn_obj,
             filename = out_h5seurat,
             assay = "mnn.reconstructed",
             overwrite = TRUE
)

# fastMNN embedding will be in adata.X
Convert(out_h5seurat, out_h5ad, overwrite = F)

plot_grid(ncol = 1,
          DimPlot(mnn_obj, reduction = "umap", group.by = 'batch', label = F)+NoAxes(),
          DimPlot(mnn_obj, reduction = "umap", group.by = 'map_ann_level1', label = F)+NoAxes(),
          DimPlot(mnn_obj, reduction = "umap", group.by = 'map_ann_level2', label = F)+NoAxes(),
          DimPlot(mnn_obj, reduction = "umap", label = F)+NoAxes()
)
ggplot2::ggsave(filename = file.path(figdir,'fastmnn_umap.pdf'),width = 10,height = 32)