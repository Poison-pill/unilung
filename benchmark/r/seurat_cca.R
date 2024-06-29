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
out_h5seurat <- file.path(outdir,'seurat_cca.h5seurat')
out_h5ad <- file.path(outdir,'seurat_cca.h5ad')

obj <- readRDS(in_rds)

seurat_list <- SplitObject(obj, split.by = "batch")

# normalize and identify HVG for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seurat_list)
save(features,seurat_list,file="/public8/lilab/student/fjwu/btit1/inte_bench/cca_rpca_feature.RData")

# find anchors
load("/public8/lilab/student/fjwu/btit1/inte_bench/cca_rpca_feature.RData")
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, dims = 1:50)

# integrate by anchors
seu_comb <- IntegrateData(anchorset = seurat_anchors)
save(seu_comb,seurat_anchors,file="/public8/lilab/student/fjwu/btit1/inte_bench/seurat_cca.RData")
save(seu_comb,file="/public8/lilab/student/fjwu/btit1/inte_bench/inte_cca.RData")

load("/public8/lilab/student/fjwu/btit1/inte_bench/inte_cca.RData")
DefaultAssay(seurat_comb) <- "integrated"
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
ggplot2::ggsave(filename = file.path(figdir,'seuratcca_umap.pdf'),width = 10,height = 32)