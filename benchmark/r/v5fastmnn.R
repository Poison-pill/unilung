library(Seurat)
library(BPCells)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(patchwork)
library(anndata)
library(cowplot)
library(SeuratDisk)
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")

indir <- "inte_bench"
outdir <- "inte_bench/inte_data"
figdir <- "inte_bench/fig"
in_rds <- file.path(indir,'raw.rds')


v4_obj <- readRDS(in_rds)

# covert a v4 assay to a v5 assay
seu_obj <- v4_obj
seu_obj[["RNA"]]<-as(object=seu_obj[["RNA"]],Class="Assay5")

# split 
seu_obj[["RNA"]] <- split(seu_obj[["RNA"]], f = seu_obj$batch)

seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize")
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
seu_obj <- ScaleData(seu_obj)
seu_obj <- RunPCA(seu_obj, npcs = 50,verbose = F)

# fastmnn
seu_obj <- IntegrateLayers(
  object = seu_obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
seu_obj <- FindNeighbors(seu_obj, reduction = "integrated.mnn", dims = 1:30)
seu_obj <- RunUMAP(seu_obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")


seu_obj <- JoinLayers(seu_obj)
saveRDS(seu_obj,file=file.path(outdir,'v5fastmnn.rds'))

# covert a v5 assay to a v4 assay
seu_obj2 <- seu_obj
seu_obj2[["RNA"]]<-as(object=seu_obj2[["RNA"]],Class="Assay")
SaveH5Seurat(seu_obj2,
             filename = file.path(outdir,'v5fastmnn.h5seurat'),
             assay = "mnn.reconstructed",
             overwrite = TRUE
)

Convert(file.path(outdir,'v5fastmnn.h5seurat'), "h5ad", overwrite = TRUE)

plot_grid(ncol = 1,
          DimPlot(seu_obj, reduction = "umap.mnn", group.by = 'batch', label = F)+NoAxes(),
          DimPlot(seu_obj, reduction = "umap.mnn", group.by = 'map_ann_level1', label = F)+NoAxes(),
          DimPlot(seu_obj, reduction = "umap.mnn", group.by = 'map_ann_level2', label = F)+NoAxes(),
          DimPlot(seu_obj, reduction = "umap.mnn", label = F)+NoAxes()
)
ggplot2::ggsave(filename = file.path(figdir,'v5mnn_umap.pdf'),width = 10,height = 32)
