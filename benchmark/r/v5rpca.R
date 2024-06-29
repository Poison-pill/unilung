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

# rpca
seu_obj<-readRDS("/public8/lilab/student/fjwu/btit1/inte_bench/v5rpca_seuobj.rds")
seu_obj <- IntegrateLayers(
  object = seu_obj, method = RPCAIntegration, k.weight = 50,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
seu_obj <- FindNeighbors(seu_obj, reduction = "integrated.rpca", dims = 1:30)
seu_obj <- RunUMAP(seu_obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

seu_obj <- JoinLayers(seu_obj)
saveRDS(seu_obj,file=file.path(outdir,'v5rpca.rds'))

# covert a v5 assay to a v4 assay
seu_obj2 <- seu_obj
seu_obj2[["RNA"]]<-as(object=seu_obj2[["RNA"]],Class="Assay")
SaveH5Seurat(seu_obj2,
             filename = file.path(outdir,'v5rpca.h5seurat'),
             assay = "RNA",
             overwrite = TRUE
)

Convert(file.path(outdir,'v5rpca.h5seurat'), "h5ad", overwrite = TRUE)

plot_grid(ncol = 1,
          DimPlot(seu_obj, reduction = "umap.rpca", group.by = 'batch', label = F)+NoAxes(),
          DimPlot(seu_obj, reduction = "umap.rpca", group.by = 'map_ann_level1', label = F)+NoAxes(),
          DimPlot(seu_obj, reduction = "umap.rpca", group.by = 'map_ann_level2', label = F)+NoAxes(),
          DimPlot(seu_obj, reduction = "umap.rpca", label = F)+NoAxes()
)
ggplot2::ggsave(filename = file.path(figdir,'v5rpca_umap.pdf'),width = 10,height = 32)
