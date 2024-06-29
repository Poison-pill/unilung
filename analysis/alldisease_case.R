library(Seurat)
library(tidyverse)
library(cowplot)

dir <- "core_analysis"


myprocess <- function(adata){
  library(reticulate)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  use_condaenv("envs/scvi")
  scanpy <- import("scanpy")
  pd <- import("pandas")
  np <- import("numpy")
  
  ad <- paste0(adata,"_inte.h5ad")
  in_ad <- scanpy$read_h5ad(file.path(dir,ad))
  count <- np$transpose(in_ad$X)
  colnames(count) <- rownames(in_ad$obs)
  rownames(count) <- rownames(in_ad$var)
  mt <- paste0(adata,"_meta.csv")
  meta <- read.csv(file.path(dir,mt),row.names = 1)
  seu <- CreateSeuratObject(counts = count, meta.data = meta)
  
  anvi <- in_ad$obsm[["X_scANVI"]]
  colnames(anvi) <- paste0("scanvi_", 1:10)
  rownames(anvi) <- colnames(seu)
  seu[["scanvi"]] <- CreateDimReducObject(embeddings = anvi, key = "scanvi_", assay = "RNA")
  svraw <- paste0(adata,"_raw.rds")
  saveRDS(seu,file.path(dir,svraw))
  print("=============== Done save raw rds ===============")
  
  sce<-NormalizeData(seu,normalization.method = "LogNormalize",scale.factor = 10000)
  sce<-FindVariableFeatures(sce,reductionselection.method = "vst",nfeatures = 2000)
  sce<-ScaleData(sce)
  sce<-RunPCA(object = sce,npcs = 30,pc.genes=VariableFeatures(sce),verbose = F)
  sce<-FindNeighbors(sce,reduction="scanvi",dims=1:10)
  sce <- RunUMAP(sce, reduction="scanvi", dims=1:10)
  print("=============== Done UMAP ===============")
  
  svpdf <- paste0(adata,"_celltype2.pdf")
  plot_grid(ncol = 2,
            DimPlot(sce,reduction = 'umap',label = F,group.by = "final_celltype")+NoAxes(),
            DimPlot(sce,reduction = 'umap',label = F,group.by = "donor_status")+NoAxes())
  ggplot2::ggsave(filename = file.path(dir,svpdf),width = 15,height = 7)
  svprocess <- paste0(adata,"_proc.rds")
  saveRDS(sce,file.path(dir,svprocess))
  return(sce)
}


myprocess("cd8t")
myprocess("cd4t")
myprocess("macro2")



#define a pseudo_cell func
pseudo_cell <- function(object,organize="orig.ident",cell_n=10,method = 'average',chunk_size=100000){
  
  library(Seurat);library(dplyr)
  
  average_expr <- function(myobject,categorie,chunk.no){
    if( ncol(myobject) >1){
      ##generate category.matrix
      category <- myobject@meta.data[,'pseudo_label']
      if (length(table(category)) >1 ) {
        category.matrix <- Matrix::sparse.model.matrix( ~ -1 + category )
        colnames(category.matrix) <- gsub("category","",colnames(category.matrix))
        colnames(category.matrix) <- paste0(colnames(category.matrix),".",chunk.no)
      }else{
        category.matrix <- matrix(data = 1, nrow = ncol(myobject),dimnames = list(Cells(myobject), paste0(categorie,"_N.1.",chunk.no)))
      }
      ##aggregate cell by average expr
      if (method == 'average') {category.matrix <- sweep(category.matrix, 2,Matrix::colSums(category.matrix), FUN = "/") }
      data.use <- GetAssayData(myobject, slot = "counts")
      aggreate.counts = data.use %*% category.matrix
    }else{
      aggreate.counts <- myobject[["RNA"]]@counts
      colnames(aggreate.counts) <- paste0(categorie,"_N.1.",chunk.no);
    }
    return(aggreate.counts)
  }
  
  DefaultAssay(object) <- "RNA"
  categories=names(table(object@meta.data[,organize]))
  data.list <- list()
  
  for (idx in 1:length(categories)) {
    categorie = categories[idx]
    subObj <- object[, object@meta.data[, organize] == categorie]
    
    if(ncol(subObj) > chunk_size){
      subObj@meta.data <- subObj@meta.data %>% mutate(cellN=colnames(subObj)) %>% group_by(get(organize)) %>% 
        mutate( chunk_label = rep(1:ceiling(n()/chunk_size), each=chunk_size, length.out=n()) ) %>% 
        tibble::column_to_rownames('cellN') %>% data.frame()
      data.list.tmp <- list() 
      for (chunk in 1:length(table(subObj$chunk_label))) {
        subObj.chunk <- subset(subObj,subset=chunk_label==chunk)
        subObj.chunk@meta.data <- subObj.chunk@meta.data %>% mutate(cellN=colnames(subObj.chunk)) %>% group_by(get(organize)) %>% 
          mutate( bins.no = rep(1:ceiling(n()/cell_n), each=cell_n, length.out=n()) ) %>% 
          mutate(pseudo_label=paste0(get(organize),"_N.",bins.no)) %>% tibble::column_to_rownames('cellN') %>% data.frame()
        agg.counts <- average_expr(subObj.chunk,categorie=categorie,chunk.no = chunk)
        data.list.tmp[[chunk]] <- CreateSeuratObject(counts = agg.counts)
        data.list.tmp[[chunk]]$orig.ident <- categorie
      }
      toRet.chunk = data.list.tmp[[chunk]]
      if (length(data.list.tmp) > 1) {toRet.chunk <- merge(data.list.tmp[[1]],data.list.tmp[2:length(data.list.tmp)])}
      data.list[[idx]] <- DietSeurat(toRet.chunk)
    }else{
      subObj@meta.data <- subObj@meta.data %>% mutate(cellN=colnames(subObj)) %>% group_by(get(organize)) %>% 
        mutate( bins.no = rep(1:ceiling(n()/cell_n), each=cell_n, length.out=n()) ) %>% 
        mutate(pseudo_label=paste0(get(organize),"_N.",bins.no)) %>% tibble::column_to_rownames('cellN') %>% data.frame()
      agg.counts <- average_expr(subObj,categorie=categorie,chunk.no = 1)
      data.list[[idx]] <- CreateSeuratObject(counts = agg.counts)
      data.list[[idx]]$orig.ident <- categorie
    }
  }
  toRet =  data.list[[idx]]
  if (length(data.list) > 1) {toRet <- merge(data.list[[1]],data.list[2:length(data.list)])}
  toRet <- DietSeurat(toRet)
  return(toRet)
}


build_pseudo <- function(dat,organize = NULL,cell_n = NULL){
  ad <- paste0(dat,"_proc.rds")
  seu <- readRDS(file.path(dir,ad))
  sc <- pseudo_cell(seu,organize = organize,cell_n = cell_n)
  saveRDS(sc,file.path(dir,paste0(dat,"_pseudo.rds")))
  write.csv(sc@assays$RNA@counts, file.path(dir,paste0(dat,"_pseudo.csv")), row.names = T)
  return(sc)
}

cd4t_small <- build_pseudo("cd4t",organize="donor_status",cell_n=10)
cd8t_small <- build_pseudo("cd8t",organize="donor_status",cell_n=10)
macro2_small <- build_pseudo("macro2",organize="donor_status",cell_n=10)






## run scFEA in bash
python core_analysis/scFEA/src/scFEA.py --data_dir core_analysis/scFEA/data --input_dir core_analysis \
--test_file cd8t_pseudo.csv \
--moduleGene_file module_gene_m168.csv \
--stoichiometry_matrix cmMat_c70_m168.csv \
--res_dir core_analysis/metabolism/cd8t \
--output_flux_file core_analysis/metabolism/cd8t/flux.csv \
--output_balance_file core_analysis/metabolism/cd8t/balance.csv \
--sc_imputation True

python core_analysis/scFEA/src/scFEA.py --data_dir core_analysis/scFEA/data --input_dir core_analysis \
--test_file cd4t_pseudo.csv \
--moduleGene_file module_gene_m168.csv \
--stoichiometry_matrix cmMat_c70_m168.csv \
--res_dir core_analysis/metabolism/cd4t \
--output_flux_file core_analysis/metabolism/cd4t/flux.csv \
--output_balance_file core_analysis/metabolism/cd4t/balance.csv \
--sc_imputation True

python core_analysis/scFEA/src/scFEA.py --data_dir core_analysis/scFEA/data --input_dir core_analysis \
--test_file macro2_pseudo.csv \
--moduleGene_file module_gene_m168.csv \
--stoichiometry_matrix cmMat_c70_m168.csv \
--res_dir core_analysis/metabolism/macro2 \
--output_flux_file core_analysis/metabolism/macro2/flux.csv \
--output_balance_file core_analysis/metabolism/macro2/balance.csv \
--sc_imputation True







# add flux as a new assay
scEFA_flux <- function(name){
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  dir <- "core_analysis"
  
  rds <- paste0(name,"_pseudo.rds")
  obj <- readRDS(file.path(dir,rds))
  rb.gene <- rownames(obj[grep("^RP[SL]",rownames(obj))])
  mt.gene <- rownames(obj[grep("^MT-",rownames(obj))])
  obj <- subset(x = obj,feature=setdiff(rownames(obj),rb.gene))
  obj <- subset(x = obj,feature=setdiff(rownames(obj),mt.gene))
  obj <- NormalizeData(obj,normalization.method = "LogNormalize",scale.factor = 10000)
  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 2000) %>% 
    ScaleData() %>% RunPCA(npcs = 30,pc.genes=VariableFeatures(obj),verbose = F)
  obj <- FindNeighbors(obj,dims=1:30)
  obj <- RunUMAP(obj, dims=1:30)
  
  flux <- paste0("metabolism/",name,"/flux.csv")
  predFlux <- read.csv(file.path(dir,flux), header = T, row.names = 1)
  Flux <- data.matrix(predFlux) %>% t()
  obj[["FLUX"]] <- CreateAssayObject(counts = Flux)
  DefaultAssay(obj) <- 'FLUX'
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = F)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 30, reduction.name = 'pca.flux', verbose = F)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = F)
  #obj <- FindClusters(obj, resolution = 0.5, verbose = F)
  obj <- RunUMAP(obj, dims = 1:30, assay = 'FLUX', reduction.name = "umap.flux", verbose = F)
  saveRDS(obj,file.path(dir,paste0(name,'_metabolism.rds')))
  plot_grid(ncol = 2,
            DimPlot(obj,reduction = "umap",group.by = "orig.ident") + ggtitle('umap of Gene') + NoAxes(),
            DimPlot(obj,reduction = 'umap.flux',group.by = "orig.ident",label = F) + ggtitle('umap of Flux') + NoAxes())
  ggplot2::ggsave(filename = file.path(dir,paste0(name,'_metabolism.pdf')),width = 15,height = 7)
}

# gene deg caculate
save_deg <- function(name){
  library(Seurat)
  library(tidyverse)
  dir <- "core_analysis"
  
  obj <- readRDS(file.path(dir,paste0(name,"_metabolism.rds")))
  if("SeuratProject" %in% obj@active.ident){
    obj@active.ident <- as.factor(obj$orig.ident)
  }
  gene.markers <- FindAllMarkers(object = obj,assay = 'RNA',only.pos = F,min.pct = 0.1,logfc.threshold = 0.25,test.use = "wilcox")
  gene_deg <- gene.markers %>% group_by(cluster) %>% top_n(n = 2, wt = abs(avg_log2FC))
  save(obj,gene.markers,gene_deg,file = file.path(dir,paste0(name,"_deg.RData")))
}

## deg volcanoplot
multiVolcanoPlot <- function(ntop,celltype=NULL){
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  dir <- "core_analysis"
  load(file.path(dir,paste0(celltype,"_deg.RData")))
  df <- gene.markers
  
  df$label <- ifelse(df$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
  topsig <- df %>% group_by(cluster) %>% distinct(gene, .keep_all = T) %>% top_n(ntop, abs(avg_log2FC))
  
  #非Top
  df$size <- case_when(!(df$gene %in% topsig$gene)~ 1, df$gene %in% topsig$gene ~ 2)
  dt <- filter(df,size==1)
  
  # 背景柱状图数据准备
  dbar <- df %>% group_by(cluster) %>% summarise_all(list(min = min, max = max)) %>% 
    select(cluster, avg_log2FC_min, avg_log2FC_max, label_min) %>% rename(label = label_min)
  

  ## plot
  ggplot()+
    geom_col(data = dbar,  # 绘制负向背景柱状图
             mapping = aes(x = cluster,y = avg_log2FC_min),
             fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
    geom_col(data = dbar, # 绘制正向背景柱状图
             mapping = aes(x = cluster,y = avg_log2FC_max),
             fill = "#dcdcdc",alpha = 0.6, width = 0.7)+
    geom_jitter(data = dt, # 除 top10外的基因散点图
                aes(x = cluster, y = avg_log2FC, color = label),
                size = 0.85,
                width =0.4)+
    geom_jitter(data = topsig, # top10
                aes(x = cluster, y = avg_log2FC, color = label),
                size = 1,
                width =0.4)+
    geom_tile(data = topsig, # 绘制中心分组
              aes(x = cluster,
                  y = 0,
                  fill = cluster),
              height=0.4,
              width=1,
              color = "white",
              alpha = 0.5,
              show.legend = F)+
    ggsci::scale_fill_simpsons() + # 自定义颜色
    scale_color_manual(name=NULL,values = c("#A73030FF","black"))+
    geom_text_repel(data = topsig,  # 筛选你想要标记的基因
                    aes(x = cluster, y = avg_log2FC, label = gene),
                    size = 2.5, 
                    max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
                    color = 'black',
                    force = 1.2,
                    arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last"))+
    labs(x="donor_status", y="Average log2FC")+
    geom_text(data=topsig, # 绘制中心分组标记图文本注释
              aes(x=cluster, y=0, label=cluster),
              size = 2.5,
              color ="white")+
    theme_minimal() + 
    theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
          axis.line.y = element_line(color = "black",size = 1.2),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.direction = "vertical",
          legend.justification = c(1,0),
          legend.text = element_text(size = 13))
  ggplot2::ggsave(filename = file.path(dir,paste0(celltype,'_deg.pdf')),width = 10,height = 6)
}
multiVolcanoPlot(ntop = 5,celltype = "cd4t")
multiVolcanoPlot(ntop = 5,celltype = "cd8t")
multiVolcanoPlot(ntop = 5,celltype = "macro2")



## metabolism module deg
metabolism.markers <- FindAllMarkers(object = obj,assay = 'FLUX',only.pos = F,min.pct = 0.1,logfc.threshold = 0.25,test.use = "wilcox")
metabolism_module <- read.csv("Human_M168_information.symbols.csv")
module <- data.frame(id=metabolism_module$X,
                     name=paste(metabolism_module$Compound_IN_name,metabolism_module$Compound_OUT_name,sep = "_"))
module$id <- gsub("_","-",module$id)
metabolism.markers$id <- metabolism.markers$gene
metabolism.markers <- left_join(metabolism.markers,module,by="id")

DotPlot(obj,unique(metabolism.markers$gene),assay = "FLUX")+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色


## Feature plot for metabolism module
## M-22 in cd4t from luad; M-114 in cd4t from asthma 
# (one module, celltype and status specific)
FeaturePlot(obj,"M-22",reduction = "umap.flux",cols = c("grey", "#A73030FF")) 
# Glycine to Glycine-OUT
FeaturePlot(obj,"M-114",reduction = "umap.flux",cols = c("grey", "#A73030FF")) 
# (E,E)-Farnesyl-PP to Geranylgeranyl-PP

## M-156 in cd4t/cd8t/macro2 from PF
# (one module, across celltype, status specific)
FeaturePlot(obj,"M-156",reduction = "umap.flux",cols = c("grey", "#A73030FF")) 
# CDP to Cytidine

## M-122 in all cancer (cd4t from sclc; macro2 from luad/lusc) 
# (module specific, across celltype and status)
FeaturePlot(obj,"M-122",reduction = "umap.flux",cols = c("grey", "#A73030FF")) 
# (GlcNAc)4 (Man)3 (Asn)1 to (Gal)2 (GlcNAc)4 (LFuc)1 (Man)3 (Neu5Ac)2 (Asn)1

## M-113 in ILD and PF (cd4t from ild; cd8t from pf) 
# (module specific, across celltype and status)
FeaturePlot(obj,"M-113",reduction = "umap.flux",cols = c("grey", "#A73030FF")) 
# Acetyl-CoA to (E,E)-Farnesyl-PP

## M-113, M-167 and M-168 in cd4t from ILD ()
# (across module, celltype and status specific)
FeaturePlot(obj,c("M-113","M-167","M168",reduction = "umap.flux",cols = c("grey", "#A73030FF"))
# Acetyl-CoA to (E,E)-Farnesyl-PP
# (E,E)-Farnesyl-PP to Cholesterol
# Cholesterol to Chenodeoxycholate


## plot in cd4t
load("cd4t_deg.RData")
FeaturePlot(obj,c("M-22","M-114",
                  "M-156",
                  "M-122",
                  "M-113",
                  "M-113","M-167","M-168"),
            reduction = "umap.flux",cols = c("grey", "#A73030FF"),ncol = 2)
ggplot2::ggsave(filename = "cd4t_meta_feat.pdf",width = 8,height = 16)

load("cd8t_deg.RData")
FeaturePlot(obj,c("M-156",
                  "M-113"),
            reduction = "umap.flux",cols = c("grey", "#A73030FF"),ncol = 2)
ggplot2::ggsave(filename = "cd8t_meta_feat.pdf",width = 8,height = 4)

load("macro2_deg.RData")
FeaturePlot(obj,c("M-156",
                  "M-122"),
            reduction = "umap.flux",cols = c("grey", "#A73030FF"),ncol = 2)
ggplot2::ggsave(filename = "macro2_meta_feat.pdf",width = 8,height = 4)


## build metabolism deg count for vln plot in cnsknowall
mat <- as.matrix(obj@assays$FLUX@counts) %>% t() %>% as.data.frame()
deg_mat <- mat[,which(colnames(mat) %in% metabolism.markers$gene)]
deg_mat <- cbind(status=obj@meta.data$orig.ident,deg_mat)
write.csv(deg_mat,"deg_mat.csv",row.names = F)




## gene deg enrichment analysis
myenrich <- function(celltype,width=NULL,height=NULL){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
  dir <- "core_analysis"
  load(file.path(dir,paste0(celltype,"_deg.RData")))
  
  go_result <- rbind()
  for (i in unique(gene.markers$cluster)) {
    deg_up <- subset(gene.markers, avg_log2FC>1 & cluster==i)
    go_up <- enrichGO(gene = deg_up$gene,
                      OrgDb = 'org.Hs.eg.db',
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01)
    if(nrow(go_up)>=2){
      up <- data.frame(go_up)[1:2,]
      up$type <- "up"
    } else {
      up <- rbind()
    }
    deg_down <- subset(gene.markers, avg_log2FC<-1 & cluster==i)
    go_down <- enrichGO(gene = deg_down$gene,
                      OrgDb = 'org.Hs.eg.db',
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01) 
    if(nrow(go_down)>=2){
      down <- data.frame(go_down)[1:2,]
      down$type <- "down"
    } else {
      down <- rbind()
    }
    all <- rbind(up,down)
    all$group=i
    all <- all[!duplicated(all$Description),]
    go_result <- rbind(go_result,all)
  }
  go_result$pvalue <- -log10(go_result$pvalue)
  ggplot(go_result,aes(pvalue,Description))+
    geom_bar(aes(y=reorder(Description,pvalue),x=pvalue,fill=go_result$type)
             ,stat='identity',width=0.5)+
    scale_fill_manual(values = c(up="#7DC5A0",down="#D58890"))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    facet_grid(group~., scale="free", space = "free_y")
  ggplot2::ggsave(filename = file.path(dir,paste0(celltype,"_enrich.pdf")),
                  width = width,height = height)
}


myenrich("cd4t",width = 8,height = 10)
myenrich("cd8t",width = 9,height = 10)
myenrich("macro2",width = 8,height = 12)


