library(Seurat)
library(tidyverse)
library(cowplot)
library(CytoTRACE)
dir <- "case_ca"
setwd(dir)


sce <- readRDS(file.path(dir,'epi_new.rds'))
sce2 <- SketchData(object = sce,ncells = 100000,method = "LeverageScore",sketched.assay = "bpcell")
sce2 <- sce2@assays$bpcell@counts
sce_sub <- subset(sce,cells=colnames(sce2))

mat <- as.matrix(sce_sub@assays$RNA@counts)
results <- CytoTRACE(mat, ncores = 8)

emb <- sce_sub@reductions[["umap"]]@cell.embeddings

celltype <- as.character(sce_sub$annotation)
names(celltype) <- rownames(sce_sub@meta.data)
plotCytoTRACE(results, phenotype = celltype,emb = emb,outputDir = paste0(dir,"/cytotrace/anno_"))

atlas <- as.character(sce_sub$sub_atlas)
names(atlas) <- rownames(sce_sub@meta.data)
plotCytoTRACE(results, phenotype = atlas, emb = emb, outputDir = paste0(dir,"/cytotrace/atlas_"))

stage <- as.character(sce_sub$stage)
names(stage) <- rownames(sce_sub@meta.data)
plotCytoTRACE(results, phenotype = stage, emb = emb, outputDir = paste0(dir,"/cytotrace/stage_"))





## monocle2
run_monocle <- function(x){
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(monocle)
  dir <- "case_ca"
  options(Seurat.object.assay.version = "v3")
  
  sce <- readRDS(file.path(dir,'epi_new.rds'))
  small <- sce[,-which(sce$annotation=="Ciliated cell")]
  small <- small[,small$sub_atlas==x]
  
  sce2 <- CreateSeuratObject(counts = small@assays$RNA@counts,meta.data = small@meta.data,class = 'Assay')
  sce2[["scanvi"]] <- small[["scanvi"]]
  sce2 <- NormalizeData(sce2,normalization.method = "LogNormalize",scale.factor = 10000)
  sce2 <- FindVariableFeatures(sce2,reductionselection.method = "vst",nfeatures = 2000)
  sce2 <- ScaleData(sce2)
  sce2 <- RunPCA(object = sce2,npcs = 30,pc.genes=VariableFeatures(sce2),verbose = F)
  sce2 <- FindNeighbors(sce2,reduction="scanvi",dims=1:10) %>% RunUMAP(reduction="scanvi", dims=1:10)
  
  sce3 <- SketchData(object = sce2,ncells = 40000,method = "LeverageScore",sketched.assay = "bpcell")
  sce4 <- sce3@assays$bpcell@counts
  sce_sub <- subset(sce,cells=colnames(sce4))
  
  mtx <- as(as.matrix(sce_sub@assays$RNA@counts),'sparseMatrix')
  pdat <- sce_sub@meta.data
  
  #add celltype annotation
  pdat$annotation <- sce_sub$annotation
  pdat$stage <- sce_sub$stage
  fdat <- data.frame(gene_short_name=row.names(sce_sub),row.names = row.names(sce_sub))
  pd <- new("AnnotatedDataFrame",data = pdat)
  fd <- new("AnnotatedDataFrame",data = fdat)
  cds <- newCellDataSet(mtx,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds,min_expr = 0.1)  
  
  #choose DEG for sort
  expre_genes <- VariableFeatures(sce_sub)
  diff <- differentialGeneTest(cds[expre_genes,],fullModelFormulaStr = "~annotation",cores = 4) %>% subset(qval<0.01)
  ordergene <- diff[order(diff$qval,decreasing = F),] %>% rownames
  
  cds <- setOrderingFilter(cds,ordergene)
  cds <- reduceDimension(cds,max_components = 2,method='DDRTree')
  cds <- orderCells(cds)
  cds_stage <- cds
  
  print("=============== monocle plot ===============")
  library(viridis)
  library(ggsci)
  plot_grid(ncol = 3,
            plot_cell_trajectory(cds,color_by = "Pseudotime",size=1,show_backbone = F,show_branch_points = T)+NoAxes()+scale_color_viridis(option="magma"),
            plot_cell_trajectory(cds,color_by = "annotation",size=1,show_backbone = F,show_branch_points = T)+NoAxes()+scale_color_jama(),
            plot_cell_trajectory(cds,color_by = "stage",size=1,show_backbone = F,show_branch_points = T)+NoAxes()+scale_color_nejm())
  ggplot2::ggsave(filename = file.path(dir,paste0("monocle/",x,"_traj.pdf")),width = 20,height = 7)
  save(sce_sub,cds,file = file.path(dir,paste0("monocle/",x,"_cds.RData")))
}

run_monocle("LUAD")
run_monocle("LUSC")
run_monocle("SCLC")

# plot cytotrace to monocle
run_cyto <- function(x){
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(monocle)
  dir <- "case_ca"
  
  load(file.path(dir,paste0("monocle/",x,"_cds.RData")))
  
  print("=============== Start cytotrace ===============")
  library(CytoTRACE)
  mat <- as.matrix(sce_sub@assays$RNA@counts)
  results <- CytoTRACE(mat, ncores = 8)
  
  monocle_meta <- data.frame(t(cds@reducedDimS))
  colnames(monocle_meta) <- c("C1", "C2")
  emb_monocle <- monocle_meta[,1:2]
  
  print("=============== Cytotrace plot ===============")
  celltype <- as.character(sce_sub$annotation)
  names(celltype) <- rownames(sce_sub@meta.data)
  plotCytoTRACE(results, phenotype = celltype,emb = emb_monocle,outputDir = paste0(dir,"/monocle/",x,"_anno_"))
  
  stage <- as.character(sce_sub$stage)
  names(stage) <- rownames(sce_sub@meta.data)
  plotCytoTRACE(results, phenotype = stage, emb = emb_monocle, outputDir = paste0(dir,"/monocle/",x,"_stage_"))
}

run_cyto("LUAD")
run_cyto("LUSC")
run_cyto("SCLC")

## multiple cancer cell trajectory
library(Seurat)
library(tidyverse)
library(cowplot)
library(monocle)
library(ggsci)
dir <- "case_ca"
options(Seurat.object.assay.version = "v3")

sce <- readRDS(file.path(dir,'epi_new.rds'))
small <- sce[,sce$annotation=="Malignant cell"]

sce2 <- CreateSeuratObject(counts = small@assays$RNA@counts,meta.data = small@meta.data,class = 'Assay')
sce2[["scanvi"]] <- small[["scanvi"]]
sce2 <- NormalizeData(sce2,normalization.method = "LogNormalize",scale.factor = 10000)
sce2 <- FindVariableFeatures(sce2,reductionselection.method = "vst",nfeatures = 2000)
sce2 <- ScaleData(sce2)
sce2 <- RunPCA(object = sce2,npcs = 30,pc.genes=VariableFeatures(sce2),verbose = F)
sce2 <- FindNeighbors(sce2,reduction="scanvi",dims=1:10) %>% RunUMAP(reduction="scanvi", dims=1:10) %>% FindClusters(resolution = 0.01)
plot_grid(ncol=2,
          DimPlot(sce2,reduction = "umap",cols = pal_lancet("lanonc")(9),
                  label = F,group.by = "seurat_clusters")+ggtitle("Cluster"),
          DimPlot(sce2,reduction = "umap",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                  label = F,group.by = "sub_atlas")+ggtitle("Cancer type"))
ggplot2::ggsave(filename = file.path(dir,"monocle/multica_umap.pdf"),width = 15,height = 7)
saveRDS(sce2,file = file.path(dir,"monocle/multica_umap.rds"))

sce3 <- SketchData(object = sce2,ncells = 30000,method = "LeverageScore",sketched.assay = "bpcell")
sce4 <- sce3@assays$bpcell@counts
sce_sub <- subset(sce,cells=colnames(sce4))

mtx <- as(as.matrix(sce_sub@assays$RNA@counts),'sparseMatrix')
pdat <- sce_sub@meta.data

#add celltype annotation
pdat$sub_atlas <- sce_sub$sub_atlas
pdat$stage <- sce_sub$stage
fdat <- data.frame(gene_short_name=row.names(sce_sub),row.names = row.names(sce_sub))
pd <- new("AnnotatedDataFrame",data = pdat)
fd <- new("AnnotatedDataFrame",data = fdat)
cds <- newCellDataSet(mtx,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)  

#choose DEG for sort
expre_genes <- VariableFeatures(sce_sub)
diff <- differentialGeneTest(cds[expre_genes,],fullModelFormulaStr = "~annotation",cores = 4) %>% subset(qval<0.01)
ordergene <- diff[order(diff$qval,decreasing = F),] %>% rownames

cds <- setOrderingFilter(cds,ordergene)
cds <- reduceDimension(cds,max_components = 2,method='DDRTree')
cds <- orderCells(cds)

print("=============== monocle plot ===============")
library(viridis)
library(ggsci)
library(cowplot)
plot_grid(ncol = 4,
          plot_cell_trajectory(cds,color_by = "Pseudotime",size=1,show_backbone = F,show_branch_points = F)+NoAxes()+scale_color_viridis(option="magma"),
          plot_cell_trajectory(cds,color_by = "sub_atlas",size=1,show_backbone = F,show_branch_points = F)+NoAxes()+scale_color_jama(),
          plot_cell_trajectory(cds,color_by = "State",size=1,show_backbone = F,show_branch_points = F)+NoAxes()+scale_color_aaas(),
          plot_cell_trajectory(cds,color_by = "stage",size=1,show_backbone = F,show_branch_points = F)+NoAxes()+scale_color_nejm())
ggplot2::ggsave(filename = file.path(dir,"monocle/multica_traj.pdf"),width = 25,height = 7)
save(sce_sub,cds,file = file.path(dir,"monocle/multica_cds.RData"))

run_cyto_multica <- function(x){
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(monocle)
  dir <- "case_ca"
  
  load(file.path(dir,paste0("monocle/",x,"_cds.RData")))
  
  print("=============== Start cytotrace ===============")
  library(CytoTRACE)
  mat <- as.matrix(sce_sub@assays$RNA@counts)
  results <- CytoTRACE(mat, ncores = 8)
  
  monocle_meta <- data.frame(t(cds@reducedDimS))
  colnames(monocle_meta) <- c("C1", "C2")
  emb_monocle <- monocle_meta[,1:2]
  
  print("=============== Cytotrace plot ===============")
  sub_atlas <- as.character(sce_sub$sub_atlas)
  names(sub_atlas) <- rownames(sce_sub@meta.data)
  plotCytoTRACE(results, phenotype = sub_atlas,emb = emb_monocle,outputDir = paste0(dir,"/monocle/",x,"_subca_"))
  
  stage <- as.character(sce_sub$stage)
  names(stage) <- rownames(sce_sub@meta.data)
  plotCytoTRACE(results, phenotype = stage, emb = emb_monocle, outputDir = paste0(dir,"/monocle/",x,"_stage_"))
}
run_cyto_multica("multica")


## monocle version in conda acid
# point 1 analysis for SCLC
library(Seurat)
library(tidyverse)
library(cowplot)
library(monocle)
library(ggsci)
dir <- "case_ca"
load(file.path(dir,"monocle/multica_cds.RData"))

beam_res=BEAM(cds,branch_point = 1,cores = 8)
beam_res=beam_res[,c("gene_short_name","pval","qval")]

tmp=plot_genes_branched_heatmap(cds[row.names(subset(beam_res,qval<1e-4)),],
                                branch_point = 1,
                                num_clusters = 4,
                                cores = 4,
                                branch_labels = c("Cell fate 1", "Cell fate 2"),
                                hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                branch_colors = c("#979797", "#F05662", "#7990C8"),
                                use_gene_short_name = T,
                                show_rownames = F,
                                return_heatmap = T 
)
save(cds,beam_res,tmp,file = file.path(dir,paste0("monocle/point1.RData")))
write.csv(beam_res,file = file.path(dir,paste0("monocle/point1_deg.csv")),row.names=T)
ggsave(tmp$ph_res,filename = file.path(dir,paste0("monocle/point1_heatmap.pdf")),width = 5,height = 6)


library(clusterProfiler)
library(org.Hs.eg.db)


gene_group <- tmp$annotation_row
gene_group$gene <- rownames(gene_group)

fa <- lapply (unique(gene_group$Cluster),function(x){
  small_gene_group <- filter(gene_group,gene_group$Cluster==x)
  go <- enrichGO(gene = unique(small_gene_group$gene),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'SYMBOL',ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.1,
                 readable = TRUE)
  if (dim(go@result)[1] != 0) {
    tt <- go@result
    }
  }
)




library(Seurat)
library(tidyverse)
library(cowplot)
library(monocle)
library(ggsci)

load(file.path(dir,"monocle/multica_cds.RData"))
sce <- readRDS(file.path(dir,"monocle/multica_umap.rds"))

subcds <- cds[,cds$State %in% c(2,3,4)]
sce2 <- sce[,colnames(subcds)]
sce2$State <- subcds$State
sce2 <- FindNeighbors(sce2,reduction="scanvi",dims=1:10)
sce2 <- RunUMAP(sce2, reduction="scanvi", dims=1:10) %>% FindClusters(resolution = 0.02)
sc1 <- sce2[,(sce2$seurat_clusters %in% c("2","3","4","5","6","8") & sce2$sub_atlas=="SCLC")]
sc2 <- sce2[,(sce2$seurat_clusters=="0" & sce2$sub_atlas=="SCLC")]
group <- data.frame(cellid=colnames(sce2),cluster=sce2$seurat_clusters)
group[colnames(sc1),"group"] <- "SCLC-1"
group[colnames(sc2),"group"] <- "SCLC-2"
group[-which(group$cellid %in% c(colnames(sc1),colnames(sc2))),"group"] <- "NSCLC"
sce2$group <- group$group


sce3 <- sce
sce3[["bpcell"]] <- CreateAssayObject(counts = sce2@assays$RNA@counts)
DefaultAssay(sce3) <- "bpcell"
sce3 <- NormalizeData(sce3,normalization.method = "LogNormalize",scale.factor = 10000)
VariableFeatures(sce3)<-VariableFeatures(sce3)
sce3 <- ScaleData(sce3)
sce3 <- RunPCA(object = sce3,npcs = 30,pc.genes=VariableFeatures(sce3),verbose = F)
sce3[["scanvinew"]] <- CreateDimReducObject(embeddings = sce2@reductions[["scanvi"]]@cell.embeddings, key = "scanvinew_", assay = "bpcell")
sce3 <- FindNeighbors(sce3,reduction="scanvinew",dims=1:10)
sce3 <- RunUMAP(sce3, reduction="scanvinew", return.model = TRUE, reduction.name="umapnew",dims=1:10) %>% FindClusters(resolution = 0.02)

sce3 <- ProjectData(object = sce3,
                    assay = "RNA",
                    full.reduction = "pca.full",
                    sketched.assay = "bpcell",
                    sketched.reduction = "pca",
                    umap.model = "umapnew",
                    dims = 1:30,
                    refdata = list(cluster_full = "seurat_clusters"))
DefaultAssay(sce3) <- "RNA"
plot_grid(ncol=2,
          DimPlot(sce3,reduction = "full.umapnew",label = F,group.by = "cluster_full")+
            ggtitle("Cluster")+scale_color_igv(),
          DimPlot(sce3,reduction = "full.umapnew",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                  label = F,group.by = "sub_atlas")+ggtitle("Cancer type"))
ggplot2::ggsave(filename = "kk_umap.pdf",width = 15,height = 7)
ggplot2::ggsave(filename = file.path(dir,"monocle/point_umap.pdf"),width = 15,height = 7)


sc1 <- sce3[,(sce3$cluster_full %in% c("2","3","4","5","6","8") & sce3$sub_atlas=="SCLC")]
sc2 <- sce3[,(sce3$cluster_full=="0" & sce3$sub_atlas=="SCLC")]
sce3_group <- data.frame(cellid=colnames(sce3),cluster=sce3$cluster_full)
sce3_group[colnames(sc1),"group"] <- "SCLC-1"
sce3_group[colnames(sc2),"group"] <- "SCLC-2"
sce3_group[-which(sce3_group$cellid %in% c(colnames(sc1),colnames(sc2))),"group"] <- "NSCLC"
sce3$group <- sce3_group$group

sclc <- sce3[,which(sce3$group %in% c("SCLC-1","SCLC-2"))]

# CytoTRACE for sclc
library(CytoTRACE)
library(ggprism)
library(ggpubr)

results <- CytoTRACE(as.matrix(sclc@assays$RNA@counts), ncores = 8)
dat <- data.frame(group=sclc$group,cytotrace=results$CytoTRACE)
ggplot(dat,aes(x=group,y=cytotrace))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=group),
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  ggtitle("boxplot")+
  scale_fill_prism(palette = "candy_bright")

t.test(dat[dat$group=="SCLC-1",2],dat[dat$group=="SCLC-2",2],var.equal=T)


# gene expression level
tmp <- sclc
Idents(tmp) <- tmp$group
sclc_marker <- FindAllMarkers(tmp,only.pos = T,min.pct = 0.1,logfc.threshold = 0.01,test.use = "wilcox")
scdeg <- sclc_marker %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

save(sce2,subcds,sce3,sclc,sclc_marker,file=file.path(dir,"monocle/point_sclcfull.RData"))


scmarker <- c("ASCL1","NEUROD1","POU2F3","YAP1")
VlnPlot(sclc,scmarker,group.by = 'group',flip = T,stack = T,split.by = "group")
ggsave("point_vlnplot1.pdf")

neuromarker <- c("DLL3","INSM1","HES6","CHGA")
VlnPlot(sclc,neuromarker,group.by = 'group',flip = T,stack = T,split.by = "group")
ggsave("point_vlnplot2.pdf")

gene_list <- c("STMN2","GNG7","CD44","BCAT1","CCND2","PAGE2","PAGE5","PAGE2B")
VlnPlot(sclc,gene_list,group.by = 'group',flip = T,stack = T,split.by = "group")
ggsave("point_vlnplot3.pdf")

# GSEA
run_gsea <- function(deg,num=NULL){
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSEABase)
  dir <- "case_ca"
  
  dir <- getwd()
  geneList <- deg$avg_log2FC 
  names(geneList) <- toupper(rownames(deg))
  geneList <- sort(geneList,decreasing = T)
  geneset <- read.gmt(file.path(dir,"c6.all.v2023.2.Hs.symbols.gmt"))  
  egmt <- GSEA(geneList, TERM2GENE=geneset, 
               minGSSize = 1,
               pvalueCutoff = 0.99,
               verbose=FALSE)
  gsea_results <- egmt@result
  gsea_results$norm_nes <- -log10(gsea_results$p.adjust)*(gsea_results$NES)
  gsea_results <- gsea_results[gsea_results$p.adjust<0.05,]
  
  library(viridis)
  ggplot(gsea_results,aes(-log10(P.adj)*(NES),Description))+
    geom_bar(aes(y=reorder(Description,norm_nes),x=norm_nes,fill=NES),stat='identity',width=0.5)+
    theme_bw()+
    scale_fill_viridis(option="rocket",direction = -1)+
    theme(panel.grid = element_blank())
  ggsave(paste0("point_gsea",num,".pdf"))
}

lc1.genes<- sclc_marker[sclc_marker$cluster=="SCLC-1",]
lc2.genes<- sclc_marker[sclc_marker$cluster=="SCLC-2",]
run_gsea(lc1.genes,1)
run_gsea(lc2.genes,2)




## pseudotime gene analysis for ad and sc in state 1,2,5
library(Seurat)
library(tidyverse)
library(monocle)
dir <- "case_ca"
load(file.path(dir,"monocle/multica_cds.RData"))

run_deg <- function(state){
  data <- cds[,cds$State %in% state]
  expre_genes <- row.names(subset(fData(data),num_cells_expressed >= 10))
  branch_deg <- differentialGeneTest(data[intersect(expre_genes,rownames(data)),],fullModelFormulaStr = "~Pseudotime",cores = 4) %>% subset(qval<0.01)
  branch_deg <- branch_deg[order(branch_deg$qval,decreasing = F),]
  return(branch_deg)
}

subcds <- cds[,cds$State %in% c(1,2,5)]
branch_deg <- run_deg(c(1,2,5))


sce <- readRDS(file.path(dir,"monocle/multica_umap.rds"))
sce2 <- sce[,colnames(subcds)]
sce2 <- FindNeighbors(sce2,reduction="scanvi",dims=1:10)
sce2 <- RunUMAP(sce2, reduction="scanvi", dims=1:10) %>% FindClusters(resolution = 0.01)
marker <- FindAllMarkers(object = sce2,only.pos = T,min.pct = 0.1,logfc.threshold = 0.25,test.use = "wilcox")
top20 <- marker %>% group_by(cluster) %>% top_n(n = 20,wt = avg_log2FC)

track_gene <- branch_deg[intersect(rownames(subcds),top20$gene),]
track_gene2 <- unique(track_gene$gene_short_name) %>% na.omit()



library(viridis)
library(ggsci)
library(cowplot)

plot_pseudotime_heatmap(subcds[track_gene2,], 
                        num_cluster = 3, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols = magma(100))

top15 <- marker %>% group_by(cluster) %>% top_n(n = 15,wt = avg_log2FC)
time_gene <- branch_deg[intersect(rownames(subcds),top15$gene),]
time_gene <- time_gene[order(time_gene$qval,decreasing = F),]
time_gene2 <- rownames(time_gene)

gene_list <- c("CXCL6","ALDH3B1","SPRR1B","KRT7","SLC25A12","SFTA2")

plot_genes_in_pseudotime(subcds[gene_list,], color_by="State", min_expr=0.5, ncol = 2)+
  scale_color_nejm()
ggplot2::ggsave(filename = "pseudo_gene.pdf",width = 14,height = 8)

plot_genes_in_pseudotime(subcds[time_gene2[1:10],], color_by="Pseudotime", min_expr=0.5, ncol = 2)+
  scale_color_gradientn(colours = c('#336699','#66CC66','#FFCC33'))
ggsave(p,filename = file.path(dir,paste0("monocle/pseudo_heatmap.pdf")),width = 5,height = 6)


## analysis for ad,inter ad-sc, and sc
sce3 <- sce2[,sce2$sub_atlas %in% c("LUAD","LUSC")]
plot_grid(ncol=2,
          DimPlot(sce3,reduction = "umap",cols = pal_lancet("lanonc")(9),
                  label = F,group.by = "seurat_clusters")+ggtitle("Cluster"),
          DimPlot(sce3,reduction = "umap",cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                  label = F,group.by = "sub_atlas")+ggtitle("Cancer type"))




inter_ad_sc <- sce3[,(sce3$seurat_clusters=="0" & sce3$sub_atlas=="LUSC") |
                      (sce3$seurat_clusters=="1" & sce3$sub_atlas=="LUAD")]
sce3_group <- data.frame(cellid=colnames(sce3),cluster=sce3$seurat_clusters)
sce3_group[which(sce3_group$cluster=="0"),"group"] <- "LUAD"
sce3_group[which(sce3$seurat_clusters %in% c("1","2","3")),"group"] <- "LUSC"
sce3_group[colnames(inter_ad_sc),"group"] <- "Transition_LUAD_LUSC"
sce3$group <- sce3_group$group

tmp <- sce3
Idents(tmp) <- tmp$group
sce3_marker <- FindAllMarkers(tmp,only.pos = T,min.pct = 0.1,logfc.threshold = 0.25,test.use = "wilcox")
cadeg <- sce3_marker %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
final_gene <- intersect(cadeg$gene,branch_deg$gene_short_name)
final_gene2 <- c("NTS","CSTA","KRT7","DSG3","FBN2","KRT8P4","GATA1","SPRR1B",
                 "PLCXD2","KRT18P43","MYH3","CXXC1","CXCL6","HBB","HBA2",     
                 "SFTPA2","SFTPA1","SFTPB","SPINK1","SFTA2","SLC25A12","CTSE","AQP1","AGR3","TFF1")

save(cds,subcds,branch_deg,sce2,sce3,sce3_marker,marker,file = file.path(dir,paste0("monocle/pseudo_gene.RData")))

DotPlot(sce3,features = final_gene2,group.by = "group")+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))+
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggplot2::ggsave(filename = "pseudo_dot.pdf",width = 10,height = 8)




results <- CytoTRACE(as.matrix(sce3@assays$RNA@counts), ncores = 8)
emb <- sce3@reductions[["umap"]]@cell.embeddings

group <- as.character(sce3$group)
names(group) <- rownames(sce3@meta.data)
plotCytoTRACE(results, phenotype = group,emb = emb,outputDir = paste0(dir,"/monocle/pseudo_"))

# subset cytotrace score for boxplot
library(ggprism)
library(ggpubr)
dat <- data.frame(group=sce3$group,cytotrace=results$CytoTRACE)

ggplot(dat,aes(x=group,y=cytotrace))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=group),
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  ggtitle("boxplot")+
  stat_compare_means(comparisons = list(c("LUAD", "Transition_LUAD_LUSC"),
                                        c("LUAD","LUSC"),
                                        c("Transition_LUAD_LUSC","LUSC")),
                     size = 6,
                     method = "t.test",label = "p.value")+
  scale_fill_prism(palette = "candy_bright")

t.test(dat[dat$group=="LUAD",2],dat[dat$group=="Transition_LUAD_LUSC",2],var.equal=T)
t.test(dat[dat$group=="Transition_LUAD_LUSC",2],dat[dat$group=="LUSC",2],var.equal=T)
t.test(dat[dat$group=="LUSC",2],dat[dat$group=="LUAD",2],var.equal=T)
               



run_cyto_multica <- function(x){
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(monocle)
  dir <- "case_ca"
  
  load(file.path(dir,paste0("monocle/",x,"_cds.RData")))
  
  print("=============== Start cytotrace ===============")
  library(CytoTRACE)
  mat <- as.matrix(sce_sub@assays$RNA@counts)
  results <- CytoTRACE(mat, ncores = 8)
  
  monocle_meta <- data.frame(t(cds@reducedDimS))
  colnames(monocle_meta) <- c("C1", "C2")
  emb_monocle <- monocle_meta[,1:2]
  
  print("=============== Cytotrace plot ===============")
  sub_atlas <- as.character(sce_sub$sub_atlas)
  names(sub_atlas) <- rownames(sce_sub@meta.data)
  plotCytoTRACE(results, phenotype = sub_atlas,emb = emb_monocle,outputDir = paste0(dir,"/monocle/",x,"_subca_"))
  
  stage <- as.character(sce_sub$stage)
  names(stage) <- rownames(sce_sub@meta.data)
  plotCytoTRACE(results, phenotype = stage, emb = emb_monocle, outputDir = paste0(dir,"/monocle/",x,"_stage_"))
}
run_cyto_multica("multica")
