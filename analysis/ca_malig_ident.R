library(reticulate)
library(Seurat)
library(tidyverse)
library(cowplot)
dir <- "case_ca"

use_condaenv("envs/scvi")
scanpy <- import("scanpy")
pd <- import("pandas")
np <- import("numpy")


in_ad <- scanpy$read_h5ad(file.path(dir,"epi_raw.h5ad"))
count <- np$transpose(in_ad$X)
colnames(count) <- rownames(in_ad$obs)
rownames(count) <- rownames(in_ad$var)
meta <- read.csv(file.path(dir,"epi_raw.csv"),row.names = 1)
seu <- CreateSeuratObject(counts = count, meta.data = meta)

anvi <- in_ad$obsm[["X_scANVI"]]
colnames(anvi) <- paste0("scanvi_", 1:10)
rownames(anvi) <- colnames(seu)
seu[["scanvi"]] <- CreateDimReducObject(embeddings = anvi, key = "scanvi_", assay = "RNA")
saveRDS(seu,file.path(dir,"epi_raw.rds"))


## infercnv
library(infercnv)
library(AnnoProbe)
library(reticulate)
dir <- "case_ca"

use_condaenv("envs/scvi")
scanpy <- import("scanpy")
np <- import("numpy")

sce <- readRDS(file.path(dir,'epi_new.rds'))
ref_epi <- scanpy$read_h5ad(file.path(dir,"ref_epi10k.h5ad"))
epi_dat <- np$transpose(ref_epi$X)
colnames(epi_dat) <- rownames(ref_epi$obs)
rownames(epi_dat) <- rownames(ref_epi$var)


dat <- as.data.frame(sce@assays$RNA@counts)
colnames(dat) <- colnames(sce)
rownames(dat) <- rownames(sce)

dat2 <- dat[intersect(rownames(dat),rownames(epi_dat)),]
epi_dat2 <- epi_dat[intersect(rownames(dat),rownames(epi_dat)),]

all_dat <- cbind(dat2,epi_dat2)
groupinfo <- data.frame(v1=colnames(all_dat),
                        v2=c(sce$annotation,rep('ref-epithelial',10000)))

geneInfor <- annoGene(rownames(all_dat),"SYMBOL",'human')
geneInfor <- geneInfor[with(geneInfor,order(chr,start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]

data <- all_dat[rownames(all_dat) %in% geneInfor[,1],]
data <- data[match(geneInfor[,1], rownames(data)),] 

print("=============== Start save files ===============")
expfile <- file.path(paste0(dir,'/infercnv'),'expfile.txt')
write.table(data,file = expfile,sep = '\t',quote = F)
groupfile <- file.path(paste0(dir,'/infercnv'),'groupfile.txt')
write.table(groupinfo,file = groupfile,sep = '\t',quote = F,col.names = F,row.names = F)
genefile <- file.path(paste0(dir,'/infercnv'),'genefile.txt')
write.table(geneInfor,file = genefile,sep = '\t',quote = F,col.names = F,row.names = F)


options("Seurat.object.assay.version" = "v3")
options(future.globals.maxSize= 5*1024^3)
options(scipen = 100)

expfile <- file.path(paste0(dir,'/infercnv'),'expfile.txt')
groupfile <- file.path(paste0(dir,'/infercnv'),'groupfile.txt')
genefile <- file.path(paste0(dir,'/infercnv'),'genefile.txt')

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expfile,
                                    annotations_file=groupfile,
                                    delim="\t",
                                    gene_order_file= genefile,
                                    ref_group_names='ref-epithelial') 
print("=============== Done create cnvdata ===============")

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir=paste0(dir,"/infercnv/cnv_all2"),
                              cluster_by_groups=TRUE,
                              no_prelim_plot=TRUE,
                              write_expr_matrix=TRUE,
                              write_phylo = TRUE,
                              denoise = TRUE,
                              HMM = FALSE,
                              output_format = "pdf")



## cnv score
cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}


dat <- readRDS("cnv_all/run.final.infercnv_obj")
cnv_score <- cnvScore(dat@expr.data)

colnames(cnv_score) <- "score"
cnv_score <- rownames_to_column(cnv_score, var='cellid')

group <- read.table("case_ca/infercnv/groupfile.txt")
colnames(group) <- c("cellid","cluster")

cnv_score <- left_join(cnv_score,group,by="cellid")


ggplot(cnv_score,aes(x=cluster, y=score,fill=cluster))+geom_violin(aes(fill=cluster),color="NA")+
  scale_fill_igv()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = .4),
        axis.ticks.x=element_line(color="black",linewidth = .4,lineend = 1),
        axis.ticks.y=element_line(color="black",linewidth = .4,lineend = 1),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 9))+
  ylab("cnv score")+xlab("celltype")
ggsave("cnv_score.pdf",width = 15,height = 6,units = "cm")




