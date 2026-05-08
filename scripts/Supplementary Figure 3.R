#########################################################################################
#### Fig.S3A
#########################################################################################
setwd("rna-analysis/analysis/")
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(corrplot)
library(openxlsx)
library(ggrepel)
library(enrichplot)
library(cowplot)

sample.re<-read.xlsx("../../ref_data/sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])

### PCA 相关性 确认样本----
data<-read.table("../gene-TPM-matrix.mm10.txt",check.names = F)
colname.new<-lapply(colnames(data),function(x){strsplit(strsplit(x,"rna-analysis/mapping_mm10/")[[1]][2],".genes.results")[[1]][1]})%>%as.character()
colnames(data)<-colname.new
colnames(data)<-sample.re[colnames(data),]
tpm<-data
#去除未表达gene
data<-data[rowSums(data)!=0,]

## 样本一致性----
corr<-cor(data,method = "spearman")
pdf("Nomal-Leukemia.corrplot.pdf",width = 6.4,height = 6)
pheatmap::pheatmap(corr,clustering_method = "ward.D2",border_color = NA,
                   color =colorRampPalette(c("skyblue", "salmon","red"))(100) ,cutree_rows = 5,cutree_cols = 5)
dev.off()
# 存在异质性，尤其是样本间
## PCA----
## pca plot 
# data<-log2(tpm+1)
samples<-colnames(tpm)
sample_use<-setdiff(samples,c("N_CD8_0","N_DP_0","L0_CD8_1", "L0_DP_1",  "L0_CD8_2", "L0_DP_2"))


data<-tpm[,sample_use] #排除L0样本
# data<-count
loc<-rowSums(data)!=0
data.filter<-data[loc,]
{
  pca<-prcomp(t(data.filter),scale. = T)
  head(pca$x)
  summary(pca)
  pca_result<-as.data.frame(pca$x[,1:3])
  sample<-colnames(data.filter)
  shape<-lapply(sample, function(x){strsplit(x,split = "_")[[1]][1]})%>%as.character()
  group<-lapply(sample, function(x){strsplit(x,split = "_")[[1]][2]})%>%as.character()
  pca_result<-cbind(pca_result,sample,shape,group)
  
  # 2d pca
  colnames(pca_result)<-c("PC1" ,   "PC2",    "PC3"  ,  "name", "sample" , "class")
  pca_result[pca_result$sample!="N",]$sample<-"L"

  # ONLY CD8
  library(ggforce)
  pca_result_cd8<-pca_result %>% filter(class=="CD8")
  pca_result_DP_N<-pca_result %>% filter(class=="DP"&sample=="N")
  plot<-rbind(pca_result_cd8,pca_result_DP_N)
  p<-ggplot(plot,aes(x=PC1,y=PC2,shape=class,color=sample))+
    geom_point(size=5)+
    geom_mark_ellipse(aes(color=sample),linetype="dashed")+
    geom_text_repel(aes(label=name,x=plot[,1],y=plot[,2]),data=plot,size=2)+
    theme_bw()+
    theme(
      plot.margin = unit(c(1,1,1,1),'lines'),
      panel.grid = element_blank())+
    ylim(c(-150,100))+
    xlim(c(-100,160))+
  labs(x="PCA 1",y="PCA 2")
  
  ggsave(plot=p,"out/Fig.S3A L-N.CD8-N.PCA.pdf",width = 5.5,height = 3.8)
}


#########################################################################################
#### Fig.S3B
#########################################################################################
library(Seurat)
library(patchwork)

seu<-readRDS("clone.inte.ccr.rds")
seu$sample<-"sample1"
seu@meta.data[seu@meta.data$orig.ident %in% c("Thy_2","DN_2"),]$sample<-"sample2"

p1<-DimPlot(seu,group.by = "sample",pt.size = 1.2)
p2<-DimPlot(seu,group.by = "seurat_clusters",label = T,repel = T)+ggtitle("cluster")
p3<-DimPlot(seu,label = T,repel = T)+ggtitle("celltype")

p1|p2|p3
ggsave(filename = "out/Fig.S3B.pdf",width = 15,height = 4)


#########################################################################################
#### Fig.S3C
#########################################################################################
DotPlot(seu,feature=c("Ptprc","Cd3d","Cd8a","Cd4","Il2ra","Kit","Cd44","Notch1","Hes1","Wnt5b","Mki67","Bcl11a","Bcl11b"
+ ),cols  = c("lightblue","red"))+coord_flip()+theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0))
ggsave(filename = "out/Fig.S3C marker.pdf",width = 4.5,height = 5)


#########################################################################################
#### Fig.S3D
#########################################################################################

library(Seurat)
library(monocle3)
library(ggplot2)
library(scales)
library(ggthemes)
library(dplyr)

theme_blank<-theme(panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = margin(1,1,1,1,unit = "cm"))


data<-subset(seu,group=="S1")
cell_metadata <- data@meta.data
mtx <- GetAssayData(data, assay = "RNA", slot = 'counts')%>%as.matrix()
gene_annotation <- data.frame(gene_short_name = rownames(mtx))
int.embed <- Embeddings(data, reduction = "umap")
rownames(gene_annotation) <- rownames(mtx)

save(cell_metadata,mtx,gene_annotation,int.embed,file = "toDELL.monocle3-need.clone.S1.rdata")


data<-subset(seu,group=="S2")
cell_metadata <- data@meta.data
mtx <- GetAssayData(data, assay = "RNA", slot = 'counts')%>%as.matrix()
gene_annotation <- data.frame(gene_short_name = rownames(mtx))
int.embed <- Embeddings(data, reduction = "umap")
rownames(gene_annotation) <- rownames(mtx)

save(cell_metadata,mtx,gene_annotation,int.embed,file = "toDELL.monocle3-need.clone.S2.rdata")


load("toDELL.monocle3-need.clone.S1.rdata")
# load("cds.rdata")

cds <- new_cell_data_set(mtx,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 15)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)

plot_cells(cds, color_cells_by="partition")

# 调整umap-采用原始的umap
cds.embed <- cds@int_colData$reducedDims$UMAP

int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="useIdent") + ggtitle('int.umap')


cds <- cluster_cells(cds)
cds <- learn_graph(cds,learn_graph_control=list(minimal_branch_len=5, geodesic_distance_ratio=0.5, euclidean_distance_ratio=1))
# cds <- learn_graph(cds)

# save(cds,file = "cds.rdata")
cds <- order_cells(cds)
S1.cds<-cds
cds<-S1.cds
library(ggrastr)

p1<-plot_cells(cds, reduction_method="UMAP", color_cells_by="useIdent",rasterize = T,
               cell_size=1.5,trajectory_graph_color="black",label_cell_groups=F,
               graph_label_size=3,label_leaves=F,label_branch_points=F) +theme_bw()+
  scale_color_manual(values = hue_pal()(9))+theme_blank

p2<-plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",label_roots=T,rasterize = T,
               cell_size=1.5,trajectory_graph_color="black",
               graph_label_size=3) +theme_bw()+theme_blank

p3<-plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",show_trajectory_graph=FALSE,rasterize = T,
               cell_size=1.5) +theme_bw()+theme_blank

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "S1.clone.monocle3.raster.pdf",width = 18,height = 5)

#### S2 ----

library(Seurat)
library(monocle3)
library(ggplot2)
library(scales)
library(ggthemes)
library(dplyr)

theme_blank<-theme(panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = margin(1,1,1,1,unit = "cm"))

load("toDELL.monocle3-need.clone.S2.rdata")
# load("cds.rdata")

cds <- new_cell_data_set(mtx,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 15)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)

plot_cells(cds, color_cells_by="partition")

# 调整umap-采用原始的umap
cds.embed <- cds@int_colData$reducedDims$UMAP

int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="useIdent") + ggtitle('int.umap')


cds <- cluster_cells(cds)
cds <- learn_graph(cds,learn_graph_control=list(minimal_branch_len=5, geodesic_distance_ratio=0.5, euclidean_distance_ratio=1))
# cds <- learn_graph(cds)

# save(cds,file = "cds.rdata")
cds <- order_cells(cds)

S2.cds<-cds



### 作图 ----
cds<-S1.cds
# library(ggrastr)

p1<-plot_cells(cds, reduction_method="UMAP", color_cells_by="useIdent",
               # rasterize = T,
               cell_size=1.5,trajectory_graph_color="black",label_cell_groups=F,
               graph_label_size=3,label_leaves=F,label_branch_points=F) +theme_bw()+
  scale_color_manual(values = hue_pal()(9))+theme_blank

p2<-plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",label_roots=T,
               # rasterize = T,
               cell_size=1.5,trajectory_graph_color="black",
               graph_label_size=3) +theme_bw()+theme_blank

p3<-plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",show_trajectory_graph=FALSE,
               # rasterize = T,
               cell_size=1.5) +theme_bw()+theme_blank

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "out/Fig.S3D S1.clone.monocle3.pdf",width = 18,height = 5)

cds<-S2.cds

p1<-plot_cells(cds, reduction_method="UMAP", color_cells_by="useIdent",
               # rasterize = T,
               cell_size=1.5,trajectory_graph_color="black",label_cell_groups=F,
               graph_label_size=3,label_leaves=F,label_branch_points=F) +theme_bw()+
  scale_color_manual(values = hue_pal()(9))+theme_blank

p2<-plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",label_roots=T,
               # rasterize = T,
               cell_size=1.5,trajectory_graph_color="black",
               graph_label_size=3) +theme_bw()+theme_blank

p3<-plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",show_trajectory_graph=FALSE,
               # rasterize = T,
               cell_size=1.5) +theme_bw()+theme_blank

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "out/Fig.S3D S2.clone.monocle3.pdf",width = 18,height = 5)







#########################################################################################
#### Fig.S3E
#########################################################################################

setwd("~/Project/YuYong_TALL/out")
library(openxlsx)
library(homologene)
library(dplyr)
library(ggsignif)

gene<-read.xlsx("41467_2025_61222_MOESM5_ESM.xlsx") %>% filter(cluster=="Stem_like") %>% pull(gene) %>% 
  human2mouse(.) %>%pull(mouseGene) %>% unique()

seu<-AddModuleScore(seu,features = list(stemness=gene),name = "stemness",assay = "RNA")
seu$stemness<-seu$stemness1
seu$stemness1<-NULL
VlnPlot(seu,features = "stemness",pt.size = 0)

df<-seu@meta.data %>% dplyr::select(ident,stemness)
ggplot(df,aes(x=ident,y=stemness))+
  stat_boxplot(geom="errorbar",width=0.2)+
  # geom_jitter(aes(x=ident,y=stemness),width = 0.1,size=0.8,alpha=0.1)+
  geom_boxplot(outlier.alpha = 0,width=0.6,aes(fill=ident))+
  geom_signif(comparisons = list(c("DN2 (11)", "DN4 (12)"),
                                 c("DN2 (11)", "DP (3)"),
                                 c("DN2 (11)", "CD8+ (other)")),
              y_position=c(0.28, 0.26, 0.24), 
              tip_length = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
              map_signif_level=TRUE)+

  # geom_boxplot(width=0.6,outlier.size = 0,outlier.alpha = 0)+
  theme_classic()+
  ylab("Stemness score")+
  xlab("")+
  theme(axis.text.x = element_text(size = 12,angle = -45,hjust = 0))+
  NoLegend()
  # scale_fill_manual(values = RColorBrewer::brewer.pal(4,"Set2"))
# ggsave(filename = "Stemness.scRNA-seq.pdf",width = 3.5,height = 4)
ggsave(filename = "out/Fig.S3E Stemness.pdf",width = 3.5,height = 4)