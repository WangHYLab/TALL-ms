#########################################################################################
#### Fig.S3A
#########################################################################################
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(homologene)

normal<-readRDS("data/normal.scenic.RDS")
normal
DimPlot(normal,label = T,group.by = "celltype_normal")

normal<-RunTSNE(normal,dims = 1:20)
normal<-RunUMAP(normal,dims = 1:15,min.dist = 2)
p1<-DimPlot(normal,label = T,reduction = "umap",group.by = "seurat_clusters",pt.size = 0.8)
p2<-DimPlot(normal,label = T,reduction = "umap",group.by = "celltype_normal",pt.size = 0.8)

ggsave(plot = p2, filename = "out/Fig.S3A normal.celltype.pdf", width = 5.5, height = 5)

p1+p2
DotPlot(normal,features = c(
  "Hhex","Cd44",
  "Kit","Il7r","Il2ra","Cd4","Cd8a","Cd8b1",
  "Bcl11a","Bcl11b","Notch1",
  "Rag1","Rag2","Gata3","Dntt","Cd247","Ptprc",
  # "Cd52",
  "Cd3d",
  # "Cd3e",
  "Tcf7",
  "Ms4a1","Cd19","Epcam","Krt18","S100a8","Cd68","Itgam"
),
cluster.idents = T,group.by = "seurat_clusters")+
  xlab("")+ylab("")+theme_bw()+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0,size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black")
  )
DimPlot(normal,label = T,reduction = "umap",group.by = "orig.ident",pt.size = 0.8)
saveRDS(normal,"normal.rds")

#########################################################################################
#### Fig.S3B
#########################################################################################
normal<-readRDS("normal.rds")

DotPlot(normal,features = c(
                            "Kit","Il7r","Il2ra","Cd4","Cd8a","Cd8b1",
                            "Bcl11b","Notch1","Rag1","Rag2","Gata3","Dntt","Cd247","Ptprc",
                            "Cd3d","Tcf7"
                            ),cols = c("lightblue","salmon"),
        cluster.idents = F)+
  xlab("")+ylab("")+theme_bw()+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.5,hjust = 0,size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black")
  )

# cluster17 isolation but expression of maker support it is DN3 cell
ggsave(filename = "out/Fig.S3B normal.marker.dotplot.pdf", width = 6.5, height = 3.1)





#########################################################################################
#### Fig.S3C-E
#########################################################################################
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(homologene)


setwd("~/Project/YuYong_TALL/out/normal_tall")

normal <- readRDS("~/Project/YuYong_TALL/out/SCENIC_withCTL/normal.scenic.RDS")
normal <- RunUMAP(normal, dims = 1:15, min.dist = 2)
DimPlot(normal, label = T, reduction = "umap", group.by = "celltype_normal", pt.size = 0.8)

{
    data.list <- readRDS("~/Project/YuYong_TALL/data.list.filtered&removeDB.rds")
    seu <- merge(data.list$DN1, c(data.list$DN2, data.list$Thy1, data.list$Thy2))

    DefaultAssay(seu) <- "RNA"
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu)
    seu <- RunUMAP(seu, dims = 1:20)

    # Calculating Cell cycle score
    seu <- JoinLayers(seu)
    seu <- CellCycleScoring(seu, s.features = s.gene, g2m.features = g2m.gene)

    # Regress out cell cycle score
    # seu$CC.Difference <- seu$S.Score - seu$G2M.Score
    # seu <- ScaleData(seu, vars.to.regress = "CC.Difference", features=VariableFeatures(seu))
    seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(seu))
    seu <- RunPCA(seu, reduction.name = "pca.ccregress")

    # Perform RPCA integration
    seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

    # options(future.globals.maxSize = 8000 * 1024^2)
    seu <- IntegrateLayers(
        object = seu,
        method = CCAIntegration,
        dims = 1:20,
        orig.reduction = "pca.ccregress",
        new.reduction = "integrated.cca",
        verbose = FALSE
    )

    # re-join layers after integration
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
    seu <- RunUMAP(seu, dims = 1:20, reduction = "integrated.cca")
    seu <- FindNeighbors(seu, dims = 1:20, reduction = "integrated.cca")
    seu <- FindClusters(seu, resolution = 0.8)
    DimPlot(seu)
    saveRDS(seu, "data/SCIN.ccregress.integrate.rds")
}
tall <- readRDS(file = "data/SCIN.ccregress.integrate.rds")

DimPlot(tall, label = T)
DotPlot(tall, features = c("Ms4a1", "Cd19", "Cd68", "Cd86", "Ptprc", "Cd3d", "Cd52", "S100a8", "Csf1r", "Klrc1"))
## remove 非T细胞亚群 15
tall <- subset(tall, seurat_clusters != "15")
DimPlot(tall, label = T)
#### 合并
tall$celltype <- "Tall"
normal$celltype <- as.character(normal$celltype_normal)


# 正确合并数据集（Seurat v5推荐方式）
seu <- merge(x = tall, y = normal, add.cell.ids = c("TALL", "Normal"), project = "TALL_vs_Normal")

# 添加分组信息（合并后立即添加，避免后续混淆）
seu$group <- ifelse(seu$orig.ident %in% c("N_DN", "N_THY"), "Normal", "TALL")
table(seu$group) # 验证分组是否正确

# 设置默认Assay
DefaultAssay(seu) <- "RNA"


# 2. 统一标准化与可变基因选择
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)


# 3. 细胞周期评分与校正
s.genes <- human2mouse(cc.genes.updated.2019$s.genes) %>%
    pull(mouseGene) %>%
    unique()
g2m.genes <- human2mouse(cc.genes.updated.2019$g2m.genes) %>%
    pull(mouseGene) %>%
    unique()
# 计算细胞周期评分（必须在分层前计算，否则会报错）
seu <- JoinLayers(seu)
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)

# 查看细胞周期分布
table(seu$Phase)
DimPlot(seu, reduction = "pca", group.by = "Phase") # 先跑一次PCA看原始周期效应

seu.bk <- seu
# 按种类TALL和normal去批次
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$group)

# 分层后重新标准化和选择可变基因
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)

seu.bk2 <- seu
seu <- ScaleData(seu,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = VariableFeatures(seu), verbose = FALSE
)


seu.bk3 <- seu
# PCA
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)


# plot
ElbowPlot(seu, ndims = 50)
# 1:20
dims.use <- 1:20


# 5. Harmony
seu.inte <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    dims = dims.use
)


# plot
DimPlot(seu.inte, reduction = "harmony", group.by = "group") + ggtitle("Harmony整合后")
DimPlot(seu.inte, reduction = "harmony", group.by = "Phase") + ggtitle("细胞周期校正后")


# 6. reduction

seu.inte <- RunUMAP(seu.inte, reduction = "harmony", dims = dims.use, n.neighbors = 30)
DimPlot(seu.inte, group.by = "group") + ggtitle("Harmony整合后")

seu.inte <- FindNeighbors(seu.inte, reduction = "harmony", dims = dims.use)
seu.inte <- FindClusters(seu.inte, resolution = 0.5) # 生成多个分辨率供选择


p1 <- DimPlot(seu.inte, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + ggtitle(label = "")
DimPlot(seu.inte, reduction = "umap", group.by = "group") + ggtitle(label = "")

seu.inte$sample <- "Normal"
seu.inte@meta.data[seu.inte@meta.data$orig.ident %in% c("DN_1", "Thy_1"), ]$sample <- "TALL-1"
seu.inte@meta.data[seu.inte@meta.data$orig.ident %in% c("DN_2", "Thy_2"), ]$sample <- "TALL-2"

cols <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F")

p2 <- DimPlot(seu.inte,
    reduction = "umap",
    group.by = "sample", label = T,
    cols = c("#8DD3C7", "salmon", "#BC80BD")
) + ggtitle(label = "")

DimPlot(seu.inte, reduction = "umap", group.by = "celltype", label = T)

ggsave(plot = p1 + p2, filename = "normal&tall.cluster.sample.umap.pdf", width = 10, height = 4.5)

## assign

DotPlot(seu.inte,
    features = c(
        "Hhex", "Cd44",
        "Kit", "Il7r", "Il2ra", "Cd4", "Cd8a", "Cd8b1",
        "Bcl11a", "Bcl11b", "Notch1",
        "Rag1", "Rag2", "Gata3", "Dntt", "Cd247", "Ptprc",
        # "Cd52",
        "Cd3d",
        # "Cd3e",
        "Tcf7",
        "Ms4a1", "Cd19"
    ),
    # cols = c("lightblue","red"),
    cluster.idents = T
) +
    coord_flip() + xlab("") + ylab("") + theme_bw() +
    theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
    )

seu.inte@meta.data$cell_type <- NA


## no T cell
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(17), ]$cell_type <- "B cell"

## DN
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(2, 16), ]$cell_type <- "DN3"
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(14), ]$cell_type <- "DN4"
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(10), ]$cell_type <- "DN1/DN2"

## DP
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(1, 7, 8, 11, 12, 13), ]$cell_type <- "DP"

## SP
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(0, 3, 4, 5, 6, 9, 15), ]$cell_type <- "CD8"

table(seu.inte$cell_type, useNA = "ifany")

table(seu.inte$cell_type, seu.inte$seurat_clusters)

Idents(seu.inte) <- seu.inte$cell_type

# UMAP

seu.inte$cell_type <- factor(seu.inte$cell_type, levels = c("DN1/DN2", "DN3", "DN4", "DP", "CD8", "B cell"))
p3 <- DimPlot(seu.inte, label = TRUE, label.size = 6, repel = TRUE, group.by = "cell_type")
DotPlot(seu.inte,
    features = c(
        # "Hhex",
        "Cd44",
        "Kit", "Il7r", "Il2ra", "Cd4", "Cd8a", "Cd8b1",
        "Bcl11a", "Bcl11b", "Notch1",
        "Rag1", "Rag2", "Gata3", "Dntt", "Cd247", "Ptprc",
        # "Cd52",
        "Cd3d",
        # "Cd3e",
        "Tcf7",
        "Ms4a1", "Cd19"
    ), group.by = "cell_type",
    # cols = c("lightblue","red"),
    cluster.idents = F
) +
    coord_flip() + xlab("") + ylab("") + theme_bw() +
    theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
    )

# ggsave(plot = p1 + p2 + p3, filename = "normal&tall.cluster&sample&ct.umap.pdf", width = 15.5, height = 4.5)


## try to re-cluster
seu.inte2 <- FindClusters(seu.inte, resolution = 0.8) # 生成多个分辨率供选择


# remove B cell
out <- subset(seu.inte, cell_type != "B cell")
out$cell_type <- factor(out$cell_type, levels = c("DN1/DN2", "DN3", "DN4", "DP", "CD8"))
p1 <- DimPlot(out, reduction = "umap", label = TRUE, label.size = 6, raster = T, pt.size = 0.9, raster.dpi = c(500, 500), group.by = "seurat_clusters") + ggtitle(label = "")

p2 <- DimPlot(out,
    reduction = "umap",
    group.by = "sample", label = T, label.size = 6, raster = T, pt.size = 0.9, raster.dpi = c(500, 500),
    cols = c("#8DD3C7", "salmon", "#BC80BD")
) + ggtitle(label = "")
p3 <- DimPlot(out, label = TRUE, label.size = 6, repel = TRUE, raster = T, raster.dpi = c(500, 500), pt.size = 0.9, group.by = "cell_type") + ggtitle(label = "")
# ggsave(plot = p1 + p2 + p3, filename = "Fig.S3C-E normal&tall.cluster&sample&ct.umap.noBcell.raster.pdf", width = 15.5, height = 4.5)

p1 <- DimPlot(out, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.9, raster.dpi = c(500, 500), group.by = "seurat_clusters") + ggtitle(label = "")

p2 <- DimPlot(out,
    reduction = "umap",
    group.by = "sample", label = T, label.size = 6, pt.size = 0.9, raster.dpi = c(500, 500),
    cols = c("#8DD3C7", "salmon", "#BC80BD")
) + ggtitle(label = "")
p3 <- DimPlot(out, label = TRUE, label.size = 6, repel = TRUE, raster.dpi = c(500, 500), pt.size = 0.9, group.by = "cell_type") + ggtitle(label = "")

ggsave(plot = p1 + p2 + p3, filename = "out/Fig.S3C-E normal&tall.cluster&sample&ct.umap.noBcell.pdf", width = 15.5, height = 4.5)

saveRDS(out, file = "data/TALL_Normal.integrated.ccregressed.Tcell.rds")


#########################################################################################
#### Fig.S3F
#########################################################################################

out<-readRDS("data/TALL_Normal.integrated.ccregressed.Tcell.rds")

feat<-c('Notch1', 'Il7r', 'Kit' , 'Bcl11b','Cd4', 'Cd8a')
p<-FeaturePlot(out,features = feat,ncol = 3)
# ggsave(plot = p, filename = "out/Fig.S3F normal&tall.featureplot.umap.pdf", width = 15.5, height = 9)
p<-FeaturePlot(out,features = feat,ncol = 3,raster = T)
ggsave(plot = p, filename = "out/Fig.S3F normal&tall.featureplot.umap.raster.pdf", width = 15.5, height = 9)



#########################################################################################
#### Fig.S3G
#########################################################################################

library(Seurat)
library(patchwork)

seu <- readRDS("data/clone.inte.ccr.rds")
seu$sample <- "sample1"
seu@meta.data[seu@meta.data$orig.ident %in% c("Thy_2", "DN_2"), ]$sample <- "sample2"

p1 <- DimPlot(seu, group.by = "sample", pt.size = 1.2)
p2 <- DimPlot(seu, group.by = "seurat_clusters", label = T, repel = T) + ggtitle("cluster")
p3 <- DimPlot(seu, label = T, repel = T) + ggtitle("celltype")

p1 | p2 | p3
ggsave(filename = "out/Fig.S3G.pdf", width = 15, height = 4)

#########################################################################################
#### Fig.S3H
#########################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
seu <- readRDS("data/clone.inte.ccr.rds")
DotPlot(seu,
    feature = c(
        # "Ptprc","Cd3d","Cd8a","Cd4","Il2ra","Kit","Cd44","Notch1","Hes1","Wnt5b","Mki67","Bcl11a","Bcl11b"
        # "Hhex", "Cd44",
        "Kit", "Il7r", "Il2ra", "Cd4", "Cd8a", "Cd8b1",
        # "Bcl11a",
        "Bcl11b", "Notch1",
        "Rag1", "Rag2", "Gata3", "Dntt", "Cd247", "Ptprc",
        # "Cd52",
        "Cd3d",
        # "Cd3e",
        "Tcf7"
        # "Ms4a1", "Cd19"
    ),
    cols = c("lightblue", "red")
) + coord_flip() + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))

ggsave(filename = "out/Fig.S3H marker.pdf", width = 4.5, height = 5)




#########################################################################################
#### Fig.S3I
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

save(cell_metadata,mtx,gene_annotation,int.embed,file = "data/toDELL.monocle3-need.clone.S1.rdata")


data<-subset(seu,group=="S2")
cell_metadata <- data@meta.data
mtx <- GetAssayData(data, assay = "RNA", slot = 'counts')%>%as.matrix()
gene_annotation <- data.frame(gene_short_name = rownames(mtx))
int.embed <- Embeddings(data, reduction = "umap")
rownames(gene_annotation) <- rownames(mtx)

save(cell_metadata,mtx,gene_annotation,int.embed,file = "data/toDELL.monocle3-need.clone.S2.rdata")


load("data/toDELL.monocle3-need.clone.S1.rdata")
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

load("data/toDELL.monocle3-need.clone.S2.rdata")
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

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "out/Fig.S3I S1.clone.monocle3.pdf",width = 18,height = 5)

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

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "out/Fig.S3I S2.clone.monocle3.pdf",width = 18,height = 5)





#########################################################################################
#### Fig.S3J
#########################################################################################

normal<-readRDS("~/Project/YuYong_TALL/out/SCENIC_withCTL/normal.scenic.RDS")
normal<-RunUMAP(normal,dims = 1:15,min.dist = 2)
DimPlot(normal,label = T,reduction = "umap",group.by = "celltype_normal",pt.size = 0.8)



tall<-readRDS("~/Project/YuYong_TALL/clone.inte.ccr.rds")
DimPlot(tall)
tall$celltype<-as.character(tall@active.ident)
normal$celltype<-as.character(normal$celltype_normal)

data<-merge(tall,normal)

data@active.assay<-"RNA"

seu <- data
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst",nfeatures=3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
# seu <- RunUMAP(seu, dims=1:20)

# Calculating Cell cycle score
seu <- JoinLayers(seu)
seu$group<-"TALL"
seu@meta.data[seu@meta.data$orig.ident %in%c("N_DN", "N_THY"), ]$group<-"Normal"
# Regress out cell cycle score
s.gene   <- human2mouse(cc.genes.updated.2019$s.genes)   %>% pull(mouseGene) %>% unique()
g2m.gene <- human2mouse(cc.genes.updated.2019$g2m.genes) %>% pull(mouseGene) %>% unique()
seu <- CellCycleScoring(seu,s.features=s.gene,g2m.features = g2m.gene)

seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features=VariableFeatures(seu))
seu <- RunPCA(seu, reduction.name = "pca.ccregress")

# Perform harmony integration
seu[["RNA"]] <- split(seu[["RNA"]], f=seu$group)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst",nfeatures=3000)
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features=VariableFeatures(seu))
seu <- RunPCA(seu, reduction.name = "pca.ccregress")

# options(future.globals.maxSize = 8000 * 1024^2)
seu.inte <- IntegrateLayers(object = seu,
                       method = HarmonyIntegration,
                       dims=1:20,
                       orig.reduction = "pca.ccregress",
                       new.reduction = "harmony",
                       verbose = FALSE)

# re-join layers after integration
# seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
seu.inte <- RunUMAP(seu.inte, dims = 1:20,)
seu.inte <- FindNeighbors(seu.inte,dims = 1:20,reduction = "harmony")
seu.inte <- FindClusters(seu.inte,resolution = 0.8)

seu.inte$celltype<-factor(seu.inte$celltype,levels = c("ETP", "DN2","DN3","DN4","DP", "CD8+",
                                                       "DN2 (11)","DN4 (12)","DP (3)","CD8+ (other)"))
p1<-DimPlot(
  object = seu.inte,   
  reduction = "umap",  
  group.by = "celltype" ,
  label = T,pt.size = 0.5
  
)

p2<-DimPlot(
  object = seu.inte,   
  reduction = "umap",  
  group.by = "group" ,
  label = T,pt.size = 0.5
  
)
ggsave(plot = p1+p2,filename = "normal&tall.celltype.pdf",width = 13,height = 5)



saveRDS(seu.inte,"normal&CLONE.inte-harmony.ccr.rds")

seu.inte<-readRDS("normal&CLONE.inte-harmony.ccr.rds")
data<-seu.inte


DotPlot(seu.inte,features = c("Ptprc","Cd3d",#"Cd3g","Cd3e",
                         "Kit",#"Hlf","Hoxa9",
                         
                         "Il2ra","Bcl11a","Cd28",
                         # "Cd44",
                         "Cd24a",
                         "Cd4","Cd8a","Cd8b1",
                         
                         "Lef1","Bcl11b","Notch1"),cols = c("lightblue","red"),
        cluster.idents = F,group.by = "celltype")+
  coord_flip()+xlab("")+ylab("")+theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0,size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black")
  )
# ggsave(filename = "normal&tall.celltype.marker.pdf",width = 5,height = 4.5)




library(openxlsx)
library(homologene)
library(dplyr)
library(ggsignif)

gene<-read.xlsx("/home/zhengjie/Project/YuYong_TALL/out/41467_2025_61222_MOESM5_ESM.xlsx") %>% filter(cluster=="Stem_like") %>% pull(gene) %>% 
  human2mouse(.) %>%pull(mouseGene) %>% unique()

data[["RNA"]] <- JoinLayers(data[["RNA"]])

data<-AddModuleScore(data,features = list(stemness=gene),name = "stemness",assay = "RNA")
data$stemness<-data$stemness1
data$stemness1<-NULL
VlnPlot(data,features = "stemness",pt.size = 0)

df<-data@meta.data %>% dplyr::select(celltype,stemness) %>% mutate(ident=celltype)
ggplot(df,aes(x=celltype,y=stemness))+
  stat_boxplot(geom="errorbar",width=0.2)+
  # geom_jitter(aes(x=ident,y=stemness),width = 0.1,size=0.8,alpha=0.1)+
  geom_boxplot(outlier.alpha = 0,width=0.6,aes(fill=celltype))+
  # geom_signif(comparisons = list(c("DN2 (11)", "DN4 (12)"),
  #                                c("DN2 (11)", "DP (3)"),
  #                                c("DN2 (11)", "CD8+ (other)")),
  #             y_position=c(0.28, 0.26, 0.24), 
  #             tip_length = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
  #             map_signif_level=TRUE)+
  
  # geom_boxplot(width=0.6,outlier.size = 0,outlier.alpha = 0)+
  theme_classic()+
  ylab("Stemness score")+
  xlab("")+
  theme(axis.text.x = element_text(size = 12,angle = -45,hjust = 0),
        plot.margin =unit(c(1, 1, 1, 1), "cm"))+
  NoLegend()
# scale_fill_manual(values = RColorBrewer::brewer.pal(4,"Set2"))
ggsave(filename = "out/Fig.S3J TALL&Normal.Stemness.inte.ccr.new.pdf", width = 5, height = 4)


# significance test
df<-data@meta.data %>% dplyr::select(celltype,stemness) 
df


cell_types <- sort(unique(df$celltype))
n_type <- length(cell_types)


p_mat <- matrix(
  NA,
  nrow = n_type,
  ncol = n_type,
  dimnames = list(cell_types, cell_types)
)


for (i in seq_len(n_type)) {
  for (j in seq_len(n_type)) {
    if (i >= j) next 
    
    g1 <- df$stemness[df$celltype == cell_types[i]]
    g2 <- df$stemness[df$celltype == cell_types[j]]
    
    res <- wilcox.test(g1, g2)
    p_mat[i, j] <- res$p.value
    p_mat[j, i] <- res$p.value  
  }
}

diag(p_mat) <- 1


upper_pos <- upper.tri(p_mat, diag = FALSE)
p_raw_vector <- p_mat[upper_pos]
p_adj_vector <- p.adjust(p_raw_vector, method = "BH")


p_mat_adj <- p_mat
p_mat_adj[upper_pos] <- p_adj_vector
p_mat_adj[lower.tri(p_mat_adj)] <- t(p_mat_adj)[lower.tri(p_mat_adj)]


print(p_mat, digits = 3)

print(p_mat_adj, digits = 3)

openxlsx::write.xlsx(list(p=p_mat %>% data.frame()%>% tibble::rownames_to_column("celltype"),
                          padj=p_mat_adj  %>% data.frame()%>% tibble::rownames_to_column("celltype")), 
                     file = "out/Fig.S3J stemness_wilcoxon-celltype.exact.xlsx")

