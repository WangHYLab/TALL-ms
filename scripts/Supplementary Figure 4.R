#########################################################################################
#### Fig.S4A
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

ggsave(plot = p2, filename = "out/Fig.S4A normal.celltype.pdf", width = 5.5, height = 5)

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
#### Fig.S4B
#########################################################################################
normal <- readRDS("normal.rds")

DotPlot(normal,
  features = c(
    "Kit", "Il7r", "Il2ra", "Cd4", "Cd8a", "Cd8b1",
    "Bcl11b", "Notch1", "Rag1", "Rag2", "Gata3", "Dntt", "Cd247", "Ptprc",
    "Cd3d", "Tcf7"
  ), cols = c("lightblue", "salmon"),
  cluster.idents = F
) +
  xlab("") + ylab("") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )
ggsave(filename = "out/Fig.S4B normal.marker.dotplot.pdf", width = 6.5, height = 3.1)



#########################################################################################
#### Fig.S4C-E
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
## remove none T cell cluster 15
tall <- subset(tall, seurat_clusters != "15")
DimPlot(tall, label = T)
#### 合并
tall$celltype <- "Tall"
normal$celltype <- as.character(normal$celltype_normal)


# merge
seu <- merge(x = tall, y = normal, add.cell.ids = c("TALL", "Normal"), project = "TALL_vs_Normal")

# add group
seu$group <- ifelse(seu$orig.ident %in% c("N_DN", "N_THY"), "Normal", "TALL")
table(seu$group) # check

# default Assay
DefaultAssay(seu) <- "RNA"


# norm and vf 
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)


# cell cycle
s.genes <- human2mouse(cc.genes.updated.2019$s.genes) %>%
    pull(mouseGene) %>%
    unique()
g2m.genes <- human2mouse(cc.genes.updated.2019$g2m.genes) %>%
    pull(mouseGene) %>%
    unique()
# cal
seu <- JoinLayers(seu)
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)

# check phase
table(seu$Phase)
DimPlot(seu, reduction = "pca", group.by = "Phase") 

seu.bk <- seu
# split to batch
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$group)

# re norm 
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

ggsave(plot = p1 + p2 + p3, filename = "out/Fig.S4C-E normal&tall.cluster&sample&ct.umap.noBcell.pdf", width = 15.5, height = 4.5)

saveRDS(out, file = "data/TALL_Normal.integrated.ccregressed.Tcell.rds")



#########################################################################################
#### Fig.S4F
#########################################################################################
out <- readRDS("TALL_Normal.integrated.ccregressed.Tcell.rds")

feat <- c("Notch1", "Il7r", "Kit", "Bcl11b", "Cd4", "Cd8a")
p <- FeaturePlot(out, features = feat, ncol = 3)
# ggsave(plot = p, filename = "normal&tall.featureplot.umap.pdf", width = 15.5, height = 9)
p <- FeaturePlot(out, features = feat, ncol = 3, raster = T)
ggsave(plot = p, filename = "out/Fig.S4F NormalAndTALL.featureplot.umap.raster.pdf", width = 15.5, height = 9)



#########################################################################################
#### Fig.S4G
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
ggsave(filename = "out/Fig.S4G.pdf", width = 15, height = 4)


#########################################################################################
#### Fig.S4H
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

ggsave(filename = "out/Fig.S4H marker.pdf", width = 4.5, height = 5)




#########################################################################################
#### Fig.S4I
#########################################################################################
source("script/funcs.R")
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(homologene)

DN_1<-process(samplename ="20221123_DN_RNA", project ="DN_1",onlyread = TRUE )
DN_2<-process(samplename ="20221127_DN_RNA", project ="DN_2",onlyread = TRUE )
Thy_1<-process(samplename ="20221123_Thy_RNA", project ="Thy_1",onlyread = TRUE )
Thy_2<-process(samplename ="20221127_Thy_RNA", project ="Thy_2",onlyread = TRUE )
DN_5<-process(samplename ="20240126_DN_RNA", project ="DN_5" ,onlyread = TRUE )
Thy_5<-process(samplename ="20240126_Thy_RNA", project ="Thy_5",onlyread = TRUE )

data<-merge(DN_1,c(DN_2,Thy_1,Thy_2,DN_5,Thy_5))
p<-VlnPlot(data,features = c('nCount_RNA','nFeature_RNA','percent.mt'),pt.size = 0,group.by = "orig.ident")
# ggsave(filename = "unfilter.vlnplot.pdf",width = 10,height = 4)

data.list<-list(DN1=DN_1,
                DN2=DN_2,
                Thy1=Thy_1,
                Thy2=Thy_2,
                DN5=DN_5,
                Thy5=Thy_5)
saveRDS(data.list, "data/data.list.filtered&removeDB.rds")

seu <- merge(data.list$DN1,c(data.list$DN2,data.list$Thy1,data.list$Thy2))

DefaultAssay(seu)<-"RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst",nfeatures=3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims=1:20)

# Calculating Cell cycle score
seu <- JoinLayers(seu)
seu <- CellCycleScoring(seu,s.features=s.gene,g2m.features = g2m.gene)

# Regress out cell cycle score
# seu$CC.Difference <- seu$S.Score - seu$G2M.Score
# seu <- ScaleData(seu, vars.to.regress = "CC.Difference", features=VariableFeatures(seu))
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features=VariableFeatures(seu))
seu <- RunPCA(seu, reduction.name = "pca.ccregress")

# Perform RPCA integration
seu[["RNA"]] <- split(seu[["RNA"]], f=seu$orig.ident)

# options(future.globals.maxSize = 8000 * 1024^2)
seu <- IntegrateLayers(object = seu,
                       method = CCAIntegration,
                       dims=1:20,
                       orig.reduction = "pca.ccregress",
                       new.reduction = "integrated.cca",
                       verbose = FALSE)

# re-join layers after integration
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
seu <- RunUMAP(seu, dims = 1:20, reduction = "integrated.cca")
seu <- FindNeighbors(seu,dims = 1:20,reduction = "integrated.cca")
seu <- FindClusters(seu,resolution = 0.8)
DimPlot(seu)
saveRDS(seu, "data/TALL.ccregress.integrate.rds")


### Use TCR to select tall cell
data <- readRDS("data/TALL.ccregress.integrate.rds")
data$barcode<-rownames(data@meta.data)


data.tcr <- readRDS("data/data.tcr.rds")
DN_1<-data.tcr$DN_1
DN_2<-data.tcr$DN_2
Thy_1<-data.tcr$Thy_1
Thy_2<-data.tcr$Thy_2
DN_1$barcode<-paste0(DN_1$barcode,"_1")
DN_1$sample<-"DN_1"
DN_2$barcode<-paste0(DN_2$barcode,"_2")
DN_2$sample<-"DN_2"
Thy_1$barcode<-paste0(Thy_1$barcode,"_3")
Thy_1$sample<-"Thy_1"
Thy_2$barcode<-paste0(Thy_2$barcode,"_4")
Thy_2$sample<-"Thy_2"
# 所有TCR 数据
TCR<-rbind(DN_1,DN_2,Thy_1,Thy_2)

overlap<-intersect(rownames(data@meta.data),TCR$barcode)
TCR.use<-TCR%>%filter(barcode %in% overlap)

data$withTCR<-"none"
data@meta.data[data@meta.data$barcode%in%overlap,]$withTCR<-"yes"
DimPlot(data,group.by = "withTCR")+ggtitle(label = paste0("cell with TCR(",round(length(overlap)/nrow(data@meta.data),4)*100,"%)"))
# ggsave(filename = "./out/data.4.withTCR.pdf",width = 5,height = 4.6)

data.use<-subset(data,withTCR=="yes")
overlap<-intersect(rownames(data.use@meta.data),TCR$barcode)
TCR.use<-TCR%>%filter(barcode %in% overlap)


freq.tcr<-TCR.use$cdr3%>%table()%>%data.frame()
cdr3_1<-'CASSDAAGARDTQYF' #DN1 THY1
cdr3_2<-'CAWSLTGGARSQNTLYF' #DN2 THY2
bc1<-TCR.use[TCR.use$cdr3==cdr3_1,]$barcode%>%unique()
bc2<-TCR.use[TCR.use$cdr3==cdr3_2,]$barcode%>%unique()

data.use$cdr3_clone<-"other"
data.use@meta.data[data.use@meta.data$barcode%in%bc1,]$cdr3_clone<-"clone1"
data.use@meta.data[data.use@meta.data$barcode%in%bc2,]$cdr3_clone<-"clone2"
DimPlot(data.use,group.by = "cdr3_clone",split.by = "orig.ident")+ggtitle(label ="cdr3 clone",
                                                                          subtitle = paste0("clone1(",round(length(bc1)/nrow(data.use@meta.data),4)*100,"%): ",cdr3_1,
                                                                                            " (DN1/Thy1)\nclone2(",round(length(bc2)/nrow(data.use@meta.data),4)*100,"%): ",cdr3_2," (DN2/Thy2)"))

table(data.use$cdr3_clone,data.use$orig.ident) %>% heatmap()

seu <- RunUMAP(data, dims = 1:20, reduction = "integrated.cca")
seu <- FindNeighbors(seu,dims = 1:20,reduction = "integrated.cca")
seu <- FindClusters(seu,resolution = 0.8)
p<-FeaturePlot(seu,features = c("Ptprc","Cd3d","Kit","Cd4",
                             "Cd8a","Il2ra","Notch1","Bcl11a",
                             "Bcl11b","Runx1","Runx3","Mki67"),cols = c("grey","red"))
# ggsave(plot=p,filename = "Clone.features.pdf",width = 14,height = 9)

## celltype annotation ---
library(celldex)
library(SingleR)
seu<-annoBySingler(seu)
seu@active.ident<-seu$seurat_clusters

## check cell marker 
# feat<-c(
#   # "Itgam","Itgax", #Myeloid
#   # "Ms4a1","Cd79a",'Iglc2',#'Iglc3',#'Ctsh', #B cell
#   "Kit","Il2ra",
#   'Hes1','Mpzl2', #DN
#   # 'Tnfrsf4','Tnfrsf9','Foxp3',# T reg
#   "Cd3d","Ptprc","Lef1",
#   "Cd4","Cd8a", #DP
#   "Mki67", # profi
#   "Bcl11a","Bcl11b",
#   "Notch1", # T-ALL
#   'Tcrg-C2','Cpa3','Trdc'# Tgd
# )
# DotPlot(seu,cluster.idents = T,features =  feat,assay = "RNA")+coord_flip()

## rename ident
seu<-RenameIdents(seu,
                  "11"="DN2 (11)",  #Kit高, Il2ra高，略带Cd8a
                  "12"="DN4 (12)",  #Kit低, Il2ra低，略带Cd8a
                  "3"="DP (3)",
                  "0"="CD8+ (other)",
                  "1"="CD8+ (other)",
                  "2"="CD8+ (other)",
                  "4"="CD8+ (other)",
                  "5"="CD8+ (other)",
                  "6"="CD8+ (other)",
                  "7"="CD8+ (other)",
                  "8"="CD8+ (other)",
                  "9"="CD8+ (other)",
                  "10"="CD8+ (other)"
                  )
saveRDS(seu, "data/clone.inte.ccr.rds")

### TCR data
DN_1_contig<-read.csv("cellranger/20221123_DN_TCR/outs/filtered_contig_annotations.csv")
DN_2_contig<-read.csv("cellranger/20221127_DN_TCR/outs/filtered_contig_annotations.csv")
Thy_1_contig<-read.csv("cellranger/20221123_Thy_TCR/outs/filtered_contig_annotations.csv")
Thy_2_contig<-read.csv("cellranger/20221127_Thy_TCR/outs/filtered_contig_annotations.csv")
data.tcr<-list(DN_1=DN_1_contig,DN_2=DN_2_contig,Thy_1=Thy_1_contig,Thy_2=Thy_2_contig)
saveRDS(data.tcr, "data/data.tcr.rds")

### run slingshot
library(slingshot)
library(RColorBrewer)
library(paletteer) 
library(viridis)
library(Seurat)
library(tradeSeq)
library(tidyverse)
library(grDevices)
library(ggthemes)
library(Seurat)
data<-readRDS("data/clone.inte.ccr.rds")


sce <- as.SingleCellExperiment(data,assay = "RNA")
sce <- slingshot(sce, clusterLabels = 'ident', reducedDim = 'UMAP', start.clus= "DN2 (11)",end.clus=c("CD8+ (other)"))
saveRDS(sce, "data/clone.slingshot.rds")

#设置颜色
cell_pal <- function(cell_vars, colors,...) {
  categories <- sort(unique(cell_vars))
  pal <- setNames(colors, categories)
  return(pal[cell_vars])
}

cell_colors <- cell_pal(sce$ident, scales::hue_pal()(9))
pdf(file = "out/Fig.S4I clone.slngshot.pdf", width = 5, height = 5)
plot(reducedDims(sce)$UMAP, col = cell_colors,pch=16, asp = 1, cex = 0.4)
lines(SlingshotDataSet(sce), lwd=2, col='black',type = "lineages")

#计算celltype坐标位置，用于图中标记
celltype_label <- data@reductions$umap@cell.embeddings%>% as.data.frame() %>%cbind(celltype = data@meta.data$ident) %>%
  group_by(celltype) %>%summarise(UMAP1 = median(umap_1),UMAP2 = median(umap_2))
index<-celltype_label
for (i in 1:nrow(index)) {text(index$celltype[i], x=index$UMAP1[i]-1, y=index$UMAP2[i])}
dev.off()


#########################################################################################
#### Fig.S4J
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

# use ori umap
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

# use ori umap
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



### PLOT
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

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "out/Fig.S4J S1.clone.monocle3.pdf",width = 18,height = 5)

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

ggsave(plot = cowplot::plot_grid(p1,p2,p3,ncol = 3),filename = "out/Fig.S4J S2.clone.monocle3.pdf",width = 18,height = 5)



#########################################################################################
#### Fig.S4K
#########################################################################################
### extracting the gene expression matrix from the Seurat
library(Seurat)
library(tidyverse)

### load data
seu <- readRDS("data/clone.inte.ccr.rds") ## scripts\Figure 2.R  in Fig.2E
# DimPlot(seu)
ctl <- readRDS("data/data.ctl.assigned.rds") ## scripts\utils\Normal-lineage.R
# DimPlot(ctl, group.by = "orig.ident")

ctl_use <- subset(ctl, orig.ident %in% c("DN", "THY"))
# DimPlot(ctl_use,label = T)
# DimPlot(ctl_use,group.by = "orig.ident",label = T)
# filter the main cell type for better visualization
main_cell <- table(ctl_use$ident) %>%
  sort() %>%
  data.frame() %>%
  filter(Freq > 50) %>%
  pull(Var1) %>%
  as.character()

ctl_use_main <- subset(ctl_use, ident %in% main_cell)
# DotPlot(ctl_use_main,features = c('Ptprc',"Cd3d","Cd3g","Cd3e","Cd4",
#                              "Cd8a","Ms4a1","Cd68","Cd34","Kit","Hes1","Il2ra"),cluster.idents = T)

normal <- ctl_use_main
normal <- RunUMAP(normal, dims = 1:30)
# DimPlot(normal,label = T,group.by = "orig.ident")
normal$orig.ident <- paste0("N_", as.character(normal$orig.ident))

seu@assays$SCT <- NULL
DefaultAssay(seu) <- "RNA"
normal@assays$SCT <- NULL
DefaultAssay(normal) <- "RNA"
seu@reductions
normal@reductions

### merge data
data <- merge(seu, normal)
data <- JoinLayers(data, assay = "RNA", layers = "counts")
saveRDS(data, "data/merge_ctl&clone.rds")
# for pyscenic
count <- LayerData(data, layer = "count")
count <- t(as.matrix(count))

write.csv(count, file = "data/for.scenic.count.ctl-tall.csv")

### Run pyscenic:  scripts\utils\pySCENIC.sh


### plot
## cell assign check for merge_ctl&clone.rds
data <- readRDS("data/merge_ctl&clone.rds")
{
  data <- readRDS("data/merge_ctl&clone.rds")

  DimPlot(data)
  normal <- subset(data, orig.ident %in% c("N_DN", "N_THY"))
  seu <- NormalizeData(normal)
  seu <- FindVariableFeatures(seu, selection.method = "vst")
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- RunUMAP(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 1.0)
  seu@active.ident

  DotPlot(seu,
    features = c(
      "Ptprc", "Cd3d", "Cd3g", "Cd3e", "Cd4", "Cd8a", "Cd8b1",
      "Kit", # "Hlf","Hoxa9",
      "Il2ra", "Cd44", "Cd24a",
      "Lef1", "Bcl11a", "Bcl11b", "Notch1",
      "Nkg7", "Klra1", "Klrb1a", "Ncam1"
    ),
    group.by = "seurat_clusters", cluster.idents = T
  ) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))


  seu@active.ident <- seu$seurat_clusters

  seu <- RenameIdents(seu,
    "16" = "ETP",
    "10" = "DN2",
    "6" = "DN2",
    "21" = "DN2",
    "20" = "DN3",
    "15" = "DN3",
    "1" = "DN3",
    "3" = "DN3",
    "18" = "DN3",
    "17" = "DN3",
    "13" = "DN4",
    "9" = "DN4",
    "12" = "DP",
    "14" = "DP",
    "19" = "DP",
    "4" = "DP",
    "22" = "DP",
    "11" = "DP",
    "0" = "DP",
    "2" = "DP",
    # "9"="DP",

    "8" = "CD8+",
    "5" = "CD8+",
    "7" = "CD8+"
  )
  # cluster17 isolation but expression of maker support it is DN3 cell
  seu$celltype_normal <- seu@active.ident
  normal <- seu

  DotPlot(normal,
    features = c(
      "Ptprc", "Cd3d", # "Cd3g","Cd3e",
      "Kit", # "Hlf","Hoxa9",

      "Il2ra", "Bcl11a", "Cd28",
      # "Cd44",
      "Cd24a",
      "Cd4", "Cd8a", "Cd8b1",
      "Lef1", "Bcl11b", "Notch1"
    ), cols = c("lightblue", "red"),
    cluster.idents = F
  ) +
    coord_flip() + xlab("") + ylab("") + theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black")
    )

  # ggsave("Normal.celltype.marker.2.pdf",width = 5,height = 4.5)

  saveRDS(seu, "data/normal.scenic.RDS")
} # code to get normal.scenic.RDS, which is to run SCENIC in normal samples(DN & Thy)

#### SCENIC regulon activity plot
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
# library(scRNAseq)
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize)
library(RcisTarget)

## 1.extract out_SCENIC.loom data
loom <- open_loom("data/aucell_output.loom")

regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons_incidMat[1:4, 1:4]
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)
close_loom(loom)

rownames(regulonAUC)
names(regulons)

## 提取细胞类型
data <- readRDS("merge_ctl&clone.rds")
normal <- readRDS("normal.scenic.RDS")

data$celltype <- as.character(data$ident)
data@meta.data[normal@meta.data %>% rownames(), ]$celltype <- paste0("N-", normal@active.ident %>% as.character())

data$celltype <- factor(data$celltype, levels = c(
  "N-ETP", "N-DN2", "N-DN3", "N-DN4", "N-DP", "N-CD8+",
  "DN2 (11)", "DN4 (12)", "DP (3)", "CD8+ (other)"
))

sub_regulonAUC <- regulonAUC[, match(colnames(data), colnames(regulonAUC))]


cellClusters <- data.frame(
  row.names = colnames(data),
  seurat_clusters = as.character(data$seurat_clusters)
)
cellTypes <- data.frame(
  row.names = colnames(data),
  celltype = data$celltype
)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4, 1:4]

# save
save(sub_regulonAUC, cellTypes, cellClusters, data,
  file = "data/for_rss_and_visual.Rdata"
)
load("data/for_rss_and_visual.Rdata")


### 画图

selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(
  rownames(cellTypes),
  cellTypes[, selectedResolution]
)

# remove extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(
  cellsPerGroup,
  function(cells) {
    rowMeans(getAUC(sub_regulonAUC)[, cells])
  }
)

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
  center = T, scale = T
))
out <- data.frame(regulonActivity_byGroup_Scaled, check.names = F, check.rows = F) %>% tibble::rownames_to_column("TF")

# write.xlsx(out,file = "regulonActivity.byGroup.Scaled.0609.xlsx")

dim(regulonActivity_byGroup_Scaled)
# [1] 209   9
plot <- na.omit(regulonActivity_byGroup_Scaled)

# pdf("~/Project/YuYong_TALL/out/SCENIC_ourData/Fig.1O TF.all.pdf", width = 5.2, height = 20)
Heatmap(
  plot,
  name = "z-score",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  # row_title_rot                = 90,
  column_names_rot = -90,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE
)

# dev.off()



## rss TF
rss <- calcRSS(
  AUC = getAUC(sub_regulonAUC),
  cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution]
)
rss <- na.omit(rss)
rssPlot <- plotRSS(rss)

rss_high_tf <- rssPlot$df %>%
  group_by(cellType) %>%
  arrange(desc(RSS)) %>%
  slice_head(n = 5) %>%
  pull(Topic) %>%
  as.character() %>%
  unique()

p1 <- plotly::ggplotly(rssPlot$plot)
# ggsave(filename = "SCENIC.TF.RSS.celltype-dotplot.0609.pdf",width = 3.5,height = 15)

rss_high_tf
rssPlot <- plotRSS(rss[rss_high_tf, ])
p2 <- plotly::ggplotly(rssPlot$plot)
# ggsave(filename = "SCENIC.TF.RSS.top10.celltype-dotplot.0609.pdf",width = 3.5,height = 6)


score <- regulonActivity_byGroup_Scaled
head(score)
df <- do.call(
  rbind,
  lapply(1:ncol(score), function(i) {
    dat <- data.frame(
      path = rownames(score),
      cluster = colnames(score)[i],
      sd.1 = score[, i],
      sd.2 = apply(score[, -i], 1, median)
    )
  })
)
df$fc <- df$sd.1 - df$sd.2
top5 <- df %>%
  group_by(cluster) %>%
  top_n(5, fc)

n <- score[c(top5$path) %>% unique(), ]

## show selected TFs
regulonsToPlot <- c(
  "Arid5b(+)", "Nfkb2(+)", "Hmga2(+)", "Tcf7l1(+)", "Twist1(+)", "Meis1(+)", "Spi1(+)",
  "Taf6(+)", "Erg(+)", "Nfe2(+)", "Jun(+)", "Pbx1(+)", "Tfap4(+)",
  "Myc(+)", "Runx1(+)",
  "Stat1(+)", "Maf(+)", "Runx3(+)", "Relb(+)",
  "Ets2(+)", "Bcl6(+)", "Klf13(+)", "Ets2(+)", "Jund(+)",
  "E2f1(+)", "Tfdp1(+)", "Hdac2(+)",
  "E2f8(+)", "Vdr(+)",
  "Runx1(+)", "Bcl11b(+)", "Myb(+)", "Myc(+)", "Dnmt1(+)", "E2f7(+)", "Tcf12(+)", "Tcf3(+)", "Ikzf1(+)", "Hoxa10(+)", "Ets1(+)",
  "Smad4(+)", "Tef(+)",
  "Gfi1b(+)", "Gata2(+)", "Tal1(+)", "Tbx21(+)",
  "Ikzf3(+)", "Emx2(+)", "Dpf1(+)", "Cebpa(+)",
  "Atf5(+)"
) %>% unique()

ht <- plot_heatmap(regulonsToPlot)
plot_heatmap <- function(regulonsToPlot) {
  regulonsToPlot <- regulonsToPlot[regulonsToPlot %in% row.names(sub_regulonAUC)]

  top2 <- df %>%
    group_by(cluster) %>%
    top_n(3, fc)
  rank <- df %>%
    group_by(cluster) %>%
    top_n(100, fc)
  n <- score[c(regulonsToPlot) %>% unique(), ]

  genes <- unique(c(regulonsToPlot, top2$path))
  genes <- unique(c(regulonsToPlot))
  existing_genes <- genes[genes %in% rownames(n)]

  at_pos <- which(rownames(n) %in% existing_genes)
  labels_vec <- rownames(n)[at_pos] # 按热图实际行顺序

  # pheatmap::pheatmap(n,cluster_rows = F,cluster_cols = F)

  # plot heatmap
  ht <- Heatmap(
    n,
    name = "Z-score",
    cluster_rows = F,
    cluster_columns = F,
    clustering_method_columns = "ward.D",
    show_column_names = TRUE,
    show_row_names = TRUE,
    col = colorRamp2(c(-2, 0, 2), c("#1A4d73", "white", "#c12c1f")),
    row_gap = unit(1.5, "mm"),
    column_gap = unit(1.5, "mm"),
    border = TRUE,
    top_annotation = columnAnnotation(
      df = data.frame(Group = c(rep("Normal", 6), rep("T-ALL", 4))),
      col = list(Group = c("T-ALL" = "#F69389", "Normal" = "lightblue"))
    ),
  )


  ht_opt$message <- FALSE
  draw(ht)
  return(ht)
}

pdf("out/Fig.S4K TF regulon.pdf", width = 5.2, height = 8)
draw(ht)
dev.off()


#########################################################################################
#### Fig.S4L
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
ggsave(filename = "out/Fig.S4L TALL&Normal.Stemness.inte.ccr.new.pdf", width = 5, height = 4)


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
                     file = "out/Fig.S4L stemness_wilcoxon-celltype.exact.xlsx")
