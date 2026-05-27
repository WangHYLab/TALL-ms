library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(homologene)


setwd("normal_tall")

normal <- readRDS("SCENIC_withCTL/normal.scenic.RDS")
normal <- RunUMAP(normal, dims = 1:15, min.dist = 2)
DimPlot(normal, label = T, reduction = "umap", group.by = "celltype_normal", pt.size = 0.8)

{
    data.list <- readRDS("data.list.filtered&removeDB.rds")
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
    saveRDS(seu, "SCIN.ccregress.integrate.rds")
}
tall <- readRDS(file = "SCIN.ccregress.integrate.rds")

DimPlot(tall, label = T)
DotPlot(tall, features = c("Ms4a1", "Cd19", "Cd68", "Cd86", "Ptprc", "Cd3d", "Cd52", "S100a8", "Csf1r", "Klrc1"))
## remove none-T cell subset 15
tall <- subset(tall, seurat_clusters != "15")
DimPlot(tall, label = T)
#### 
tall$celltype <- "Tall"
normal$celltype <- as.character(normal$celltype_normal)


seu <- merge(x = tall, y = normal, add.cell.ids = c("TALL", "Normal"), project = "TALL_vs_Normal")

seu$group <- ifelse(seu$orig.ident %in% c("N_DN", "N_THY"), "Normal", "TALL")
table(seu$group) 

DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)

# ccgene
s.genes <- human2mouse(cc.genes.updated.2019$s.genes) %>%
    pull(mouseGene) %>%
    unique()
g2m.genes <- human2mouse(cc.genes.updated.2019$g2m.genes) %>%
    pull(mouseGene) %>%
    unique()

seu <- JoinLayers(seu)
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)

table(seu$Phase)
DimPlot(seu, reduction = "pca", group.by = "Phase") 

seu.bk <- seu
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$group)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)

# CCS regress
# only vf
seu.bk2 <- seu
seu <- ScaleData(seu,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = VariableFeatures(seu), verbose = FALSE
)

# backup before PCA
seu.bk3 <- seu
# PCA
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)


# select dimensions for integration based on ElbowPlot
ElbowPlot(seu, ndims = 50)
dims.use <- 1:20


## Harmony
seu.inte <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    dims = dims.use
)

# joint layers 
# seu.inte[["RNA"]] <- JoinLayers(seu.inte[["RNA"]])

# DimPlot(seu.inte, reduction = "harmony", group.by = "group") + ggtitle("Harmony整合后")
# DimPlot(seu.inte, reduction = "harmony", group.by = "Phase") + ggtitle("细胞周期校正后")

# reduce dimensionality and cluster on integrated data
seu.inte <- RunUMAP(seu.inte, reduction = "harmony", dims = dims.use, n.neighbors = 30)
DimPlot(seu.inte, group.by = "group") + ggtitle("Harmony整合后")
# 寻找邻居并聚类
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

# assign cell type based on marker expression and cluster distribution
## Bcell
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(17), ]$cell_type <- "B cell"
## DN
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(2, 16), ]$cell_type <- "DN3"
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(14), ]$cell_type <- "DN4"
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(10), ]$cell_type <- "DN1/DN2"
## DP
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(1, 7, 8, 11, 12, 13), ]$cell_type <- "DP"
## SP
seu.inte@meta.data[seu.inte@meta.data$seurat_clusters %in% c(0, 3, 4, 5, 6, 9, 15), ]$cell_type <- "CD8"

# verify annotation
table(seu.inte$cell_type, useNA = "ifany")
# # check annotation consistency with clusters
# table(seu.inte$cell_type, seu.inte$seurat_clusters)

Idents(seu.inte) <- seu.inte$cell_type

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

ggsave(plot = p1 + p2 + p3, filename = "normal&tall.cluster&sample&ct.umap.pdf", width = 15.5, height = 4.5)


## try to recluster cluster10
# seu.inte2 <- FindClusters(seu.inte, resolution = 0.8) # 生成多个分辨率供选择


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
ggsave(plot = p1 + p2 + p3, filename = "normal&tall.cluster&sample&ct.umap.noBcell.raster.pdf", width = 15.5, height = 4.5)



# table(seu.inte$seurat_clusters,seu.inte$celltype) %>% pheatmap::pheatmap(.,scale = "row")

# save result
# saveRDS(seu.inte, file = "TALL_Normal.integrated.ccregressed.rds")
saveRDS(out, file = "TALL_Normal.integrated.ccregressed.Tcell.rds")
