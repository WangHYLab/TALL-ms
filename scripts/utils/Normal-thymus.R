###############################################################################################
# analysis for normal thymus data
###############################################################################################
normal <- readRDS("normal.scenic.RDS")

## normal cell plot
normal <- RunUMAP(normal, dims = 1:15, min.dist = 2)
p1 <- DimPlot(normal, label = T, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.8)
p2 <- DimPlot(normal, label = T, reduction = "umap", group.by = "celltype_normal", pt.size = 0.8)

ggsave(plot = p2, filename = "normal.celltype.pdf", width = 5.5, height = 5)
DotPlot(normal,
    features = c(
        # "Ptprc","Cd3d",#"Cd3g","Cd3e",
        # "Kit",#"Hlf","Hoxa9",
        #
        # "Il2ra","Bcl11a","Cd28",
        # # "Cd44",
        # "Cd24a",
        # "Cd4","Cd8a","Cd8b1",
        #
        # "Lef1","Bcl11b","Notch1"

        # "Hhex","Cd44",
        "Kit", "Il7r", "Il2ra", "Cd4", "Cd8a", "Cd8b1",
        # "Spi1", "Cd44",
        "Bcl11a", "Bcl11b", "Notch1",
        "Rag1", "Rag2", "Gata3", "Dntt", "Cd247", "Ptprc",
        # "Cd52",
        "Cd3d",
        # "Cd3e",
        "Tcf7"
    ),
    cluster.idents = F
) +
    xlab("") + ylab("") + theme_bw() +
    theme(
        axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
    )

# cluster17 isolation but expression of maker support it is DN3 cell
ggsave(filename = "normal.marker.dotplot.pdf", width = 6.5, height = 3.1)
saveRDS(normal, "normal.rds")


### monocle3 analysis
library(Seurat)
library(monocle3)
library(ggplot2)
library(scales)
library(ggthemes)
library(dplyr)
library(ggrastr)
theme_blank <- theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = margin(1, 1, 1, 1, unit = "cm"))


data <- normal
cell_metadata <- data@meta.data
mtx <- GetAssayData(data, assay = "RNA", slot = "counts") %>% as.matrix()
gene_annotation <- data.frame(gene_short_name = rownames(mtx))
int.embed <- Embeddings(data, reduction = "umap")
rownames(gene_annotation) <- rownames(mtx)

save(cell_metadata, mtx, gene_annotation, int.embed, file = "toDELL.monocle3-need.S5.rdata")

## Run in another server
load("toDELL.monocle3-need.S5.rdata")
# load("cds.rdata")

cds <- new_cell_data_set(mtx,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 15)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)

# plot_cells(cds, color_cells_by="partition")

# use Seurat UMAP embedding to replace monocle3's default embedding, for better visualization
cds.embed <- cds@int_colData$reducedDims$UMAP

int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed
# p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype_normal") + ggtitle('int.umap')


# cds <- cluster_cells(cds)
cds <- cluster_cells(cds,resolution = 1e-5, partition_qval = 0 )
# plot_cells(cds, color_cells_by = "cluster")
cds <- learn_graph(cds, learn_graph_control = list(minimal_branch_len = 2))
# plot_cells(cds, color_cells_by = "cluster")

# save(cds,file = "cds.rdata")
cds <- order_cells(cds)


p1 <- plot_cells(cds,
    reduction_method = "UMAP", color_cells_by = "celltype_normal", rasterize = T,
    cell_size = 1.5, trajectory_graph_color = "black", label_cell_groups = F,
    graph_label_size = 4, label_leaves = F, label_branch_points = F
) + theme_bw() +
    scale_color_manual(values = hue_pal()(9)) + theme_blank

p2 <- plot_cells(cds,
    reduction_method = "UMAP", color_cells_by = "pseudotime", label_roots = T, rasterize = T,
    cell_size = 1.5, trajectory_graph_color = "black",
    graph_label_size = 4
) + theme_bw() + theme_blank

ggsave(plot = cowplot::plot_grid(p1, p2, ncol = 2), filename = "normal.monocle3.raster.pdf", width = 12, height = 5)




###############################################################################################
# analysis for merge tall-clone and normal
###############################################################################################
normal <- readRDS("~/Project/YuYong_TALL/out/SCENIC_withCTL/normal.scenic.RDS")
normal <- RunUMAP(normal, dims = 1:15, min.dist = 2)
DimPlot(normal, label = T, reduction = "umap", group.by = "celltype_normal", pt.size = 0.8)



tall <- readRDS("~/Project/YuYong_TALL/clone.inte.ccr.rds")
DimPlot(tall)
tall$celltype <- as.character(tall@active.ident)
normal$celltype <- as.character(normal$celltype_normal)

data <- merge(tall, normal)

data@active.assay <- "RNA"

seu <- data
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
# seu <- RunUMAP(seu, dims=1:20)

# Calculating Cell cycle score
seu <- JoinLayers(seu)
seu$group <- "TALL"
seu@meta.data[seu@meta.data$orig.ident %in% c("N_DN", "N_THY"), ]$group <- "Normal"
# Regress out cell cycle score
s.gene <- human2mouse(cc.genes.updated.2019$s.genes) %>%
    pull(mouseGene) %>%
    unique()
g2m.gene <- human2mouse(cc.genes.updated.2019$g2m.genes) %>%
    pull(mouseGene) %>%
    unique()
seu <- CellCycleScoring(seu, s.features = s.gene, g2m.features = g2m.gene)

seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(seu))
seu <- RunPCA(seu, reduction.name = "pca.ccregress")

# Perform harmony integration
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$group)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(seu))
seu <- RunPCA(seu, reduction.name = "pca.ccregress")

# options(future.globals.maxSize = 8000 * 1024^2)
seu.inte <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    dims = 1:20,
    orig.reduction = "pca.ccregress",
    new.reduction = "harmony",
    verbose = FALSE
)

# re-join layers after integration
# seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
seu.inte <- RunUMAP(seu.inte, dims = 1:20, )
seu.inte <- FindNeighbors(seu.inte, dims = 1:20, reduction = "harmony")
seu.inte <- FindClusters(seu.inte, resolution = 0.8)

seu.inte$celltype <- factor(seu.inte$celltype, levels = c(
    "ETP", "DN2", "DN3", "DN4", "DP", "CD8+",
    "DN2 (11)", "DN4 (12)", "DP (3)", "CD8+ (other)"
))
p1 <- DimPlot(
    object = seu.inte,
    reduction = "umap",
    group.by = "celltype",
    label = T, pt.size = 0.5
)

p2 <- DimPlot(
    object = seu.inte,
    reduction = "umap",
    group.by = "group",
    label = T, pt.size = 0.5
)
ggsave(plot = p1 + p2, filename = "normal&tall.celltype.pdf", width = 13, height = 5)



saveRDS(seu.inte, "normal&CLONE.inte-harmony.ccr.rds")

data <- seu.inte


DotPlot(seu.inte,
    features = c(
        "Ptprc", "Cd3d", # "Cd3g","Cd3e",
        "Kit", # "Hlf","Hoxa9",

        "Il2ra", "Bcl11a", "Cd28",
        # "Cd44",
        "Cd24a",
        "Cd4", "Cd8a", "Cd8b1",
        "Lef1", "Bcl11b", "Notch1"
    ), cols = c("lightblue", "red"),
    cluster.idents = F, group.by = "celltype"
) +
    coord_flip() + xlab("") + ylab("") + theme_bw() +
    theme(
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")
    )
ggsave(filename = "normal&tall.celltype.marker.pdf", width = 5, height = 4.5)




library(openxlsx)
library(homologene)
library(dplyr)
library(ggsignif)

gene <- read.xlsx("/home/zhengjie/Project/YuYong_TALL/out/41467_2025_61222_MOESM5_ESM.xlsx") %>%
    filter(cluster == "Stem_like") %>%
    pull(gene) %>%
    human2mouse(.) %>%
    pull(mouseGene) %>%
    unique()

data[["RNA"]] <- JoinLayers(data[["RNA"]])

data <- AddModuleScore(data, features = list(stemness = gene), name = "stemness", assay = "RNA")
data$stemness <- data$stemness1
data$stemness1 <- NULL
VlnPlot(data, features = "stemness", pt.size = 0)

df <- data@meta.data %>%
    dplyr::select(celltype, stemness) %>%
    mutate(ident = celltype)
ggplot(df, aes(x = celltype, y = stemness)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    # geom_jitter(aes(x=ident,y=stemness),width = 0.1,size=0.8,alpha=0.1)+
    geom_boxplot(outlier.alpha = 0, width = 0.6, aes(fill = celltype)) +
    # geom_signif(comparisons = list(c("DN2 (11)", "DN4 (12)"),
    #                                c("DN2 (11)", "DP (3)"),
    #                                c("DN2 (11)", "CD8+ (other)")),
    #             y_position=c(0.28, 0.26, 0.24),
    #             tip_length = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
    #             map_signif_level=TRUE)+

    # geom_boxplot(width=0.6,outlier.size = 0,outlier.alpha = 0)+
    theme_classic() +
    ylab("Stemness score") +
    xlab("") +
    theme(
        axis.text.x = element_text(size = 12, angle = -45, hjust = 0),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    NoLegend()
# scale_fill_manual(values = RColorBrewer::brewer.pal(4,"Set2"))
ggsave(filename = "TALL&Normal.Stemness.inte.ccr.pdf", width = 5, height = 4)
