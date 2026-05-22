########################################################################################################
#### Ellen data analysis - ScienceImmunology2022
########################################################################################################
setwd("ref_data/GSE183026_ScienceImmunology2022/out")
library(Seurat)
library(tidyverse)
library(homologene)
library(ggplot2)
library(SingleR)
library(celldex)
library(patchwork)

get_data(){
  ## ATO data
  data<-Read10X_h5("ref_data/GSE183026_ScienceImmunology2022/GSE183026_ATO_filtered_gene_bc_matrices_h5.h5")
  seu<-CreateSeuratObject(data,min.cells = 3,min.features = 1200,project = "ATO")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  
  VlnPlot(seu,features =c("nCount_RNA","nFeature_RNA","percent.mt") ,pt.size = 0)
  seu<-subset(seu,nCount_RNA<27000&nFeature_RNA<4400&percent.mt<9)
  ATO<-seu
  
  ## Thy data
  data2<-Read10X("ref_data/GSE183026_ScienceImmunology2022/GSM4072329/")
  meta<-read.csv("ref_data/GSE183026_ScienceImmunology2022/GSE137165_meta.csv")
  tag<-read.csv("ref_data/GSE183026_ScienceImmunology2022/GSM4072329/tags.csv",row.names = 1)
  table(meta$sample_id)
  seu<-CreateSeuratObject(data2,min.cells = 3,min.features = 1300,project = "Thy")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  
  cell.use<-intersect(meta$cell_barcode,tag %>% rownames())
  seu[["HTO"]] <- CreateAssayObject(counts = t(x = tag[colnames(seu), 1:5]))
  
  # Normalize data
  seu <- NormalizeData(seu)
  seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
  # Demultiplex data
  seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)
  
  # RidgePlot(seu, assay = "HTO", features = c("subpop-1", "subpop-2"), ncol = 2)
  # HTOHeatmap(seu, assay = "HTO")
  # FeatureScatter(seu, feature1 = "subpop-1", feature2 = "subpop-CTRL")
  Thy <- subset(seu, nFeature_RNA>6800|hash.ID =="subpop-CTRL", invert = TRUE)
  Thy <- subset(Thy, idents = "Doublet", invert = TRUE)
  
  
  data<-merge(ATO,Thy)
  table(data$orig.ident)
  save(Thy,ATO,file = "Ellen.data.rdata")
  saveRDS(data,"Ellen.data.rds")
}
# get_data()
s.gene   <- human2mouse(cc.genes.updated.2019$s.genes)   %>% pull(mouseGene) %>% unique()
g2m.gene <- human2mouse(cc.genes.updated.2019$g2m.genes) %>% pull(mouseGene) %>% unique()

## integrate ----
load("Ellen.data.rdata")
data<-readRDS("Ellen.data.rds")
data <- NormalizeData(data)
data <- FindVariableFeatures(data,nfeatures = 3000)
# data<-CellCycleScoring(data,
#                  s.features   = human2mouse(cc.genes.updated.2019$s.genes)   %>% pull(mouseGene) %>% unique(),
#                  g2m.features = human2mouse(cc.genes.updated.2019$g2m.genes) %>% pull(mouseGene) %>% unique() )


data <- ScaleData(data)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap",group.by = "orig.ident")
ggsave(filename = "1.umap.uncellcycle.unintegrate.pdf",width = 5.5,height =5)



data2 <- IntegrateLayers(object = data, method = CCAIntegration,dims=1:20, verbose = F,orig.reduction = "pca",new.reduction = "cca1")
data2 <- RunUMAP(data2, dims = 1:30,reduction = "cca1",reduction.name = "cca1_umap")
DimPlot(data2, reduction = "cca1_umap",group.by = "orig.ident")
ggsave(filename = "2.umap.uncellcycle.integrate.pdf",width = 5.5,height =5)

## new wf :https://github.com/satijalab/seurat/issues/8273
load("Ellen.data.rdata")
seu <- merge(ATO,Thy)
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
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features=rownames(seu))
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
seu <- RunUMAP(seu, dims = 1:20, reduction = "integrated.cca",
               n.neighbors = 50,
               # dens.lambda = 3,n.epochs = 300,
               min.dist = 0.02)
seu <- FindNeighbors(seu,dims = 1:20,reduction = "integrated.cca")
seu <- FindClusters(seu,resolution = 0.8)
DimPlot(seu)

FeaturePlot(seu,features = c("Cd34","Kit","Hlf","Gcnt2","Il2ra","Bcl11b","Lef1"))
saveRDS(seu,"data.ccregress.integrate.rds")
saveRDS(seu,"data.ccregress-allgene.integrate.rds")

#### 注释和筛选 ----
data<-readRDS("data.ccregress.integrate.rds")
annoBySingler<-function(data,fineOrMain="fine",addClusterToAnno=TRUE,check.cluster=TRUE){
  ref.se<-ImmGenData()
  # matrix <- GetAssayData(data, slot = 'data')
  matrix <- LayerData(data,layer = "data")

  singler.cluster <-SingleR(matrix,ref = ref.se,
                            labels = ref.se$label.fine,
                            clusters = data@meta.data$seurat_clusters)

  celltype = data.frame(ClusterID = rownames(singler.cluster),celltype = singler.cluster$labels,stringsAsFactors = F)
  ## load to meta.data
  data@meta.data$celltype = 'NA'
  if(addClusterToAnno){
    for (i in 1:nrow(celltype)) {
      if(check.cluster){
        data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- paste0("c",i,"_",celltype$celltype[i])
      }else{
        data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- paste0("c",i-1,"_",celltype$celltype[i])
      }
      
    }
  }else{
    for (i in 1:nrow(celltype)) {
      data@meta.data[which(data@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]
    }
  }
  
  return(data)
}

# # 细胞注释：  
# data.2<-RunUMAP(data, dims = 1:20, reduction = "integrated.cca",
#               n.neighbors = 50,n.components = 5,
#               # dens.lambda = 3,n.epochs = 300,
#               min.dist = 0.02)

data<-annoBySingler(data)
data.1<-FindClusters(data,resolution = 2.5)
DimPlot(data.1,label = T)
data2<-subset(data.1,seurat_clusters!=24&seurat_clusters!=7&seurat_clusters!=9&seurat_clusters!=23&seurat_clusters!=25)
DimPlot(data2,label = T)
data2 <- RunUMAP(data2, dims = 1:20, reduction = "integrated.cca",
               n.neighbors = 50,
               # dens.lambda = 3,n.epochs = 300,
               min.dist = 0.02)
data2 <- FindNeighbors(data2,dims = 1:20,reduction = "integrated.cca")
data2 <- FindClusters(data2,resolution = 1.0)
DimPlot(data2,label = T,repel = T)
data2<-annoBySingler(data2)
p1<-DimPlot(data2,label = T,repel = T,group.by = "celltype",pt.size = 0.8)
p2<-DimPlot(data2,label = T,repel = T,group.by = "orig.ident",pt.size = 0.8)+ggtitle("")
p3<-DimPlot(data2,label = T,repel = T,pt.size = 0.8)+ggtitle("")
pdf("Ellen.umap.out.pdf",width = 18,height = 5)
p2|p3|p1
dev.off()


FeaturePlot(data2,features = c("Kit","Il2ra","Bcl11b","Lef1"),ncol = 4)
ggsave(filename = "Ellen.umap.features.pdf",width = 13,height = 3)
data2<-RenameIdents(data2,
                    "11"='ETP',
                    "0"="ETP",
                    "1"="ETP",
                    "2"="ETP",
                    "3"="ETP",
                    "4"="DN2",
                    "5"="DN2",
                    "7"="DN2",
                    "8"="DN2",
                    "9"="DN2",
                    "6"="DN3",
                    "10"="DN3"
                    )
DimPlot(data2)
DotPlot(data2,features = c("Ptprc","Kit","Bcl11a","Il2ra","Bcl11b","Lef1","Cd3g","Cd4","Cd8a"))

saveRDS(data2, "for_SCIN.rds")

########################################################################################################
#### Scin expression with Ellen
########################################################################################################
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(homologene)

data.ellen <- readRDS("ref_data/GSE183026_ScienceImmunology2022/out/for_SCIN.rds") # which is finished in utils/Ellen-data.R
data.ellen$orig.ident <- paste0("Ellen.", data.ellen$orig.ident)
data.ellen$ori_celltype <- data.ellen@active.ident


data.tlineage <- readRDS("../data.Tlineage.assign.rds") # which is finished in scripts/utils/Normal-lineage.R
data.tlineage$ori_celltype <- data.tlineage@active.ident


data.tlineage$Scin_scaled <- data.tlineage@assays$SCT@scale.data["Scin", ]
ggplot(data.tlineage@meta.data, aes(x = ori_celltype, y = Scin_scaled)) +
  geom_boxplot()

VlnPlot(data.ellen, features = "Scin", assay = "RNA")
VlnPlot(data.tlineage, features = "Scin", assay = "SCT")


### integrated ----
data <- merge(data.ellen, data.tlineage)

seu <- NormalizeData(data)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:20)

seu@meta.data[is.na(seu@meta.data$celltype_singleR_ImmGenData), ]$celltype_singleR_ImmGenData <- "none"
seu$celltype <- seu$ori_celltype
seu@meta.data[seu@meta.data$celltype_singleR_ImmGenData == "Stem cells (SC.MPP34F)", ]$celltype <- "MPP"

seu$celltype <- factor(gsub("T\\.", "", seu$celltype), levels = c("HSC", "MLP", "MPP", "ETP", "DN2", "DN3", "DN4", "ISP", "DP", "NKT", "Tgd"))
VlnPlot(seu, features = "Scin", group.by = "celltype")
seu <- JoinLayers(seu)

MTX <- LayerData(seu, assay = "RNA", layer = "data")
seu$Scin_expression <- MTX["Scin", ]

# VlnPlot(seu,features = "Scin_expression",group.by = "celltype")

df <- seu@meta.data %>% dplyr::select(celltype, Scin_expression)
df <- df %>%
  group_by(celltype) %>%
  mutate(Scin = mean(Scin_expression)) %>%
  dplyr::select(celltype, Scin) %>%
  unique.data.frame()
df$celltype <- factor(df$celltype, levels = levels(seu@meta.data$celltype))

# plot
p <- ggplot(seu@meta.data, aes(x = celltype, y = Scin_expression)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1, color = "grey20") +
  geom_point(data = df, aes(x = celltype, y = Scin), color = "salmon") +
  geom_line(data = df, aes(x = celltype, y = Scin), group = 1, color = "salmon") +
  theme_classic()

ggsave(filename = "out/Fig.S7A Scin.expression.withEllen.pdf", height = 3, width = 6)


########################################################################################################
#### TF scenic analysis with Ellen data
########################################################################################################
setwd("~/Project/YuYong_TALL/out/SCENIC_mergeCTL_Ellen")

## 1.extract out_SCENIC.loom data
loom <- open_loom("aucell_output.loom")

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

## 2. load data
{
  seu <- readRDS("clone.inte.ccr.rds")
  DimPlot(seu)
  ctl <- readRDS("data.ctl.assigned.rds") # which is finished in utils/Normal-lineage.R
  DimPlot(ctl, group.by = "orig.ident")
  ctl_use <- subset(ctl, orig.ident %in% c("DN", "THY"))
  DimPlot(ctl_use, label = T)
  DimPlot(ctl_use, group.by = "orig.ident", label = T)
  main_cell <- table(ctl_use$ident) %>%
    sort() %>%
    data.frame() %>%
    filter(Freq > 50) %>%
    pull(Var1) %>%
    as.character()

  ctl_use_main <- subset(ctl_use, ident %in% main_cell)
  DotPlot(ctl_use_main, features = c(
    "Ptprc", "Cd3d", "Cd3g", "Cd3e", "Cd4",
    "Cd8a", "Ms4a1", "Cd68", "Cd34", "Kit", "Hes1", "Il2ra"
  ), cluster.idents = T)

  normal <- ctl_use_main
  normal <- RunUMAP(normal, dims = 1:30)
  DimPlot(normal, label = T, group.by = "orig.ident")
  normal$orig.ident <- paste0("N_", as.character(normal$orig.ident))

  seu@assays$SCT <- NULL
  DefaultAssay(seu) <- "RNA"
  normal@assays$SCT <- NULL
  DefaultAssay(normal) <- "RNA"
  seu@reductions
  normal@reductions

  data <- merge(seu, normal)
  data <- JoinLayers(data, assay = "RNA", layers = "counts")
  saveRDS(data, "merge_ctl&clone.rds")
} # code to get merge_ctl&clone.rds
data <- readRDS("merge_ctl&clone.rds")
{
  data <- readRDS("merge_ctl&clone.rds")

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

  saveRDS(seu, "normal.scenic.RDS")
} # code to get normal.scenic.RDS, which is to run SCENIC in normal samples(DN & Thy)
normal <- readRDS("normal.scenic.RDS")
Ellen <- readRDS("ref_data/GSE183026_ScienceImmunology2022/out/for_SCIN.rds") # # which is finished in utils/Ellen-data.R


data$celltype <- as.character(data$ident)
data@meta.data[normal@meta.data %>% rownames(), ]$celltype <- paste0("N-", normal@active.ident %>% as.character())

data$celltype %>% table()

Thy <- subset(Ellen, orig.ident == "Thy")
Thy$celltype <- paste0("Ellen-", Thy@active.ident %>% as.character())
Thy@meta.data[Thy@meta.data$orig.ident == "Thy", ]$orig.ident <- "Ellen_Thy"
Thy$celltype %>% table()

seu <- merge(data, Thy)

seu$orig.ident %>% table()
seu$celltype %>% table()

seu$celltype



## check expression --
# seurat.data<-seu
# # DimPlot(seu)
# seurat.data<-JoinLayers(seurat.data,assay = "RNA")
# seurat.data<-NormalizeData(seurat.data)
#
# p1<-DotPlot(Thy,features = c("Ptprc","Kit","Bcl11a","Il2ra","Bcl11b","Lef1","Cd3g","Cd4","Cd8a"))+
#   coord_flip()+xlab("")+ylab("")+theme_bw()+
#   theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0,size = 12,color="black"),
#         axis.text.y = element_text(size = 12,color="black")
#   )+ggtitle("Ellen-Thy")


seu@meta.data[seu@meta.data$celltype == "Ellen-ETP", ]$celltype <- "N-ETP"
seu@meta.data[seu@meta.data$celltype == "Ellen-DN2", ]$celltype <- "N-DN2"
seu@meta.data[seu@meta.data$celltype == "Ellen-DN3", ]$celltype <- "N-DN3"

seu$celltype <- factor(seu$celltype, levels = c(
  "N-ETP", "N-DN2", "N-DN3", "N-DN4", "N-DP", "N-CD8+",
  "DN2 (11)", "DN4 (12)", "DP (3)", "CD8+ (other)"
))

seu <- JoinLayers(seu, assay = "RNA")
seu <- NormalizeData(seu)



DotPlot(seu,
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

saveRDS(seu, "scenic_merge3data.rds")
seurat.data <- seu

## 3. view

sub_regulonAUC <- regulonAUC[, match(colnames(seurat.data), colnames(regulonAUC))]
# dim(sub_regulonAUC)
# seurat.data
# identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(
  row.names = colnames(seurat.data),
  seurat_clusters = as.character(seurat.data$seurat_clusters)
)
cellTypes <- data.frame(
  row.names = colnames(seurat.data),
  celltype = seurat.data$celltype
)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4, 1:4]

# save
save(sub_regulonAUC, cellTypes, cellClusters, seurat.data,
  file = "for_rss_and_visual.Rdata"
)
load("for_rss_and_visual.Rdata")

## plot ----
# TF activity mean

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


### 4.2. rss查看特异TF
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
  "Twist1(+)", "Meis1(+)", "Spi1(+)", "Gfi1b(+)", "Erg(+)", "Tcf7l2(+)", "Nfe2(+)", "Hmga2(+)", "Jun(+)", "Pbx1(+)", "Tfap4(+)",
  "Otx2(+)", "Myc(+)", "Runx1(+)",
  "Stat1(+)", "Relb(+)", "Maf(+)", "Runx3(+)",
  "Ets2(+)", "Bcl6(+)", "Klf13(+)",
  "E2f1(+)", "Tfdp1(+)", "Hdac2(+)",
  "E2f8(+)", "Vdr(+)",
  "Tal1(+)", "Runx1(+)", "Bcl11b(+)", "Myb(+)", "Myc(+)", "Dnmt1(+)", "E2f7(+)", "Tcf12(+)", "Tcf3(+)", "Ikzf1(+)", "Tcf7l1(+)", "Hoxa10(+)", "Ets1(+)",
  "Smad4(+)", "Tef(+)", "Gata2(+)",
  "Emx2(+)", "Ikzf3(+)", "Dpf1(+)", "Cebpa(+)",
  "Atf5(+)", "Tbx21(+)"
) %>% unique()
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

pheatmap::pheatmap(n, cluster_rows = F, cluster_cols = F)

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

pdf("out/Fig.2F TF regulon.withEllen.pdf", width = 5.2, height = 8)
draw(ht)
dev.off()
