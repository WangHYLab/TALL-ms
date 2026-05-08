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

saveRDS(data2,"for_SCIN.rds")