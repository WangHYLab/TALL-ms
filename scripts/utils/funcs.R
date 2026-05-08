process<-function(samplename,project = "pj"){
  data<-Read10X_h5(paste0('/home/zhengjie/Project/YuYong_TALL/cellranger/',samplename,'/outs/filtered_feature_bc_matrix.h5'))
  data<-CreateSeuratObject(data,project = project,min.cells = 5,min.features = 500)
  
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  data[["percent.rp"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")

  data <- subset(data,subset = nFeature_RNA < 6000 & percent.mt < 10)
  # data <- NormalizeData(data)
  # data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
  # data <- ScaleData(data) 
  data <- SCTransform(data)
  data <- RunPCA(data)
  data <- RunUMAP(data,dims = 1:30)
  data <- FindNeighbors(data)
  data <- FindClusters(data) 
  # p<-ElbowPlot(data)
  # ggsave(p,filename = paste0(samplename,"_elbowplot.pdf"),width = 5,height = 5)
  
  data<-removeDoublets_new(data,remove=TRUE,SCT=TRUE)
  return(data)
}

library(DoubletFinder)
removeDoublets <- function(data,remove=TRUE) {
  sweep <- paramSweep_v3(data,  sct = F)
  sweep <- summarizeSweep(sweep, GT = FALSE)
  bcmvn <- find.pK(sweep)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  # find nExp
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  DoubletRate = ncol(data)*8*1e-6 #https://www.jianshu.com/p/6770c6a05287
  nExp_poi <- round(DoubletRate * ncol(data@assays$RNA@data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  # find Doublet
  data <-doubletFinder_v3(
      data,
      PCs = 1:10,
      pN = 0.25,
      pK = mpK,
      nExp = nExp_poi,
      reuse.pANN = FALSE,
      sct = F
    )
  if(remove){
    data@meta.data[,grep('pANN',names(data@meta.data))]<-NULL
    data@meta.data$DoubletFinder<-data@meta.data[,grep('DF.',names(data@meta.data))]
    data@meta.data[,grep('DF.',names(data@meta.data))]<-NULL
    print(paste0('Origin cell : ',dim(data)[2]))
    data<-subset(data,DoubletFinder=="Singlet")
    print(paste0('Singlet cell : ',dim(data)[2]))
  }
  return(data)
}

removeDoublets_new <- function(data,remove=TRUE,strict=FALSE,SCT=TRUE) {
  library(DoubletFinder)
  sweep <- paramSweep(data, PCs = 1:30, sct = SCT)
  sweep.sum <- summarizeSweep(sweep, GT = FALSE)
  bcmvn <- find.pK(sweep.sum)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  # find nExp
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  DoubletRate = ncol(data)*8*1e-6 #https://www.jianshu.com/p/6770c6a05287
  nExp_poi <- round(DoubletRate * nrow(data@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  data <- doubletFinder(data, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = SCT)
  DF.use<-colnames(data@meta.data)[grep("DF.classifications",colnames(data@meta.data))]
  data@meta.data$DoubletFinder<-data@meta.data[,DF.use]
  data@meta.data[,DF.use]<-NULL
  
  if(strict){
    pANN.use<-colnames(data@meta.data)[grep("pANN_",colnames(data@meta.data))]
    data <- doubletFinder(data, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = pANN.use, sct = SCT)
    DF.use<-colnames(data@meta.data)[grep("DF.classifications",colnames(data@meta.data))]
    data@meta.data$DoubletFinder<-data@meta.data[,DF.use]
    data@meta.data[,grep('DF.',names(data@meta.data))]<-NULL
  }
  
  if(remove){
    data@meta.data[,grep('pANN',names(data@meta.data))]<-NULL
    print(paste0('Origin cell : ',dim(data)[2]))
    data<-subset(data,DoubletFinder=="Singlet")
    print(paste0('Singlet cell : ',dim(data)[2]))
  }
  return(data)
}

cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

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