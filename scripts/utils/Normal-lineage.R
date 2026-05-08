# this script is to extract T lineage from normal hematopoiesis and do trajectory analysis

library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
#### 数据读取、过滤----
process<-function(samplename,project = "pj",onlyread=FALSE){
  data<-Read10X_h5(paste0('cellranger/',samplename,'/outs/filtered_feature_bc_matrix.h5'))
  data<-CreateSeuratObject(data,project = project,min.cells = 5,min.features = 500)
  
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  data[["percent.rp"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")
  if(onlyread){
    return(data)
  }
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


DN  <- process(samplename ="20240126_DN_RNA", project ="DN" ,onlyread = TRUE )
THY <- process(samplename ="20240126_Thy_RNA", project ="THY",onlyread = TRUE )
BM  <- process(samplename ="20240325_BM",project = "BM",onlyread = TRUE)
LK  <- process(samplename ="20240325_LK",project = "LK",onlyread = TRUE)
LSK <- process(samplename ="20240325_LSK",project = "LSK",onlyread = TRUE)

data.raw<-merge(DN,c(THY,BM,LK,LSK))
data.raw$orig.ident<-factor(data.raw$orig.ident,levels = c("DN","THY","LSK","LK","BM"))

VlnPlot(data.raw,features = c('nCount_RNA','nFeature_RNA'),group.by = "orig.ident",pt.size = 0)
saveRDS(data.raw,"data.ctl.raw.rds")

# DN  <- process(samplename ="20240126_DN_RNA", project ="DN")
# THY <- process(samplename ="20240126_Thy_RNA", project ="THY")
data.list<-readRDS("../data.list.filtered&removeDB.rds")
THY<-data.list$Thy5
THY$orig.ident<-"THY"
DN<-data.list$DN5
DN$orig.ident<-"DN"

BM  <- process(samplename ="20240325_BM",project = "BM")
LK  <- process(samplename ="20240325_LK",project = "LK")
LSK <- process(samplename ="20240325_LSK",project = "LSK")
data.list<-list("DN"=DN,"THY"=THY,"LSK"=LSK,"LK"=LK,"BM"=BM)
saveRDS(data.list,"../data.ctl.list.filtered&removeDB.rds")

#### 数据整合 ####
data<-merge(DN,c(THY,LSK,LK,BM))  #调整顺序，用于scvelo
data$orig.ident<-factor(data$orig.ident,levels = c("DN","THY","LSK","LK","BM"))
data.bk<-data

VlnPlot(data.raw,features = c('nCount_RNA','nFeature_RNA','percent.mt',"percent.rp"),ncol = 2,group.by = "orig.ident",pt.size = 0)
ggsave(filename = "../out_ctl/data.ctl.qc.pdf",width = 6,height = 5)

data <- SCTransform(data.bk,return.only.var.genes = FALSE,variable.features.n = 3000)
data <- RunPCA(data)
{
  Seurat::ElbowPlot(data,ndims = 50)
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
  pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
  ElbowPlot(data,ndims = 50)$data %>% ggplot() +
    geom_point(aes(x = dims,y = stdev)) +
    geom_vline(xintercept = pc.use, color = "darkred") +
    theme_bw() + labs(title = "Elbow plot: quantitative approach")
  ggsave(filename = "../out_ctl/data.ctl.elbowplot.pdf",width = 5.5,height =5)
}
data <- RunUMAP(data, dims = 1:25)
data <- FindNeighbors(data, dims = 1:25)
data <- FindClusters(data, resolution = 0.5)




saveRDS(data,"../data.ctl.rds")


#### 画图   ####
{
  marker<- c("Ptprc","Cd3d","Cd3g","Cd3e",
           "Klrb1c","Itgam","Itgax","Ly76","Ly6a",
           "Cd4","Cd8a","Ccr9","Ccr7",
           "Kit","Cd34","Cd44","Cd2","Cd5","Cd7","Il2ra","Cd24a",
           "Foxp3","Gzmb","Lef1","Sell",
           "Nkg7","Ncam1",'Gnly',
           'Fcgr4',
           "Cd79a","Ms4a1","Cd19",
           "Notch1",
           "Retnlg","S100a8","S100a9",
           "Cd86","Mrc1","Cd68"
           )
}

p1<-DimPlot(data, reduction = "umap",group.by = "seurat_clusters",label = T)+theme(plot.margin = margin(1,0,0,0,"cm"))  
p2<-DimPlot(data, reduction = "umap",group.by = "orig.ident")+theme(plot.margin = margin(1,0,0,0,"cm"))  
p3<-DimPlot(data, reduction = "umap",group.by = "seurat_clusters",label = T,split.by = "orig.ident")+theme(plot.margin = margin(0,2,0,0.5,"cm"))  
p4<-DotPlot(data,features =marker,cols = c("skyblue","salmon"),
            # assay = "RNA",
            # scale = F,
            cluster.idents = T)+theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust=0.5))
ggsave(plot = cowplot::plot_grid(cowplot::plot_grid(p1,p2,ncol = 2),p3,p4,ncol = 1),filename = "../out_ctl/data.ctl.SCT-vf3k.dim30.pdf",width = 11,height =16)

p3<-DimPlot(data, reduction = "umap",group.by = "seurat_clusters",label = T,split.by = "orig.ident")+theme(plot.margin = margin(0,2,0,0.5,"cm"))  
ggsave(filename = "../out_ctl/data.ctl.SCT-vf3k.dim30.cluster-sample.pdf",width = 20,height =5)


#### singleR 注释  ####
{
  library(celldex)
  library(SingleR)

  ref<-celldex::ImmGenData()
  mtx<-LayerData(data,assay = "SCT")
  sr<-SingleR(mtx,ref =ref ,clusters = data$seurat_clusters,labels = ref$label.fine)
  sr$labels
  data$celltype_singleR_ImmGenData<-""
  data$celltype_singleR_ImmGenData_withCluster<-""
  for (cl in rownames(sr)) {
    print(cl)
    data@meta.data[data@meta.data$seurat_clusters==cl,]$celltype_singleR_ImmGenData<-sr[cl,"labels"]
    data@meta.data[data@meta.data$seurat_clusters==cl,]$celltype_singleR_ImmGenData_withCluster<-paste0("C",cl,"_",sr[cl,"labels"])
  }
}

p1<-DimPlot(data, reduction = "umap",group.by = "celltype_singleR_ImmGenData",label = T,label.size = 3,repel = T)+
  theme(plot.margin = margin(1,0,0,0,"cm")) +ggtitle("celltype_singleR_ImmGenData")+BoldTitle()

## 比例
df<-table(data$orig.ident,data$celltype_singleR_ImmGenData)%>%as.data.frame()
colnames(df)<-c("Group","Celltype","Num")
p2<-ggplot(df,aes(x=Group,y=Num,fill=Celltype))+
  geom_bar(stat = "identity",width = 0.8,position = "fill",color="white",linewidth = 0.6)+
  theme_classic()+
  ylab("Ratio")+xlab("Sample")+ggtitle("celltype-sample ratio")+BoldTitle()+
  theme(axis.text.x = element_text(angle = -45,hjust = 0,vjust=0.5),
        plot.margin = margin(1,5,0,4,"cm"),
        axis.text.x.bottom = element_text(size = 10,face = "bold"))

ggsave(plot = plot_grid(p1,p2),filename = "../out_ctl/data.ctl.SCT-vf3k.dim30.celltype_singleR_ImmGenData&ratio.pdf",width = 19,height =8)
#### marker and plot by hand  ####
{
  allgene<-data@assays$SCT@data%>%rownames()
  tcr.gd<-c(allgene[grep('Tcr[dg]',allgene)]%>%sort(),allgene[grep('Tr[dg]v',allgene)]%>%sort())
  tcr.ab<-c(allgene[grep('Tcr[ab]',allgene)]%>%sort(),allgene[grep('Tr[ab]v',allgene)]%>%sort())
  tgd<-c("Cd7","Trdv5","Trdv2-2")
  immu<-c("Ptprc","Cd53")
  lsk<-c("Atxn1","Ly6a")
  t<-c("Cd3d","Tcf7"
       # "Cd3g","Cd3e",
       )
  hsc<-c("Runx1","Mllt3","Hoxa9","Mecom","Hlf",
         # "Spink2",
         
         # "Epor",
         "Mpl","Kit","Il3ra"
         )
  cd48<-c("Cd4","Cd8a","Cd8b1")
  dn234<-c(
    # "Bcl11a","Notch1","Runx3",
    "Bcl11b","Gata3","Hes1","Kit","Il2ra")
  nkt<-c("Klrb1c","Ccl5")
  isp<-c("Cd24a","Tcrb")
  
  neu<-c("Dstn","Wfdc21","Retnlg","Csf3r","Cxcr2","Fpr2","Il1r2")
  ifim<-allgene[grep("If",allgene)]
  
  mac<-c("C1qa","C1qc","Mrc1")
  mono<-c('F13a1','Csf1r','Ms4a6c','Ly86','Slpi','Ass1')
  dc<-c("Cd209a",'Siglech',
        'Cox6a2',"Bcr", #浆细胞样
        'Irf8')
  
  tn<-c("Sell","Lef1")
  # c9<-c('Cd28','Cd5','Satb1','Cd2','Tcf7')
  c26<-c("Trbv16","Hist1h2ak")
  
  cyc<-c("Mki67")
  MEP<-c("Snca","Alas2","Gypa","Hba-a2","Gata1") # Megakaryocyte-erythroid progenitor
  Ery<-c("Car1","Cited4","Klf1","Phf10","Gata1","Hba-a2") # cellmarker2 来源红细胞marker
  MPP<-c('Cd34', # cellmarker2 来源MPP marker
         "Ly6a",
         'Kit',
         "Fcgr3","Fcgr4",#'Cd16',
         'Fcgr2b', #'Cd32',
         'Flt3', #CD135
         'Gcnt2',
         'Hlf',
         'Il7r')
  MPP<-c('Cd34', "Ly6a",'Kit','Flt3','Gcnt2','Hlf',
         "Gata2")
  GMP<-c("Fcgr3",'Fcgr2b')
  CLP<-c("Mme", #CD10
         'Cd38')
  CMP<-c("Thy1", #CD90
         # 'Cd38',
         "Mme",
         "Cd34"
         )
  MDP<-c('Csf1r', #样本来源ImmGen :Sca1- Flt3+ MCSFR+ c-Kithi
    'Mpo',
        'Ctsg',
        'Elane',
        'Ms4a3',"Dmkn",
        'Ly6c2')
  
  b<-c("Cd19","Cd79a","Ms4a1")
  b_frf<-c('Cd74','Iglc2') #"Cr2",
  plasma<-c("Sdc1")
  
  {
    data@active.ident<-data$celltype_singleR_ImmGenData_withCluster%>%factor()
    data<-RenameIdents(data,
                       'C6_Stem cells (SC.MPP34F)'="C6: HSC (MPP)",
                       'C2_Stem cells (SC.MPP34F)'="C2: HSC (MPP)",
                       'C4_Stem cells (SC.ST34F)'="C4: HSC (MPP)",
                       'C5_Stem cells (SC.ST34F)'="C5: HSC (MPP.Flt3-)",
                       
                     'C7_Stem cells (SC.CMP.DR)'="C7: HSC (CMP)",
                     'C19_Stem cells (GMP)'="C19: HSC (GMP)",
                     'C3_Stem cells (SC.MDP)'="C3: HSC (MDP)",
                     'C24_Stem cells (SC.MEP)'="C24: HSC (MEP)", #很有可能是红细胞，因为Kit没怎么表达

                     
                     'C20_T cells (T.DN3A)'="C20: T (DN3A)",
                     'C28_T cells (T.DN3A)'="C28: T (DN3A)",
                     'C1_T cells (T.DN3A)'="C1: T (DN3A)",
                     'C26_T cells (T.DN3B)'="C26: T (DN3B)",
                     
                     'C23_T cells (T.DP)'="C23: T (DP)",
                     'C0_T cells (T.DP)'="C0: T (DP)",
                     'C9_T cells (T.DP69+)'="C9: T (DP)",
                     'C22_T cells (T.ISP)'="C22: T (DP)",
                     
                     'C11_T cells (T.ISP)'="C11: T (ISP)", #T.ISP.Th#1, Immature Single-Positive
                     'C10_T cells (T.ISP)'="C10: T (ISP)",
                     
                     'C14_Tgd (Tgd.imm.VG1+)'="C14: T (Tgd)",
                     'C17_NKT (NKT.4-)'="C17: T (NKT)",
                     
                     'C13_Monocytes (MO.6C+II-)'="C13: Mo",
                     'C21_Monocytes (MO.6C+II-)'="C21: Mø",
                     'C25_DC (DC.PDC.8+)'="C25: DC (PDC)",
                     
                     'C18_B cells (proB.CLP)'="C18: Neu (pro)",
                     'C12_Neutrophils (GN.ARTH)'="C12: Neu (ARTH)", #关节炎arthritic
                     'C8_Neutrophils (GN)'="C8: Neu",

                     
                     'C16_B cells (preB.FrD)'="C16: B (preB)",
                     'C15_B cells (B.FrF)'="C15: B (FrF)",     #"快速活化和增殖"（Fast Activating and Proliferating）B细胞
                     'C27_B cells (B.FrF)'="C27: B (Plasma)"
                     )
  }
  # gene<-deg.hsc%>%filter(p_val_adj<0.05&avg_log2FC>0.5&pct.1>0.5)%>%group_by(cluster)%>%slice_head(n = 40)
  
  DotPlot(data,
          features =c(immu,hsc,MPP,GMP,MDP,MEP,CMP,t,dn234,c26,cd48,isp,tn,tgd,nkt,mono,mac,dc,neu,b,b_frf,plasma)%>%unique(),
          # features=c("Bcl11a","Bcl11b",MPP,GMP,MDP,MEP,CLP,CMP)%>%unique(),
          cols = c("lightblue","red"),
          # col.max = 1,col.min = 0,
          # scale = F,
          # group.by = "ident",
          cluster.idents = F)+
    # coord_flip()+
    xlab("")+ylab("")+
    # theme(axis.text.x = element_text(angle = -45,hjust = 0,vjust=0.5))+
    theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust=0.5))
  ggsave(filename = "../out_ctl/data.ctl.SCT-vf3k.dim30.celltype.marker.dotplot.pdf",width = 15,height = 6)
  
  data$ident<-data@active.ident
  p1<-DimPlot(data,group.by = "ident", reduction = "umap",label = T,label.size = 3,repel = T,pt.size = 0.8)+
    theme(plot.margin = margin(1,0,0,0,"cm")) +ggtitle("celltype-cluster")+BoldTitle()
  
  ## 比例
  df<-table(data$orig.ident,data$ident)%>%as.data.frame()
  colnames(df)<-c("Group","Celltype","Num")
  p2<-ggplot(df,aes(x=Group,y=Num,fill=Celltype))+
    geom_bar(stat = "identity",width = 0.8,position = "fill",color="white",linewidth = 0.6)+
    theme_classic()+
    ylab("Ratio")+xlab("Sample")+ggtitle("celltype-sample ratio")+BoldTitle()+
    theme(axis.text.x = element_text(angle = -45,hjust = 0,vjust=0.5),
          plot.margin = margin(1,5,0,4,"cm"),
          axis.text.x.bottom = element_text(size = 10,face = "bold"))
  ggsave(plot = p2,filename = "../out_ctl/data.ctl.SCT-vf3k.dim30.ratio.pdf",width = 10,height =8)
  
  ggsave(plot = plot_grid(p1,p2,ncol=1),filename = "../out_ctl/data.ctl.SCT-vf3k.dim30.celltype&ratio.pdf",width = 11,height =16)
  
}

saveRDS(data,"data.ctl.assigned.rds")


#### change cluster

library(slingshot)
library(RColorBrewer)
library(paletteer) 
library(viridis)
# BiocManager::install("tradeSeq")
library(tradeSeq)
library(tidyverse)
library(Seurat)
data.bk<-readRDS("./data.ctl.assigned.rds")
ccgene<-cc.genes.updated.2019
library(homologene)
s<-homologene(ccgene$s.genes,inTax = "9606",outTax = "10090")$`10090`
g2m<-homologene(ccgene$g2m.genes,inTax = "9606",outTax = "10090")$`10090`


DimPlot(data.bk,label = T,label.size = 3)

data<-subset(data.bk,ident%in%c(
  "C6: HSC (MPP)","C2: HSC (MPP)","C4: HSC (MPP)","C5: HSC (MPP.Flt3-)",
  "C20: T (DN3A)","C28: T (DN3A)","C1: T (DN3A)","C26: T (DN3B)",
  "C23: T (DP)","C0: T (DP)","C9: T (DP)","C22: T (DP)",
  "C11: T (ISP)", "C10: T (ISP)",
  "C14: T (Tgd)","C17: T (NKT)"))

gene<-data@assays$SCT$counts%>%rownames()

data<-RunUMAP(data,dims = 1:20,reduction.name = "umap.re")

## 重新SCT
data2 <- SCTransform(data,return.only.var.genes = FALSE,variable.features.n = 3000)
data3 <-data2
marker<-c(  'Cd34', "Ly6a",'Kit','Flt3','Gcnt2','Hlf',"Gata2",
            "Ptprc","Cd53",
            "Runx1","Mllt3","Hoxa9","Mecom","Hlf","Mpl","Kit","Il3ra",
            "Cd3d","Tcf7","Cd3g","Cd3e",
            'Notch1',
            "Bcl11b","Gata3","Hes1","Kit","Il2ra",
            "Cd4","Cd8a","Cd8b1",
            "Cd24a",
            "Klrb1c","Ccl5",
            "Cd7","Trdv5","Trdv2-2",
            "Sell","Lef1")
VariableFeatures(data3)<-c(setdiff(VariableFeatures(data.bk),c(s,g2m)),marker)%>%unique()
data <- RunPCA(data3)
{
  Seurat::ElbowPlot(data,ndims = 50)
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
  pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
  ElbowPlot(data,ndims = 50)$data %>% ggplot() +
    geom_point(aes(x = dims,y = stdev)) +
    geom_vline(xintercept = pc.use, color = "darkred") +
    theme_bw() + labs(title = "Elbow plot: quantitative approach")
  # ggsave(filename = "../out_ctl/data.ctl.elbowplot.pdf",width = 5.5,height =5)
}

data <- RunUMAP(data, dims = 1:11,n.neighbors = 30,min.dist = 0.3) #默认
data <- RunUMAP(data, dims = 1:11,n.neighbors = 30,min.dist = 0.3,reduction.name = "umap.dist0.4") #默认
DimPlot(data,reduction = "umap",group.by = "ident",label = T,label.size = 4)

data <- FindNeighbors(data, dims = 1:11)
data <- FindClusters(data, resolution = 0.5)
DimPlot(data,reduction = "umap",group.by = c("ident","seurat_clusters","orig.ident"),label = T,label.size = 4)

data<-FindSubCluster(data,graph.name = "SCT_snn",cluster = 1,resolution = 0.5,subcluster.name = "DN_subcluster")
DimPlot(data,reduction = "umap",group.by = c("ident","DN_subcluster","orig.ident"),label = T,label.size = 4)


table(data$DN_subcluster,data$ident)%>%heatmap()
DotPlot(data,
        features =marker%>%unique(),
        cols = c("lightblue","red"),
        group.by = "celltype_singleR_ImmGenData_withCluster",
        cluster.idents = F)+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust=0.5))

# singleR annotaiton
{
  library(celldex)
  library(SingleR)
  
  ref<-celldex::ImmGenData()
  mtx<-LayerData(data,assay = "SCT")
  sr<-SingleR(mtx,ref =ref ,clusters = data$DN_subcluster,labels = ref$label.fine)
  sr$labels
  data$celltype_singleR_ImmGenData<-""
  data$celltype_singleR_ImmGenData_withCluster<-""
  for (cl in rownames(sr)) {
    print(cl)
    data@meta.data[data@meta.data$DN_subcluster==cl,]$celltype_singleR_ImmGenData<-sr[cl,"labels"]
    data@meta.data[data@meta.data$DN_subcluster==cl,]$celltype_singleR_ImmGenData_withCluster<-paste0("C",cl,"_",sr[cl,"labels"])
  }
  }



saveRDS(data,"./data.Tlineage.use.rds")
data@active.ident<-data$celltype_singleR_ImmGenData_withCluster%>%factor()
data<-RenameIdents(data,
                   "C4_Stem cells (MLP)"      ="HSC",         
                   "C3_Stem cells (SC.ST34F)"    ="HSC",        
                   "C2_Stem cells (SC.ST34F)"    ="HSC",       
                   "C5_Stem cells (SC.ST34F)"  ="HSC",
                   "C17_Stem cells (SC.MPP34F)"    ="ETP", 
                   "C1_4_T cells (T.DN3A)"   ="T.DN2",
                   "C1_5_T cells (T.DN3A)"   ="T.DN3",             
                   "C1_2_T cells (T.DN3A)"   ="T.DN3",            
                   "C1_6_T cells (T.DN3A)"   ="T.DN3",
                   "C1_1_T cells (T.DN3A)"   ="T.DN3",
                   "C1_0_T cells (T.DN3A)"   ="T.DN3",             
                   "C1_3_T cells (T.DN3A)"   ="T.DN3",    
                   "C1_7_T cells (T.DN3A)"   ="T.DN3",             
                   "C13_T cells (T.DN3A)"   ="T.DN3", 
                   "C6_T cells (T.DN3B)"     ="T.DN4",               
                   "C7_T cells (T.ISP)"      ="T.DN4",                
                   "C15_T cells (T.ISP)"     ="T.ISP",                 
                   "C11_T cells (T.DP)"         ="T.DP",         
                   "C12_T cells (T.DPbl)"       ="T.DP",        
                   "C16_T cells (T.DPbl)"        ="T.DP",        
                   "C8_T cells (T.DP69+)"        ="T.DP",       
                   "C0_T cells (T.DP)"           ="T.DP", 
                   "C9_Tgd (Tgd.imm.VG1+)"   ="Tgd", 
                   "C14_Tgd (Tgd.imm.vg2+)"  ="Tgd", 
                   "C10_T cells (T.8MEM.OT1.D45.LISOVA)" ="NKT"
                   
                   )
data$Tlineage<-data@active.ident
p<-DimPlot(data,group.by = c("orig.ident","DN_subcluster","Tlineage","ident"),ncol = 2,label = T,label.size = 4)

DotPlot(data,
        features = c("Ptprc","Cd3d","Cd34","Hlf","Gcnt2","Kit",'Notch1',"Il2ra","Hes1","Mki67","Hist1h2ak","Cd8a","Cd4","Cd24a","Trdv5","Trdv2-2","Cpa3","Klrb1c","Ccl5","Klra1"
                     ),
        # features =marker%>%unique(),
        # features =c("Cd3d","Cd4","Cd8a","Kit","Il2ra","Cd44","Cd24a","Gzmb","Foxp3","Tcrb","Sell","Lef1","Ncr1"),
        # features = c(gene[grep("Trg",gene)],gene[grep("Tcrg",gene)],"Cd7","Trdv5","Trdv2-2"),
        # features = c(gene[grep("Trb",gene)],gene[grep("Tcrb",gene)]),
        cols = c("skyblue","red"),
        cluster.idents = F)+
  xlab("")+ylab("")+coord_flip()+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust=0.5))

data<-PrepSCTFindMarkers(data)
deg<-FindAllMarkers(data,only.pos = T)
deg.sig<-deg%>%filter(p_val_adj<0.05&pct.1>0.5&avg_log2FC>1)

saveRDS(data,"./data.Tlineage.assign.rds")




