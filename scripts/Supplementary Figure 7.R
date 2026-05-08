#########################################################################################
#### Fig.S7A
#########################################################################################

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(homologene)

data.ellen<-readRDS("ref_data/GSE183026_ScienceImmunology2022/out/for_SCIN.rds") # which is finished in utils/Ellen-data.R
data.ellen$orig.ident<-paste0("Ellen.",data.ellen$orig.ident)
data.ellen$ori_celltype<-data.ellen@active.ident


data.tlineage<-readRDS("../data.Tlineage.assign.rds")  # which is finished in scripts/utils/Normal-lineage.R
data.tlineage$ori_celltype<-data.tlineage@active.ident


data.tlineage$Scin_scaled<-data.tlineage@assays$SCT@scale.data["Scin",]
ggplot(data.tlineage@meta.data,aes(x=ori_celltype,y=Scin_scaled))+
  geom_boxplot()

VlnPlot(data.ellen,features = "Scin",assay = "RNA")
VlnPlot(data.tlineage,features = "Scin",assay = "SCT")


### integrated ----
data<-merge(data.ellen,data.tlineage)

seu <- NormalizeData(data)
seu <- FindVariableFeatures(seu, selection.method = "vst",nfeatures=3000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims=1:20)

seu@meta.data[is.na(seu@meta.data$celltype_singleR_ImmGenData),]$celltype_singleR_ImmGenData<-"none"
seu$celltype<-seu$ori_celltype
seu@meta.data[seu@meta.data$celltype_singleR_ImmGenData=="Stem cells (SC.MPP34F)",]$celltype<-"MPP"

seu$celltype<-factor(gsub("T\\.", "", seu$celltype),levels = c("HSC",'MLP','MPP',"ETP","DN2","DN3","DN4","ISP","DP","NKT","Tgd"))
VlnPlot(seu,features = "Scin",group.by = "celltype")
seu <- JoinLayers(seu)

MTX<-LayerData(seu,assay = "RNA",layer = "data")
seu$Scin_expression<-MTX["Scin",]

# VlnPlot(seu,features = "Scin_expression",group.by = "celltype")

df<-seu@meta.data %>% dplyr::select(celltype,Scin_expression)
df<-df %>% group_by(celltype) %>% mutate(Scin=mean(Scin_expression)) %>% dplyr::select(celltype,Scin) %>% unique.data.frame()
df$celltype <- factor(df$celltype, levels = levels(seu@meta.data$celltype))

# plot
p<-ggplot(seu@meta.data, aes(x=celltype, y=Scin_expression)) +
  stat_boxplot(geom="errorbar", width=0.2) +            
  geom_boxplot(outlier.alpha=0) +                       
  geom_jitter(width = 0.2, alpha=0.2, size=1, color="grey20") +
  geom_point(data=df, aes(x=celltype, y=Scin), color="salmon") +   
  geom_line(data=df, aes(x=celltype, y=Scin), group=1,color="salmon") +       
  theme_classic()

ggsave(filename = "out/Fig.S7A Scin.expression.ellen-tlineage.withMean.pdf",height = 3,width = 6)