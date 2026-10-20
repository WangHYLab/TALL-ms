#########################################################################################
#### Fig.S7C
#########################################################################################
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
# library(fanyi)
library(corrplot)
library(openxlsx)
library(ggrepel)
library(enrichplot)
library(cowplot)
library(eulerr)

#### 差异表达 ####
run_DESeq2<-function(key1,key2,out_postfix,tpm,raw_count,write=TRUE){
  # id2gene<-read.table("/home/zhengjie/DATA/ensmusg2symbol.GRCm39.txt",header = T,sep = "\t",row.names = 1)
  id2gene<-read.table("data/mm10.id2symbol.txt",header = T,sep = "\t",row.names = 1)
  id2description<-read.table("data/ensmusg2description.GRCm39.txt",header = T,sep = "\t",row.names = 1,check.names = F,quote = "")
  sample.re<-read.xlsx("data/sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])
  
  # sample<-colnames(raw_count)
  # dis<-sample[grep(key1,sample)]
  # ctl<-sample[grep(key2,sample)]
  dis<-key1
  ctl<-key2
  data <- raw_count[,c(dis,ctl)]
  expr<-tpm[,c(dis,ctl)]
  
  condition<-c(rep("dis",length(dis)),rep("ctl",length(ctl)))
  coldata<-data.frame(row.names = colnames(data), condition)
  print(coldata)
  dds<-DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)
  # dds$condition <- relevel(dds$condition,"ctl")
  dds<-DESeq(dds)
  res<-results(dds)
  out<-res[order(res$padj),] %>%data.frame()
  out$gene<-id2gene[rownames(out),]
  out$description<-id2description[rownames(out),]
  
  out<-out%>%dplyr::select("gene","description","baseMean","log2FoldChange", "lfcSE","stat" ,"pvalue","padj")
  out<-cbind(out,expr[rownames(out),])
  if(write){
    write.xlsx(out,file = paste0("data/",out_postfix,".DESeq2.withTPM.xlsx"),overwrite = TRUE)
  }
  return(out)
}
run_gsea<-function(deseq2_res,file_postfix="",write_xlsx=TRUE,merge_save=FALSE){
  library(GSEABase) 
  library(clusterProfiler)
  library(DOSE)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(openxlsx)
  
  geneList<-deseq2_res$stat
  names(geneList)=deseq2_res$gene
  geneList=sort(geneList,decreasing = T)
  
  message("Run GSEA in MouseHallmark...")
  gmt.mh<-read.gmt('data/mh.all.v2023.2.Mm.symbols.gmt')
  egmt.mh <- GSEA(geneList, TERM2GENE=gmt.mh,eps = 0,seed = 1,nPermSimple=10000,pvalueCutoff = 1)
  out.mh<-egmt.mh@result
  message(paste0("Sig path: ",sum(out.mh$pvalue<0.05%>%as.numeric())))
  
  message("Run GSEA in Reactome...")
  gmt.reactome<-read.gmt('data/m2.cp.reactome.v2023.2.Mm.symbols.gmt')
  egmt.reactome <- GSEA(geneList, TERM2GENE=gmt.reactome,eps = 0,seed = 1,nPermSimple=10000,pvalueCutoff = 1)
  out.reactome<-egmt.reactome@result
  message(paste0("Sig path: ",sum(out.reactome$pvalue<0.05%>%as.numeric())))
  
  message("Run GSEA in GO...")
  gmt.go<-read.gmt('data/m5.go.v2023.2.Mm.symbols.gmt')
  egmt.go <- GSEA(geneList, TERM2GENE=gmt.go,eps = 0,seed = 1,nPermSimple=10000,pvalueCutoff = 1)
  out.go<-egmt.go@result
  message(paste0("Sig path: ",sum(out.go$pvalue<0.05%>%as.numeric())))
  
  if (write_xlsx) {
    write.xlsx(x = out.mh,file = paste0(file_postfix,".gsea.MouseHallmark.xlsx"),overwrite = T)
    write.xlsx(x = out.reactome,file = paste0(file_postfix,".gsea.Reactome.xlsx"),overwrite = T)
    write.xlsx(x = out.go,file = paste0(file_postfix,".gsea.GO.xlsx"),overwrite = T)
  }
  if(merge_save){
    write.xlsx(x = list(
      hallmark=out.mh,
      reactome=out.reactome,
      go=out.go
    ),file = paste0(file_postfix,".gsea.Hallmark&Reactome&GO.xlsx"),overwrite = T)
    
  }
  
  return(list(mh=egmt.mh,go=egmt.go,reactome=egmt.reactome))
}


### 方法来源于RNA-seq.mm10.R

raw_count<-read.table("../gene-count-matrix.txt",check.names = F)
colname.new<-lapply(colnames(raw_count),function(x){strsplit(strsplit(x,"/home/zhengjie/Project/YuYong_TALL/VCTB-VCAB-DN2-analysis/mapping/")[[1]][2],".genes.results")[[1]][1]})%>%as.character()
colnames(raw_count)<-colname.new
# colnames(raw_count)<-sample.re[colnames(raw_count),]
#去除未表达gene
# raw_count<-raw_count[rowSums(raw_count)!=0,]
raw_count<-round(as.matrix(raw_count))

data<-read.table("../gene-TPM-matrix.txt",check.names = F)
colname.new<-lapply(colnames(data),function(x){strsplit(strsplit(x,"/home/zhengjie/Project/YuYong_TALL/VCTB-VCAB-DN2-analysis/mapping/")[[1]][2],".genes.results")[[1]][1]})%>%as.character()
colnames(data)<-colname.new
# colnames(data)<-sample.re[colnames(data),]
tpm<-data


sample<-colnames(raw_count)


# ~/Project/YuYong_TALL/ov-DN2-analysis
raw_count2<-read.table("ov-DN2-analysis/gene-count-matrix.txt",check.names = F)
colname.new<-lapply(colnames(raw_count2),function(x){strsplit(strsplit(x,"/home/zhengjie/Project/YuYong_TALL/ov-DN2-analysis/mapping/")[[1]][2],".genes.results")[[1]][1]})%>%as.character()
colnames(raw_count2)<-colname.new
# colnames(raw_count2)<-sample.re[colnames(raw_count2),]
#去除未表达gene
# raw_count2<-raw_count2[rowSums(raw_count2)!=0,]
raw_count2<-round(as.matrix(raw_count2))

data2<-read.table("ov-DN2-analysis/gene-TPM-matrix.txt",check.names = F)
colname.new<-lapply(colnames(data2),function(x){strsplit(strsplit(x,"/home/zhengjie/Project/YuYong_TALL/ov-DN2-analysis/mapping/")[[1]][2],".genes.results")[[1]][1]})%>%as.character()
colnames(data2)<-colname.new
# colnames(data2)<-sample.re[colnames(data2),]
tpm2<-data2


sample2<-colnames(raw_count2)


tpm<-cbind(tpm,tpm2)
######## PCA  
dis<-c("VCAB-DKO-1", "VCAB-DKO-2", "VCAB-DKO-3")
ctl<-c("VCAB1","VCAB2","VCAB3","DN2-1" , "DN2-2" , "DN2-3")

expr<-tpm[,c(dis,ctl)]

data<-expr
loc<-rowSums(data)!=0
data.filter<-data[loc,]

pca<-prcomp(t(data.filter),scale. = T)
head(pca$x)
summary(pca)
pca_result<-as.data.frame(pca$x[,1:3])
sample<-colnames(data.filter)
# shape<-sample_info$condition
group<-c("DKO", "DKO", "DKO", "VCAB", "VCAB" , "VCAB", "WT", "WT", "WT")
pca_result<-cbind(pca_result,sample,group)

# 2d pca
colnames(pca_result)<-c("PC1" ,   "PC2",    "PC3"  ,  "sample", "group")
# write.csv(pca_result,"DKO&VCAB&WT.PCA.stat.csv")
library(ggplot2)
library(ggrepel)
p1<-ggplot(pca_result,aes(x=PC1,y=PC2,color=group))+
  geom_point(size=5)+
  
  geom_text_repel(aes(label=sample,x=pca_result[,1],y=pca_result[,2]),data=pca_result,size=3)+
  theme_classic()+
  theme(
    plot.margin = unit(c(1,1,1,1),'lines'))+
  labs(x="PCA 1 (37.04%)",y="PCA 2 (23.73%)")+
  scale_color_manual(values = c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69"))
p1
ggsave("out/Fig.S7C DKO&VCAB&WT.PCA.pdf", width = 4.9, height = 3.8)


#########################################################################################
#### Fig.S7D
#########################################################################################

### VCAB vs WT DN2
dis<-c("VCAB1","VCAB2","VCAB3")
ctl<-c("DN2-1" , "DN2-2" , "DN2-3")
res.DKO_VCAB<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "VCAB vs. DN2-normal.all",tpm = tpm,raw_count = raw_count,write = TRUE)
# top100差异基因热图
deg<-res.DKO_VCAB%>%filter(padj<0.05)%>%filter(abs(log2FoldChange)>1)
deg.DKO_VCAB<-deg
write.xlsx(deg.DKO_VCAB,file = "out/VCAB vs. DN2-normal.all.DEG.padj0.05&lfc1.xlsx",quote = FALSE)

# pdf(file = "VCAB vs. DN2-normal.all.DEG.padj0.05&lfc.top50.pdf",height = 10,width = 5)
# deg[1:50,c("gene",dis,ctl)] %>% distinct(gene,.keep_all = T) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% pheatmap::pheatmap(.,scale = "row",cutree_rows = 2,cutree_cols = 2)
# dev.off()

VCAB_DN2.GSEA<-run_gsea(deseq2_res = res.DKO_VCAB,file_postfix = "VCAB vs. DN2-normal",write_xlsx = TRUE,merge_save = TRUE)
saveRDS(VCAB_DN2.GSEA,"out/VCAB-DN2.GSEA.rds")




### DKO vs WT DN2
dis<-c("VCAB-DKO-1", "VCAB-DKO-2", "VCAB-DKO-3")
ctl<-c("DN2-1" , "DN2-2" , "DN2-3")
res.DKO_VCAB<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "VCAB-DKO vs. DN2-normal.all",tpm = tpm,raw_count = raw_count,write = TRUE)
# top100差异基因热图
deg<-res.DKO_VCAB%>%filter(padj<0.05)%>%filter(abs(log2FoldChange)>1)
deg.DKO_VCAB<-deg
write.xlsx(deg.DKO_VCAB, file = "out/VCAB-DKO vs. DN2-normal.all.DEG.padj0.05&lfc1.xlsx", quote = FALSE)

# pdf(file = "VCAB-DKO vs. DN2-normal.all.DEG.padj0.05&lfc.top50.pdf",height = 10,width = 5)
# deg[1:50,c("gene",dis,ctl)] %>% distinct(gene,.keep_all = T) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("gene") %>% pheatmap::pheatmap(.,scale = "row",cutree_rows = 2,cutree_cols = 2)
# dev.off()



DKO_DN2.GSEA<-run_gsea(deseq2_res = res.DKO_VCAB,file_postfix = "VCAB-DKO vs. DN2-normal",write_xlsx = TRUE,merge_save = TRUE)
saveRDS(DKO_DN2.GSEA, "out/VCAB-DN2.GSEA.rds")

library(enrichplot)
plot_gsea<-function(res,path_id,prefix){
  
  title<-paste0(path_id,"\n",
                "NES = ",res@result %>% filter(ID==path_id)%>% pull(NES) %>% signif(.,2),"\t",
                "p-value = ",res@result %>% filter(ID==path_id)%>% pull(pvalue) %>% signif(.,2),"\t",
                "adjusted p-value = ",res@result %>% filter(ID==path_id)%>% pull(p.adjust) %>% signif(.,2)
  )
  p<-gseaplot(res,
           geneSetID = path_id,
           by = "runningScore",color = "darkred",
           title = title)
  ggsave(plot=p,filename = paste0(prefix,path_id,".pdf"),width = 5.5,height = 3.2)
  return(p)
}

path_id<-"GOBP_T_CELL_DIFFERENTIATION"
res<-VCAB_DN2.GSEA$go

path<-c("GOBP_T_CELL_DIFFERENTIATION", # 1
        "GOBP_DNA_REPLICATION", # 3
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT")


p11<-plot_gsea(VCAB_DN2.GSEA$go,path_id ="GOBP_T_CELL_DIFFERENTIATION" ,prefix = "VCAB vs WT")
p12<-plot_gsea(DKO_DN2.GSEA$go,path_id ="GOBP_T_CELL_DIFFERENTIATION" ,prefix = "DKO vs WT")

p21<-plot_gsea(VCAB_DN2.GSEA$go,path_id ="GOBP_DNA_REPLICATION" ,prefix = "VCAB vs WT")
p22<-plot_gsea(DKO_DN2.GSEA$go,path_id ="GOBP_DNA_REPLICATION" ,prefix = "DKO vs WT")

p31<-plot_gsea(VCAB_DN2.GSEA$mh,path_id ="HALLMARK_E2F_TARGETS" ,prefix = "VCAB vs WT")
p32<-plot_gsea(DKO_DN2.GSEA$mh,path_id ="HALLMARK_E2F_TARGETS" ,prefix = "DKO vs WT")

p41<-plot_gsea(VCAB_DN2.GSEA$mh,path_id ="HALLMARK_G2M_CHECKPOINT" ,prefix = "VCAB vs WT")
p42<-plot_gsea(DKO_DN2.GSEA$mh,path_id ="HALLMARK_G2M_CHECKPOINT" ,prefix = "DKO vs WT")

library(patchwork)
cowplot::plot_grid(p11,p12,p21,p22,p31,p32,p41,p42,ncol = 2)
ggsave(plot = cowplot::plot_grid(p11,p12,p21,p22,p31,p32,p41,p42,ncol = 2),filename = "out/Fig.S7D VCAB&DKO-WT.GESA.pdf",width = 9,height = 10)



#########################################################################################
#### Fig.S7E
#########################################################################################
## 绘图
deg.DKO_VCAB_DN<-openxlsx::read.xlsx("VCAB-DKO vs. VCAB+DN2-normal.all.DESeq2.withTPM.xlsx")
gene_for_plot<-c('Spi1','Hoxa9','Hoxa7','Erg','Scin','Lmo2','Mef2c','Cd34','Mycn','Cd2','Zap70','Cd3d','Ptcra','Lef1','Rag1','Cd8a','Cd8b1','Carmil2'
)           
plot<-deg.DKO_VCAB_DN %>% filter(gene %in%gene_for_plot ) %>% tibble::column_to_rownames("gene") %>% dplyr::select("DN2-1" , "DN2-2" , "DN2-3",
                                                                                                                "VCAB1","VCAB2","VCAB3",
                                                                                                                "VCAB-DKO-1", "VCAB-DKO-2", "VCAB-DKO-3")

bk <- c(seq(-1.5, 1.5, by = 0.01)) %>% unique()
col_palette <- colorRampPalette(
  colors = RColorBrewer::brewer.pal(n = 5,name = "Set1")
)(length(bk) - 1)


bk_low  <- seq(-1.5, -0.8, by = 0.03)
bk_mid  <- seq(-0.8, 0.2, by = 0.005)
bk_high <- seq(0.2, 1.5, by = 0.05)
bk <- unique(c(bk_low, bk_mid, bk_high))

col_palette <- colorRampPalette(
  colors = c("blue","white","red", "#AF0003",  "#7F0000", "#3D0000")
)(length(bk) - 1)



col_anno <- data.frame(
  Group = factor(c(rep("DN2", 3), rep("VCAB", 3), rep("VCAB-DKO", 3))),
  row.names = colnames(plot[gene_for_plot, ])
)
ann_color <- list(
  Group = c("VCAB-DKO" = "#8dd3c7", "VCAB" = "#ffffb3", "DN2" = "#bebada")
)

pdf(file = "out/Fig.S7E VCAB-DKO vs. VCAB+DN2-normal.heatmap.color-manul.pdf",height = 7,width = 5)

# 4. 绘制热图
pheatmap::pheatmap(
  mat = plot[gene_for_plot, ],
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  gaps_row = 9,          # 行分隔线，和你原有设置一致
  gaps_col = c(3, 6),    # 列分隔线，三组之间加白色竖线
  color = col_palette,
  breaks = bk,
  border_color = "white",
  annotation_col = col_anno,
  annotation_colors = ann_color,
  show_rownames = TRUE,
  show_colnames = TRUE
)


dev.off()



#########################################################################################
#### Fig.S7F
#########################################################################################
setwd("cuttag-analysis")
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(corrplot)

anno_peak<-function(bed_file,name){
  
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  #过滤BED文件 
  data<-read.table(bed_file)
  data<-data%>%filter(V1%in%paste0("chr",c(1:19,"X","Y","M")))
  # data$V1<-paste0("chr",data$V1)
  # data[data$V1=="MT",]$V1<-"M"
  data$V1%>%unique()
  write.table(data,file = file.path(dirname(bed_file),paste0(strsplit(basename(bed_file),".bed")[[1]],".clean.bed")),quote = F,col.names = F,row.names = F,sep = "\t")
  
  peak=readPeakFile(peakfile = file.path(dirname(bed_file),paste0(strsplit(basename(bed_file),".bed")[[1]],".clean.bed")))
  
  
  # GenomeInfoDb::seqlevels(txdb)
  # peak <- diffloop::addchr(peak)
  print(GenomeInfoDb::seqlevels(peak))
  print(GenomeInfoDb::seqlevels(peak)%in%GenomeInfoDb::seqlevels(txdb))
  
  out = annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=5000,annoDb = "org.Mm.eg.db")
  write.xlsx(as.data.frame(out),file = paste0(dirname(bed_file),"/",name,".xlsx"))
  write.table(as.data.frame(out),file =  paste0(dirname(bed_file),"/",name,".bed"),quote = F,col.names = F,row.names = F,sep = "\t")
  df<-as.data.frame(out)
  df$end<-df$end+1
  write.table(df,file =   paste0(dirname(bed_file),"/",name,".end+1.bed"),quote = F,col.names = F,row.names = F,sep = "\t")
  return(out)
}

require(TxDb.Mmusculus.UCSC.mm10.knownGene)

## get file (Bcl11a_bound_enhancers_annotated.bed) by run scripts/utils/Bcl11a-enhancer.sh
Bcl11a<-anno_peak("cuttag-analysis/peakCalling_mm10/B_merge_macs3_summits.bed",name = "B_merge")
Bcl11a_Enhancer<-anno_peak("ref_data/GSE283392/peakCalling/Bcl11a_bound_enhancers_annotated.bed",name = "Bcl11a_Enhancer")
### enhancer promoter的比例
Bcl11a<-read.xlsx("cuttag-analysis/peakCalling_mm10/B_merge.xlsx")
Bcl11a_E<-read.xlsx("ref_data/GSE283392/peakCalling/Bcl11a_Enhancer.xlsx")

Promoter<-Bcl11a %>% filter(annotation=="Promoter") %>% nrow()
Enhancer<-Bcl11a_E %>% filter(annotation!="Promoter")  %>% nrow()


peak_enhancer<-read.table("ref_data/GSE283392/peakCalling/Bcl11a_bound_enhancers.bed")

intersect(Bcl11a$V4,peak_enhancer$V4)


Bcl11a %>% nrow() #所有peak
num_promoter<-Bcl11a[Bcl11a$annotation=="Promoter",] %>% nrow()  #所有promoter
peak_enhancer %>% nrow()

num_co_peak<-intersect(peak_enhancer$V4,Bcl11a %>% filter(annotation=="Promoter") %>% pull(V4)) %>% length()

num_enhancer<-setdiff(peak_enhancer$V4,Bcl11a %>% filter(annotation=="Promoter") %>% pull(V4)) %>% length()

other<-setdiff(Bcl11a$V4,c(Bcl11a[Bcl11a$annotation=="Promoter",]$V4,
                           peak_enhancer$V4) %>% unique()) %>% length()

library(ggplot2)

# create df
data <- data.frame(
  category = c("Promoter", "Enhancer", "Other"),
  value = c(num_promoter,
            num_enhancer,
            other)
)

data$category<-factor(data$category,levels = c("Promoter", "Enhancer", "Other"))
# plot piechart
ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1, color = "white") +  
  coord_polar("y", start = 0) +           
  theme_void() +                          
  labs(fill = "Categories")   +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3,"Set2")))
ggsave('out/Fig.S7F Bcl11a.peak.ratio.pdf',width = 5,height = 5)


#########################################################################################
#### Fig.S7G&H
#########################################################################################

## plot Promoter and enhancer GO pathway
## enhancer annotation
setwd("cuttag-analysis")
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(corrplot)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

## get file (Bcl11a_bound_enhancers_annotated.bed) by run scripts/utils/Bcl11a-enhancer.sh
Bcl11a<-anno_peak("cuttag-analysis/peakCalling_mm10/B_merge_macs3_summits.bed",name = "B_merge")
Bcl11a_Enhancer<-anno_peak("ref_data/GSE283392/peakCalling/Bcl11a_bound_enhancers_annotated.bed",name = "Bcl11a_Enhancer")

Bcl11a<-read.xlsx("cuttag-analysis/peakCalling_mm10/B_merge.xlsx")
Bcl11a_E<-read.xlsx("ref_data/GSE283392/peakCalling/Bcl11a_Enhancer.xlsx")

## BCL11A: enhancer and promoter GO
Bcl11a %>% nrow()
Bcl11a_E %>% nrow()

Promoter<-Bcl11a %>% filter(annotation=="Promoter") %>% pull(SYMBOL) %>% unique()
Enhancer<-Bcl11a_E %>% filter(annotation!="Promoter") %>% pull(SYMBOL) %>% unique()

intersect(Promoter,Enhancer)
write.xlsx(list(overlap=data.frame(gene=intersect(Promoter,Enhancer)),
                Promoter=Bcl11a %>% filter(annotation=="Promoter"),
                Enhancer=Bcl11a_E %>% filter(annotation!="Promoter")),file = "overlap.Promoter-Enhancer.xlsx")

library(clusterProfiler)
library(org.Mm.eg.db)

res<-compareCluster(list(Bcl11a_Promoter=Promoter,
                    Bcl11a_Enhancer=Enhancer),
              fun = "enrichGO",OrgDb='org.Mm.eg.db',keyType ="SYMBOL",ont="BP")
out<-res@compareClusterResult
write.xlsx(out,file = "Bcl11a.Promoter-Enhancer.GO-BP.enrich.xlsx")

## plot
df<-read.xlsx("Bcl11a.Promoter-Enhancer.GO-BP.enrich.xlsx",sheet = "plot")

df.enhancer<-df %>% filter(Class=="Enhancer")
df.promoter<-df %>% filter(Class=="Promoter")


df<- df.promoter %>% arrange(qvalue)
df$Description<-make.unique(df$Description, sep = "_")
df$Description<-factor(df$Description,levels = df$Description)


p1<-ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description))+
  geom_bar(stat = "identity",width = 0.7,fill="salmon")+
  geom_text(
    aes(x = -log10(qvalue) + 0.1),  # 在柱子右侧添加微小偏移
    hjust = 0,                      # 左对齐文本（位于柱子右侧）
    size = 5,                       # 调整文本大小
    nudge_x = 0.05,                 # 水平偏移量（根据数据范围调整）
    check_overlap = FALSE            # 可选：避免标签重叠
  ) +
  theme_classic()+
  ylab("Enriched GO terms")+
  # xlab("-Log10(q value)")+
  xlab(expression(-Log[10]~ (italic(q)~value))) +  # 关键修改：设置斜体q和下标10
  ggtitle("Bcl11a Promoter")+
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_text(size = 12),  # 移除y轴标题
    axis.title.x = element_text(size = 12), 
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),   # 移除y轴刻度标签
    axis.ticks.y = element_blank(),   # 移除y轴刻度线
    plot.margin = margin(0.5,5,0.5,0.5,"cm")
  ) +
  scale_y_discrete(expand = c(0.05, 0.4))  # 调整y轴边距


df<- df.enhancer %>% arrange(qvalue)
df$Description<-make.unique(df$Description, sep = "_")
df$Description<-factor(df$Description,levels = df$Description)


p2<-ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description))+
  geom_bar(stat = "identity",width = 0.7,fill="lightblue")+
  geom_text(
    aes(x = -log10(qvalue) + 0.1),  # 在柱子右侧添加微小偏移
    hjust = 0,                      # 左对齐文本（位于柱子右侧）
    size = 5,                       # 调整文本大小
    nudge_x = 0.05,                 # 水平偏移量（根据数据范围调整）
    check_overlap = FALSE            # 可选：避免标签重叠
  ) +
  theme_classic()+
  ylab("Enriched GO terms")+
  # xlab("-Log10(q value)")+
  xlab(expression(-Log[10]~ (italic(q)~value))) +  # 关键修改：设置斜体q和下标10
  ggtitle("Bcl11a Enhancer")+
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_text(size = 12),  # 移除y轴标题
    axis.title.x = element_text(size = 12), 
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),   # 移除y轴刻度标签
    axis.ticks.y = element_blank(),   # 移除y轴刻度线
    plot.margin = margin(0.5,5,0.5,0.5,"cm")
  ) +
  scale_y_discrete(expand = c(0.05, 0.4))  # 调整y轴边距
p1/p2

ggsave(plot = p1/p2,filename = "Fig.S7G-H Bcl11a.Enhancer&Promoter.GO.pdf",width = 5.5,height = 8)  



#########################################################################################
#### Fig.S7I
#########################################################################################
setwd("cuttag-analysis")
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
library(openxlsx)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

anno_peak<-function(name){
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  
  # filter
  data<-read.table(paste0("peakCalling_mm10/",name,"_macs3_summits.bed"))
  data<-data%>%filter(V1%in%paste0("chr",c(1:19,"X","Y","M")))
  # data[data$V1=="MT",]$V1<-"M"
  data$V1%>%unique()
  write.table(data,file = paste0("peakCalling_mm10/clean_bed/",name,"_macs3_summits.clean.bed"),quote = F,col.names = F,row.names = F,sep = "\t")
  
  peak=readPeakFile(peakfile = paste0("peakCalling_mm10/clean_bed/",name,"_macs3_summits.clean.bed"))

  print(GenomeInfoDb::seqlevels(peak))
  print(GenomeInfoDb::seqlevels(peak)%in%GenomeInfoDb::seqlevels(txdb))
  
  out = annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=5000,annoDb = "org.Mm.eg.db")
  write.xlsx(as.data.frame(out),file = paste0("peakCalling_mm10/anno_bed/",name,".xlsx"))
  write.table(as.data.frame(out),file = paste0("peakCalling_mm10/anno_bed/",name,".bed"),quote = F,col.names = F,row.names = F,sep = "\t")
  df<-as.data.frame(out)
  df$end<-df$end+1
  write.table(df,file = paste0("peakCalling_mm10/anno_bed/",name,".end+1.bed"),quote = F,col.names = F,row.names = F,sep = "\t")
  return(out)
}
B1<-anno_peak("B1")
B2<-anno_peak("B2")
B3<-anno_peak("B3")

# homer data fomate for findMotifs.pl
process<-function(peak,name){
  df<-as.data.frame(peak)
  out<-df%>%dplyr::select(SYMBOL)
  colnames(out)<-"Acc"
  write.table(out,file = paste0("cuttag-analysis/further_analysis/motif/",name,"-findMotifs.txt"),sep = "\t",quote = F,row.names = F)
}

process(B1,"B1")
process(B2,"B2")
process(B3,"B3")

## merge samples
gene<-intersect(as.data.frame(B1)%>%pull(SYMBOL),intersect(as.data.frame(B2)%>%pull(SYMBOL),as.data.frame(B3)%>%pull(SYMBOL)))
df<-data.frame(Acc=gene)
write.table(df,file = "cuttag-analysis/further_analysis/motif/B123-findMotifs.txt",sep = "\t",quote = F,row.names = F)

# HOMER motif enrichment

# mkdir cuttag-analysis/further_analysis/motif/B123-findMotifs
# findMotifs.pl cuttag-analysis/further_analysis/motif/B123-findMotifs.txt \
# 	mouse \
# 	cuttag-analysis/further_analysis/motif/B123-findMotifs -len 8,10,12 -p 30



#########################################################################################
#### Fig.S7J
#########################################################################################

# setwd("~/Project/YuYong_TALL/out")
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(patchwork)


Runx1_mm<-read.xlsx("ref_data/GSE218147/peakCalling_mm10/Runx1.xlsx")
# Runx1_mm<-read.xlsx("ref_data/GSE218147_DN1/peakCalling/DN1_Runx1.xlsx")
Bcl11a_mm1<-read.xlsx("cuttag-analysis/peakCalling/anno_bed/B1.xlsx")
Bcl11a_mm2<-read.xlsx("cuttag-analysis/peakCalling/anno_bed/B2.xlsx")
Bcl11a_mm3<-read.xlsx("cuttag-analysis/peakCalling/anno_bed/B3.xlsx")
Bcl11a_mm<-intersect(intersect(Bcl11a_mm1$SYMBOL,Bcl11a_mm2$SYMBOL),Bcl11a_mm3$SYMBOL)
PU1_mm<-read.xlsx("ref_data/GSE31235/peakCalling/PU.1.xlsx")
Bcl11b_mm<-read.xlsx("ref_data/GSE110305_bcl11b/peakCalling/Bcl11b_mm.xlsx")
gene.mm.all<-c(Runx1_mm$SYMBOL,Bcl11a_mm,PU1_mm$SYMBOL) %>% unique()


library(eulerr)
library(ggpubr)
library(ggvenn)

set2_1<-list(  
  Bcl11a=Bcl11a_mm,
  Runx1=Runx1_mm$SYMBOL %>% unique()
)
set2_2<-list(  
  Bcl11a=Bcl11a_mm,
  PU.1=PU1_mm$SYMBOL %>% unique()
)
set2_3<-list(  
  Bcl11a=Bcl11a_mm,
  Bcl11b=Bcl11b_mm$SYMBOL %>% unique()
)
set4<-list(  
  Bcl11a=Bcl11a_mm,
  Bcl11b=Bcl11b_mm$SYMBOL %>% unique(),
  Runx1=Runx1_mm$SYMBOL %>% unique(),
  PU.1=PU1_mm$SYMBOL %>% unique()
)
p1<-ggvenn(set2_1,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("All gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))


gene<-intersect(set2_1$Bcl11a,set2_1$Runx1)
library(clusterProfiler)
library(org.Mm.eg.db)
gene.id<-bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb =org.Mm.eg.db )
gene.id<-gene.id$ENTREZID %>% unique()

res<-enrichGO(gene.id,OrgDb =org.Mm.eg.db,ont = "ALL",pvalueCutoff = 0.05 )

df<-res@result
openxlsx::write.xlsx(df, file = "data/Bcl11a&Runx1.peak.venn.GO.xlsx")



#####################  
# 绘制富集图
path<-c('GO:1903131',
        'GO:1903706',
        'GO:0030098',
        'GO:0045785',
        'GO:0051098',
        'GO:1902105',
        'GO:0010720',
        'GO:0030217',
        'GO:0002696',
        'GO:0050863')

df<-df %>% filter(ID %in% path)


df$Description<-make.unique(df$Description, sep = "_")
df$Description<-factor(df$Description,levels = df$Description)

# df$Category<-""
# df[df$NES>0,]$Category<-"Up-regulated"
# df[df$NES<0,]$Category<-"Down-regulated"
# df<-df %>% arrange(desc(Category),qvalues)
df<-df %>% arrange(pvalue)

df$Description<-factor(df$Description,levels = df$Description)


ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description))+
  geom_bar(stat = "identity",width = 0.7,fill="salmon")+
  geom_text(
    aes(x = -log10(qvalue) + 0.1),  # 在柱子右侧添加微小偏移
    hjust = 0,                      # 左对齐文本（位于柱子右侧）
    size = 5,                       # 调整文本大小
    nudge_x = 0.05,                 # 水平偏移量（根据数据范围调整）
    check_overlap = FALSE            # 可选：避免标签重叠
  ) +
  theme_classic()+
  ylab("Enriched GO terms")+
  # xlab("-Log10(q value)")+
  xlab(expression(-Log[10]~ (italic(q)~value))) +  # 关键修改：设置斜体q和下标10
  # ggtitle(title,subtitle = subtitle)+
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_text(size = 12),  # 移除y轴标题
    axis.title.x = element_text(size = 12), 
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),   # 移除y轴刻度标签
    axis.ticks.y = element_blank(),   # 移除y轴刻度线
    plot.margin = margin(0.5,5,0.5,0.5,"cm")
  ) +
  scale_y_discrete(expand = c(0.05, 0.4))  # 调整y轴边距
ggsave(filename = "out/Fig.S7J Bcl11a&Runx1.peak.venn.GO.pdf",width = 5.8,height = 6)  


#########################################################################################
#### Fig.S7K
#########################################################################################

bcl11a.deg<-read.xlsx("../VCTB-VCAB-DN2-analysis/further_analysis/VCAB-DKO vs. VCAB.DESeq2.withTPM.xlsx")
bcl11a.deg$sig<-"ns"
bcl11a.up<-bcl11a.deg %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% pull(gene) %>% unique()
bcl11a.dn<-bcl11a.deg %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% pull(gene) %>% unique()

gene2des<-bcl11a.deg %>% dplyr::select(gene,description) %>% distinct(gene,.keep_all = T) %>% tibble::column_to_rownames("gene")

Runx1_mm$SYMBOL
Runx1.deg<-read.xlsx("Venn.xlsx",sheet = "Runx1",startRow = 2,rowNames = T)
runx1.up<-Runx1.deg %>% filter(Note=="OE.UP") %>% rownames()
runx1.dn<-Runx1.deg %>% filter(Note=="OE.DN") %>% rownames()
gene<-intersect(Bcl11a_mm,Runx1_mm$SYMBOL)

library(eulerr)

df<-data.frame(Gene=gene,
               Description=gene2des[gene,]
               # Bcl11a.UP=gene %in% bcl11a.up,
               # Bcl11a.DN=gene %in% bcl11a.dn,
               # Runx1.UP=gene %in% runx1.up,
               # Runx1.DN=gene %in% runx1.dn
               )
df$Bcl11a<-""
df[df$Gene %in% bcl11a.up,]$Bcl11a<-'Up'
df[df$Gene %in% bcl11a.dn,]$Bcl11a<-'Down'
df$Runx1<-""
df[df$Gene %in% runx1.up,]$Runx1<-'Up'
df[df$Gene %in% runx1.dn,]$Runx1<-'Down'


write.xlsx(df,file = "Bcl11a&Runx1.overlap-1450.gene.withDEGinfo.xlsx")



library(VennDiagram)
library(ggplot2)

# 假设你的基因集，改成你自己的向量

df <- read.xlsx("~/Project/YuYong_TALL/cuttag-analysis/Bcl11a&Runx1.overlap-1450.gene.withDEGinfo.xlsx")

#### Bcl11a OE ----

Bcl11a.OE<-read.xlsx("~/Project/YuYong_TALL/ov-analysis/analysis/EGFP.OV-CTL.DESeq2.withTPM.xlsx")
setwd("~/Project/YuYong_TALL/out")


# 读取基因集 (假设都存在文件里，每个文件一列gene symbol)
bcl11a_runx1_overlap <- df$Gene
bcl11a_OE_up <- Bcl11a.OE %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% pull(gene)
bcl11a_OE_down <- Bcl11a.OE %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% pull(gene)
runx1_OE_up <- Runx1.deg %>% filter(Note=="OE.UP") %>% rownames()

# 只考虑共结合的1450基因作为背景
bg <- bcl11a_runx1_overlap

# 两个核心集合：
# 协同激活：Bcl11a OE 上调 ∩ Runx1 OE 上调
coactivation <- intersect(intersect(bcl11a_OE_up, bg), runx1_OE_up)

# 抢占抑制：Bcl11a OE 下调 ∩ Runx1 OE 上调
cofactor_theft <- intersect(intersect(bcl11a_OE_down, bg), runx1_OE_up)

# 输出数量
cat("协同激活 (KO下调 & OE上调):", length(coactivation), "genes\n")
cat("抢占抑制 (KO上调 & OE上调):", length(cofactor_theft), "genes\n")

# Fisher 精确检验函数
do_fisher <- function(setA, setB, background){
  # 先限制集合在背景里，去重
  setA <- unique(intersect(setA, background))
  setB <- unique(intersect(setB, background))
  
  a <- length(intersect(setA, setB))  # 交集
  b <- length(setA) - a               # 仅在A
  c <- length(setB) - a               # 仅在B
  d <- length(background) - (a + b + c)  # 背景剩余
  
  if (d < 0) {
    stop("背景集太小或集合没限制到背景，检查输入")
  }
  
  mat <- matrix(c(a,b,c,d), nrow=2)
  return(fisher.test(mat))
}


# 分别检验协同激活和抢占抑制
test_coactivation <- do_fisher(bcl11a_OE_up, runx1_OE_up, bg)
test_theft <- do_fisher(bcl11a_OE_down, runx1_OE_up, bg)

cat("协同激活 p值:", test_coactivation$p.value, "\n")
cat("抢占抑制 p值:", test_theft$p.value, "\n")

# 绘制韦恩图
pdf(file = "Bcl11a-OE&Runx1-OE.overlap-deg.co-activation.pdf",width = 4,height = 4)
venn_coactivation <- draw.pairwise.venn(
  area1 = length(intersect(bcl11a_OE_up, bg)),
  area2 = length(intersect(runx1_OE_up, bg)),
  cross.area = length(coactivation),
  category = c("Bcl11a OE Up", "Runx1 OE Up"),
  fill = c("#A390E4", "#87CEEB"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
grid.text(
  label = paste0("Bcl11a&Runx1 co-activation p=", signif(test_coactivation$p.value, 3)),
  x = 0.5, y = 0.1,  # 位置（左偏，对应Group_A）
  gp = gpar(fontsize = 10, fontface = "bold", col = "seagreen4")
)
dev.off()

pdf(file = "Bcl11a-OE&Runx1-OE.overlap-deg.cofactor theft.pdf",width = 4,height = 4)

venn_theft <- draw.pairwise.venn(
  area1 = length(intersect(bcl11a_OE_down, bg)),
  area2 = length(intersect(runx1_OE_up, bg)),
  cross.area = length(cofactor_theft),
  category = c("Bcl11a OE Down", "Runx1 OE Up"),
  fill = c("#E6A8D7", "#87CEEB"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)

grid.text(
  label = paste0("Bcl11a&Runx1 cofactor theft p=", signif(test_theft$p.value, 3)),
  x = 0.5, y = 0.1,  # 位置（左偏，对应Group_A）
  gp = gpar(fontsize = 10, fontface = "bold", col = "seagreen4")
)
dev.off()

#### 用新的Runx1的表达数据 ----
# 下载自https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154304

Runx1.deg.new<-read.xlsx("DEseq2.count&tpm.padj0.05.xlsx") %>% data.frame()

bcl11a_runx1_overlap <- df$Gene
bcl11a_KO_up <- bcl11a.deg %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% pull(gene)
bcl11a_KO_down <- bcl11a.deg %>% filter(padj<0.05) %>% filter(log2FoldChange< 0) %>% pull(gene)
runx1_DKO_down <- Runx1.deg.new %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% pull(gene)

# 只考虑共结合的1450基因作为背景
bg <- bcl11a_runx1_overlap

# 两个核心集合：
# 协同激活：Bcl11a KO 下调 ∩ Runx1 OE 上调
coactivation <- intersect(intersect(bcl11a_KO_down, bg), runx1_DKO_down)

# 抢占抑制：Bcl11a KO 上调 ∩ Runx1 OE 上调
cofactor_theft <- intersect(intersect(bcl11a_KO_up, bg), runx1_DKO_down)

# 输出数量
cat("协同激活 (KO下调 & DKO下调):", length(coactivation), "genes\n")
cat("抢占抑制 (KO上调 & DKO下调):", length(cofactor_theft), "genes\n")


# 分别检验协同激活和抢占抑制
test_coactivation <- do_fisher(bcl11a_KO_down, runx1_DKO_down, bg)
test_theft <- do_fisher(bcl11a_KO_up, runx1_DKO_down, bg)

cat("协同激活 p值:", test_coactivation$p.value, "\n")
cat("抢占抑制 p值:", test_theft$p.value, "\n")

# 绘制韦恩图
pdf(file = "Bcl11a-DKO&Runx1-DKO.overlap-deg.co-activation.pdf",width = 4,height = 4)
venn_coactivation <- draw.pairwise.venn(
  area1 = length(intersect(bcl11a_KO_down, bg)),
  area2 = length(intersect(runx1_OE_up, bg)),
  cross.area = length(coactivation),
  category = c("Bcl11a KO Down", "Runx1 DKO Down"),
  fill = c("#A390E4", "#87CEEB"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
grid.text(
  label = paste0("Bcl11a&Runx1 co-activation p=", signif(test_coactivation$p.value, 3)),
  x = 0.5, y = 0.1,  # 位置（左偏，对应Group_A）
  gp = gpar(fontsize = 10, fontface = "bold", col = "seagreen4")
)
dev.off()

pdf(file = "Bcl11a-DKO&Runx1-DKO.overlap-deg.cofactor theft.pdf",width = 4,height = 4)

venn_theft <- draw.pairwise.venn(
  area1 = length(intersect(bcl11a_KO_up, bg)),
  area2 = length(intersect(runx1_OE_up, bg)),
  cross.area = length(cofactor_theft),
  category = c("Bcl11a KO Up", "Runx1 DKO Down"),
  fill = c("#E6A8D7", "#87CEEB"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)

grid.text(
  label = paste0("Bcl11a&Runx1 cofactor theft p=", signif(test_theft$p.value, 3)),
  x = 0.5, y = 0.1,  # 位置（左偏，对应Group_A）
  gp = gpar(fontsize = 10, fontface = "bold", col = "seagreen4")
)
dev.off()

#### 重新给出表格 1450基因----
df <- read.xlsx("~/Project/YuYong_TALL/cuttag-analysis/Bcl11a&Runx1.overlap-1450.gene.withDEGinfo.xlsx")

Bcl11a.OE<-read.xlsx("~/Project/YuYong_TALL/ov-analysis/analysis/EGFP.OV-CTL.DESeq2.withTPM.xlsx")

bcl11a_OE_up <- Bcl11a.OE %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% pull(gene)
bcl11a_OE_down <- Bcl11a.OE %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% pull(gene)
df$Bcl11a.OE<-""
df[df$Gene %in% bcl11a_OE_up,]$Bcl11a.OE<-'Up'
df[df$Gene %in% bcl11a_OE_down,]$Bcl11a.OE<-'Down'

Runx1.deg.new<-read.xlsx("DEseq2.count&tpm.padj0.05.xlsx") %>% data.frame()
Runx1_DKO_up <- Runx1.deg.new %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% pull(gene)
Runx1_DKO_down <- Runx1.deg.new %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% pull(gene)

df$Runx1.DKO<-""
df[df$Gene %in% Runx1_DKO_up,]$Runx1.DKO<-'Up'
df[df$Gene %in% Runx1_DKO_down,]$Runx1.DKO<-'Down'



## 增加promoter enhancer 
Bcl11a%>% dim()
Bcl11a_E%>% dim()
Runx1 %>% dim()
Runx1_E%>% dim()

Bcl11a_promoter<-Bcl11a %>% filter(annotation=="Promoter") %>% pull(SYMBOL) %>% unique()
Bcl11a_enhancer<-Bcl11a_E %>% filter(annotation!="Promoter")%>% pull(SYMBOL) %>% unique()

Runx1_promoter<-Runx1 %>% filter(annotation=="Promoter") %>% pull(SYMBOL) %>% unique()
Runx1_enhancer<-Runx1_E %>% filter(annotation!="Promoter")%>% pull(SYMBOL) %>% unique()

df$Bcl11a.Promoter<-""
df[df$Gene %in% Bcl11a_promoter,]$Bcl11a.Promoter<-'Yes'
df$Bcl11a.Enhancer<-""
df[df$Gene %in% Bcl11a_enhancer,]$Bcl11a.Enhancer<-'Yes'

df$Runx1.Promoter<-""
df[df$Gene %in% Runx1_promoter,]$Runx1.Promoter<-'Yes'
df$Runx1.Enhancer<-""
df[df$Gene %in% Runx1_enhancer,]$Runx1.Enhancer<-'Yes'



write.xlsx(df,file = "Bcl11a&Runx1.overlap-1450.gene.withDEGinfo.2.xlsx")


###韦恩图 三者 Bcl11a-KO, Runx1-DKO, Runx-OE


Bcl11a_KO_up<-df %>% filter(Bcl11a=="Up") %>% pull(Gene)
Bcl11a_KO_dn<-df %>% filter(Bcl11a=="Down") %>% pull(Gene)

Runx1_DKO_up<-df %>% filter(Runx1.DKO=="Up") %>% pull(Gene)
Runx1_DKO_dn<-df %>% filter(Runx1.DKO=="Down") %>% pull(Gene)

Runx_OE_up<-df %>% filter(Runx1=="Up") %>% pull(Gene)
Runx_OE_dn<-df %>% filter(Runx1=="Down") %>% pull(Gene)

library(ggVennDiagram)

# 创建基因列表（示例名称，请替换为实际变量）
gene_lists <- list(
  Bcl11a_KO_up = Bcl11a_KO_up,
  Bcl11a_KO_dn = Bcl11a_KO_dn,
  Runx1_DKO_up = Runx1_DKO_up,
  Runx1_DKO_dn = Runx1_DKO_dn,
  Runx_OE_up = Runx_OE_up,
  Runx_OE_dn = Runx_OE_dn
)

# 绘制韦恩图
# ggVennDiagram(gene_lists, 
#               label_alpha = 0,  # 隐藏区域标签背景
#               set_size = 4,     # 集合名称字号
#               label_size = 3    # 交集数字字号
# ) + 
#   scale_fill_gradient(low = "white", high = "red")  # 颜色渐变
# ggsave("Bcl11a&Runx1.venn_6sets.pdf", width=10, height=8, dpi=300)

plot<-df %>% filter(Bcl11a!=""|Runx1!=""|Runx1.DKO!="") %>% dplyr::select(Gene,Bcl11a,Runx1,Runx1.DKO) %>% tibble::column_to_rownames("Gene")

plot[plot=="Up"]<-1
plot[plot=="Down"]<- -1
plot[plot==""]<-NA
plot.bk<-plot
plot <- data.frame(lapply(plot, as.numeric)) 
rownames(plot)<-rownames(plot.bk)
plot<-plot%>% arrange(Runx1.DKO,Runx1,Bcl11a)
colnames(plot)<-c("Bcl11a-KO","Runx1-OE","Runx1-DKO")

# pdf("Bcl11a&Runx1.3group.union.htp.pdf", width=3.5, height=8)

# pheatmap::pheatmap(plot,cluster_cols = F,cluster_rows = F)
# dev.off()


plot<-df %>% filter(Bcl11a!=""&(Runx1!=""|Runx1.DKO!="")) %>% dplyr::select(Gene,Bcl11a,Runx1,Runx1.DKO) %>% tibble::column_to_rownames("Gene")


plot[plot=="Up"]<-1
plot[plot=="Down"]<- -1
plot[plot==""]<-NA
plot.bk<-plot
plot <- data.frame(lapply(plot, as.numeric)) 
rownames(plot)<-rownames(plot.bk)
plot<-plot%>% arrange(Runx1.DKO,Runx1,Bcl11a)
colnames(plot)<-c("Bcl11a-KO","Runx1-OE","Runx1-DKO")

# pdf("Bcl11a&Runx1.3group.overlap.htp.pdf", width=3.5, height=10)
# pheatmap::pheatmap(plot,cluster_cols = F,cluster_rows = F)
# dev.off()


#### 需要Runx1-OE的上下调和Bcl11a-KO的上下调做韦恩图 ----

library(ggVennDiagram)

# 创建基因列表（示例名称，请替换为实际变量）
gene_lists <- list(
  Bcl11a_KO_up = Bcl11a_KO_up,
  Bcl11a_KO_dn = Bcl11a_KO_dn,
  Runx_OE_up = Runx_OE_up,
  Runx_OE_dn = Runx_OE_dn
)
a<-euler(gene_lists)


pdf("out/Fig.S7K Bcl11a&Runx1.venn4.elipse.pdf", width=5, height=4)
plot(a, 
     labels = TRUE,        # 显示标签
     counts = TRUE,        # 显示数量
     fill = RColorBrewer::brewer.pal(n = 8,name = "Paired")[c(1,3,5,7)],
     alpha = 0.5)
dev.off()
