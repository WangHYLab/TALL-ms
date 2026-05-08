#########################################################################################
#### Figure 5A 
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

## function to run DESeq2, get compare.DESeq2.withTPM.xlsx
run_DESeq2<-function(key1,key2,out_postfix,tpm,raw_count,write=TRUE){
  # id2gene<-read.table("/home/zhengjie/DATA/ensmusg2symbol.GRCm39.txt",header = T,sep = "\t",row.names = 1)
  id2gene<-read.table("/home/zhengjie/DATA/mm10.id2symbol.txt",header = T,sep = "\t",row.names = 1)
  id2description<-read.table("/home/zhengjie/DATA/ensmusg2description.GRCm39.txt",header = T,sep = "\t",row.names = 1,check.names = F,quote = "")
  sample.re<-read.xlsx("../../ref_data/sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])
  
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
    write.xlsx(out,file = paste0("./",out_postfix,".DESeq2.withTPM.xlsx"),overwrite = TRUE)
  }
  return(out)
}
#### VCAB-DKO vs VCAB ####
dis<-c("VCAB-DKO-1", "VCAB-DKO-2", "VCAB-DKO-3")
ctl<-c("VCAB-1" , "VCAB1", "VCAB2" ,"VCAB3")
res.DKO_VCAB<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "VCAB-DKO vs. VCAB.all",tpm = tpm,raw_count = raw_count,write = TRUE)

### selected degs for plot
gene1<-c(
  # Lymphocyte activation and differentiation signature: 
   'Lck','Cd3e', 'Zap70', 'Cd4','Nfil3', 'Tcf3',  'Cd8a','Ctla4', 'Dtx1','Ly6d', 'Cd3d','Ccr9'
)
gene2<-c(
  # Apoptosis signature: ,
   'Bid', 'Prkcd', 'Fyn','Dapk1',  'Rela', 'Bbc3', 'Creb3l1','Dapk3','Bax'
)
gene3<-c(
  # Cell adhension signature: 
  'Cxcl12', 'Nck2','Cd63', 'Ccr2', 'Rock1','Syk', 'Fut7', 'Smad7', 'Ccl2',   'Afdn', 'Itgb3'
)
gene4<-c(
  # Progenitors signature: 
  'Hoxa9',  'Spi1', 'Scin', 'Cd48', 'Cd7', 'Bcl11a', 'Egr1', 'Cd34', 'Nfam1', 'Hoxa7', 'Btk', 'Erg'
)

gene<-c(gene1,gene2,gene3,gene4)
deg<-read.xlsx("VCAB-DKO vs. VCAB.DESeq2.withTPM.xlsx")
library(tidyverse)
df<-deg %>% filter(gene %in% gene)
df<-expr[gene,c( "VCAB3","VCAB1","VCAB2","VCAB-DKO-2","VCAB-DKO-1","VCAB-DKO-3")]

pdf("out/Fig.5A VCAB-DKO vs. VCAB.heatmap.pdf",width = 4.2,height = 8)
pheatmap::pheatmap(df,scale = "row",cutree_cols = 2,cluster_cols = F,cluster_rows = F,gaps_col = 3,treeheight_row = 30,treeheight_col = 20)
dev.off()


#########################################################################################
#### Figure 5B
#########################################################################################
## IGV view of bigwig files:
# CUT&Tag Bcl11A(this study)
# ChIP Runx1(GSE218147) 
# ChIP H3K27ac(GSE283392)
## a demo pipline in scripts/utils/ChIP-code.sh
# note: the demo is for tcf7


#########################################################################################
#### Figure 5C
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


## RUN HOMER motif analysis in shell:
# findMotifs.pl cuttag-analysis/further_analysis/motif/B123-findMotifs.txt \
# 	mouse \
# 	cuttag-analysis/further_analysis/motif/B123-findMotifs -len 8,10,12 -p 30



#########################################################################################
#### Figure 5D
#########################################################################################

setwd("cuttag-analysis")
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(patchwork)


Runx1_mm<-read.xlsx("ref_data/GSE218147/peakCalling_mm10/Runx1.xlsx")
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

p2<-ggvenn(set2_2,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("All gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))
p3<-ggvenn(set2_3,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("All gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))
p4<-ggvenn(set4,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("All gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))

## PROMOTER
Bcl11a_mm.promoter<-intersect(intersect(Bcl11a_mm1%>% filter(annotation=="Promoter") %>% pull(SYMBOL),
                                        Bcl11a_mm2%>% filter(annotation=="Promoter") %>% pull(SYMBOL)),
                              Bcl11a_mm3%>% filter(annotation=="Promoter") %>% pull(SYMBOL))
Runx1_mm.promoter<-Runx1_mm %>% filter(annotation=="Promoter") %>% pull(SYMBOL) %>% unique()
Bcl11b_mm.promoter<-Bcl11b_mm %>% filter(annotation=="Promoter") %>% pull(SYMBOL) %>% unique()
PU1_mm.promoter<-PU1_mm %>% filter(annotation=="Promoter") %>% pull(SYMBOL) %>% unique()

set2_1<-list(  
  Bcl11a=Bcl11a_mm.promoter,
  Runx1=Runx1_mm.promoter
)
set2_2<-list(  
  Bcl11a=Bcl11a_mm.promoter,
  PU.1=PU1_mm.promoter
)
set2_3<-list(  
  Bcl11a=Bcl11a_mm.promoter,
  Bcl11b=Bcl11b_mm.promoter
)

set4<-list(  
  Bcl11a=Bcl11a_mm.promoter,
  Bcl11b=Bcl11b_mm.promoter,
  Runx1=Runx1_mm.promoter,
  PU.1=PU1_mm.promoter
)
p1.2<-ggvenn(set2_1,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("Annotated in Promoter gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))
p2.2<-ggvenn(set2_2,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("Annotated in Promoter gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))
p3.2<-ggvenn(set2_3,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("Annotated in Promoter gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))
p4.2<-ggvenn(set4,fill_alpha = 0.7,stroke_alpha = 0.5,
           fill_color = RColorBrewer::brewer.pal(8,"Set3")
)+ggtitle("Annotated in Promoter gene")+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))
ggsave(plot=(p1|p2|p3)/(p1.2|p2.2|p3.2)/(p4|p4.2),filename = "out/Fig.5D Venn of Bcl11a-Runx1-PU.1-Bcl11b.pdf",width = 18,height = 18)
