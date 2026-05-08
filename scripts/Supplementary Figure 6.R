#########################################################################################
#### Fig.S6A
#########################################################################################

select.path<-read.xlsx("VCAB-DKO vs. VCAB.DESeq2.lfc05&padj005.GO.p005.xlsx",sheet = "plot")
select.path$Category<-select.path$Direction

df<- select.path %>% group_by(Category) %>% arrange(qvalue) %>% slice_head(n=50)
df$Description<-make.unique(df$Description, sep = "_")
df$Description<-factor(df$Description,levels = df$Description)

df[df$Category=="UP",]$Category<-"Up-regulated"
df[df$Category=="Down",]$Category<-"Down-regulated"


ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description,fill=Category))+
  geom_bar(stat = "identity",width = 0.7)+
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
  scale_fill_manual(values = c("#80a492","salmon"))+
  scale_y_discrete(expand = c(0.05, 0.4))  # 调整y轴边距

ggsave(filename = "out/Fig.S6A VCAB-DKO vs. VCAB.GO.pdf",width = 8,height = 5)  


#########################################################################################
#### Fig.S6D
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


peak_enhancer<-read.table("/home/zhengjie/Project/YuYong_TALL/ref_data/GSE283392/peakCalling/Bcl11a_bound_enhancers.bed")

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
ggsave('out/Fig.S6D Bcl11a.peak.ratio.pdf',width = 5,height = 5)


#########################################################################################
#### Fig.S6E
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

ggsave(plot = p1/p2,filename = "Fig.S6E Bcl11a.Enhancer&Promoter.GO.pdf",width = 5.5,height = 8)  


#########################################################################################
#### Fig.S6F
#########################################################################################

## IGV view of bigwig files:
# CUT&Tag Bcl11A(this study)
# ChIP Runx1(GSE218147) 
# ChIP H3K27ac(GSE283392)
## a demo pipline in scripts/utils/ChIP-code.sh
# note: the demo is for tcf7
