#########################################################################################
#### Fig.S9B
#########################################################################################
setwd("ov-analysis/analysis/")

library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(corrplot)
library(openxlsx)
library(org.Mm.eg.db)
library(eulerr)


deg_OV_egfp_OV_td<-read.xlsx("OV.EGFP-TD.DESeq2.withTPM.xlsx") %>% dplyr::select(gene,log2FoldChange,padj)

OV_egfp_td.up<-deg_OV_egfp_OV_td %>% filter(log2FoldChange>0.1) %>% filter(padj<0.05) %>% pull(gene)
OV_egfp_td.dn<-deg_OV_egfp_OV_td %>% filter(log2FoldChange< -0.1) %>% filter(padj<0.05) %>% pull(gene)

deg_ctl_egfp_ctl_td<-read.xlsx("CTL.EGFP-TD.DESeq2.withTPM.xlsx") %>% dplyr::select(gene,log2FoldChange,padj)

ctl_egfp_td.up<-deg_ctl_egfp_ctl_td %>% filter(log2FoldChange>0.1) %>% filter(padj<0.05) %>% pull(gene)
ctl_egfp_td.dn<-deg_ctl_egfp_ctl_td %>% filter(log2FoldChange< -1) %>% filter(padj<0.05) %>% pull(gene)


## venn
library(eulerr)

ggvenn::ggvenn(list(
  `CTL.EGFP-TD.up`=ctl_egfp_td.up,
  `CTL.EGFP-TD.down`=ctl_egfp_td.dn,
  `OV.EGFP-TD.up`=OV_egfp_td.up,
  `OV.EGFP-TD.down`=OV_egfp_td.dn),
  fill_alpha = 0.7,stroke_alpha = 0.5,
  fill_color = RColorBrewer::brewer.pal(8,"Set3")[2:5]
)+
  theme(plot.margin = margin(0.5,2,0.5,2,"cm"))

# ggsave(filename = "OV-CTL.EGFP-TD.deg.lfc0.1padj0.05.venn.pdf",width = 8,height = 8)

intersect(OV_egfp_td.dn,ctl_egfp_td.up)
intersect(OV_egfp_td.dn,ctl_egfp_td.dn)
intersect(OV_egfp_td.up,ctl_egfp_td.dn)
intersect(OV_egfp_td.up,ctl_egfp_td.up) %>% write.csv("tmp.csv",quote = F)


#### enrichment  ----
up<-setdiff(OV_egfp_td.up,c(ctl_egfp_td.dn,ctl_egfp_td.up))

dn<-setdiff(OV_egfp_td.dn,c(ctl_egfp_td.dn,ctl_egfp_td.up))

library(clusterProfiler)
library(org.Mm.eg.db)

dn.go<-enrichGO(dn,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,pAdjustMethod = "BH")
up.go<-enrichGO(up,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,pAdjustMethod = "BH")
up.res<-up.go@result
dn.res<-dn.go@result

up.sig<-up.res %>% filter(pvalue<0.05) %>% filter(p.adjust<0.05)
dn.sig<-dn.res %>% filter(pvalue<0.05) %>% filter(p.adjust<0.05)

write.xlsx(list(up=up.sig,dn=dn.sig),file = "OV.EGFP-TD.DEG.overlap.GO.enrichment.padj0.05.0603.xlsx",quote = FALSE)

## selected pathway 

path<-read.xlsx("OV.EGFP-TD.DEG.overlap.GO.enrichment.padj0.05.0603.xlsx",sheet = "plot")
path.up<-path %>% filter(Class=="UP") %>% arrange(qvalue)

path.dn_T<-path %>% filter(Category=="T-cell differentiation and activation") %>% arrange(qvalue)
path.dn_MB<-path %>% filter(Category %in% c("Myeloid","B cell")) %>% arrange(qvalue)


plot_enrich<-function(df,fill="grey10",title,subtitle){
  df$Description<-factor(df$Description,levels = df$Description)
  ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description))+
    geom_bar(stat = "identity",width = 0.7,fill=fill)+
    geom_text(
      aes(x = -log10(qvalue) + 0.1),  
      hjust = 0,                     
      size = 5,                       
      nudge_x = 0.05,                 
      check_overlap = FALSE           
    ) +
    theme_classic()+
    ylab("Enriched GO terms")+
    # xlab("-Log10(q value)")+
    xlab(expression(-Log[10]~ (italic(q)~value))) +  
    ggtitle(title,subtitle = subtitle)+
    theme(
      axis.line.y = element_blank(),
      axis.title.y = element_text(size = 12),  
      axis.title.x = element_text(size = 12), 
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),  
      axis.ticks.y = element_blank(),  
      plot.margin = margin(0.5,5,0.5,0.5,"cm")
    ) +
    scale_y_discrete(expand = c(0.05, 0.4))  
  
}
# p1<-plot_enrich(path.up,fill = "grey20",title = "Bcl11a OE;DN2 vs. CTL;DN2",subtitle = "Up-regulated genes")
# p2<-plot_enrich(path.dn_T,fill = "salmon",title = "Bcl11a OE;DN2 vs. CTL;DN2",subtitle = "Down-regulated genes")
# p3<-plot_enrich(path.dn_MB,fill = "#aba3fb",title = "Bcl11a OE;DN2 vs. CTL;DN2",subtitle = "Down-regulated genes")
# ggsave(plot = p1+p2+p3,filename = "OV.EGFP-TD.deg.overlap.GO.0604.pdf",width = 15,height = 4)  
 
df<-path %>% group_by(Category)  %>% arrange(qvalue) %>% slice_head(n=50)
df$Description<-factor(df$Description,levels = df$Description)

df[df$Category=="NULL",]$Category<-"Up-regulated"
df[df$Category=="B cell",]$Category<-"B cell & Myeloid"
df[df$Category=="Myeloid",]$Category<-"B cell & Myeloid"
df$Category<-factor(df$Category,levels = c("T-cell differentiation and activation",
                                           "B cell & Myeloid",
                                           "Up-regulated"))
ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description,fill=Category))+
  geom_bar(stat = "identity",width = 0.7)+
  geom_text(
    aes(x = -log10(qvalue) + 0.1),  
    hjust = 0,                      
    size = 5,                       
    nudge_x = 0.05,                 
    check_overlap = FALSE            
  ) +
  theme_classic()+
  ylab("Enriched GO terms")+
  # xlab("-Log10(q value)")+
  xlab(expression(-Log[10]~ (italic(q)~value))) +  
  # ggtitle(title,subtitle = subtitle)+
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),   
    plot.margin = margin(0.5,5,0.5,0.5,"cm")
  ) +
  scale_y_discrete(expand = c(0.05, 0.4))  
ggsave(filename = "out/Fig.S9B OV.EGFP-TD.deg.overlap.GO.pdf",width = 8,height = 8)  




#########################################################################################
#### Fig.S9C
#########################################################################################
setwd("ov-analysis/analysis/")

library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(corrplot)
library(openxlsx)
library(org.Mm.eg.db)
library(eulerr)

sample.re<-read.xlsx("sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])

data<-read.table("../gene-TPM-matrix.txt",check.names = F)
colnames(data)
colname.new<-lapply(colnames(data),function(x){strsplit(strsplit(x,"ov-analysis/mapping/")[[1]][2],"_")[[1]][1]})%>%as.character()
colnames(data)<-colname.new
colnames(data)<-sample.re[colnames(data),]

tpm<-data

loc<-rowSums(tpm)!=0
data.filter<-tpm[loc,]
data.filter<-log2(data.filter+0.1)

id2gene<-read.table("mm10.id2symbol.txt",header = T,sep = "\t",row.names = 1)

data.filter$gene<-id2gene[rownames(data.filter),]

df<-data.filter %>% distinct(gene,.keep_all = T) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("gene")

data<-read.xlsx("TD.OV-CTL.DESeq2.withTPM.xlsx")
# deg1<-data %>% filter(log2FoldChange>1&padj<0.05) %>% arrange(desc(log2FoldChange)) %>% head(30)%>% pull(gene) 
# deg2<-data %>% filter(log2FoldChange< -1&padj<0.05) %>% arrange(desc(log2FoldChange)) %>% head(30)%>% pull(gene) 
# deg<-c(deg1,deg2)

#### enrichment----

up<-data %>% filter(log2FoldChange>0.1&padj<0.05) %>% pull(gene)

dn<-data %>% filter(log2FoldChange< -0.1&padj<0.05)%>% pull(gene)

library(clusterProfiler)
library(org.Mm.eg.db)

# up.id<-bitr(geneID = up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
# dn.id<-bitr(geneID = dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)

dn.go<-enrichGO(dn,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,pAdjustMethod = "BH")
up.go<-enrichGO(up,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP",pvalueCutoff = 1,pAdjustMethod = "BH")
up.res<-up.go@result
dn.res<-dn.go@result

up.sig<-up.res %>% filter(pvalue<0.05) %>% filter(p.adjust<0.05)
dn.sig<-dn.res %>% filter(pvalue<0.05) %>% filter(p.adjust<0.05)

write.xlsx(list(up=up.sig,dn=dn.sig),file = "TD.OV-CTL.DEG.overlap.GO.enrichment.padj0.05.0607.xlsx",quote = FALSE)

## selected pathway

path<-read.xlsx("TD.OV-CTL.DEG.overlap.GO.enrichment.padj0.05.0607.xlsx",sheet = "plot")


df<- path %>% group_by(Category) %>% arrange(qvalue) %>% slice_head(n=50)
df$Description<-make.unique(df$Description, sep = "_")
df$Description<-factor(df$Description,levels = df$Description)

df[df$Category=="UP",]$Category<-"Up-regulated"
df[df$Category=="Down",]$Category<-"Down-regulated"


ggplot(df,aes(y=Description,x=-log10(qvalue),label=Description,fill=Category))+
  geom_bar(stat = "identity",width = 0.7)+
  geom_text(
    aes(x = -log10(qvalue) + 0.1),  
    hjust = 0,                      
    size = 5,                       
    nudge_x = 0.05,                
    check_overlap = FALSE            
  ) +
  theme_classic()+
  ylab("Enriched GO terms")+
  # xlab("-Log10(q value)")+
  xlab(expression(-Log[10]~ (italic(q)~value))) +  #
  # ggtitle(title,subtitle = subtitle)+
  theme(
    axis.line.y = element_blank(),
    axis.title.y = element_text(size = 12),  
    axis.title.x = element_text(size = 12), 
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),   
    plot.margin = margin(0.5,5,0.5,0.5,"cm")
  ) +
  scale_y_discrete(expand = c(0.05, 0.4))  

ggsave(filename = "out/Fig.S9C TD.OV-CTL.deg.overlap.GO.pdf",width = 6.5,height = 6)  

