#########################################################################################
#### Figure 4A
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
run_DESeq2 <- function(key1, key2, out_postfix, tpm, raw_count, write = TRUE) {
    # id2gene<-read.table("/home/zhengjie/DATA/ensmusg2symbol.GRCm39.txt",header = T,sep = "\t",row.names = 1)
    id2gene <- read.table("data/mm10.id2symbol.txt", header = T, sep = "\t", row.names = 1)
    id2description <- read.table("data/ensmusg2description.GRCm39.txt", header = T, sep = "\t", row.names = 1, check.names = F, quote = "")
    sample.re <- read.xlsx("data/sample.xlsx", sheet = "bulk") %>% tibble::column_to_rownames(names(.)[1])

    # sample<-colnames(raw_count)
    # dis<-sample[grep(key1,sample)]
    # ctl<-sample[grep(key2,sample)]
    dis <- key1
    ctl <- key2
    data <- raw_count[, c(dis, ctl)]
    expr <- tpm[, c(dis, ctl)]

    condition <- c(rep("dis", length(dis)), rep("ctl", length(ctl)))
    coldata <- data.frame(row.names = colnames(data), condition)
    print(coldata)
    dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~condition)
    # dds$condition <- relevel(dds$condition,"ctl")
    dds <- DESeq(dds)
    res <- results(dds)
    out <- res[order(res$padj), ] %>% data.frame()
    out$gene <- id2gene[rownames(out), ]
    out$description <- id2description[rownames(out), ]

    out <- out %>% dplyr::select("gene", "description", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
    out <- cbind(out, expr[rownames(out), ])
    if (write) {
        write.xlsx(out, file = paste0("./", out_postfix, ".DESeq2.withTPM.xlsx"), overwrite = TRUE)
    }
    return(out)
}
#### VCAB-DKO vs VCAB ####
dis <- c("VCAB-DKO-1", "VCAB-DKO-2", "VCAB-DKO-3")
ctl <- c("VCAB-1", "VCAB1", "VCAB2", "VCAB3")
res.DKO_VCAB <- run_DESeq2(key1 = dis, key2 = ctl, out_postfix = "VCAB-DKO vs. VCAB.all", tpm = tpm, raw_count = raw_count, write = TRUE)

### selected degs for plot
gene1 <- c(
    # Lymphocyte activation and differentiation signature:
    "Lck", "Cd3e", "Zap70", "Cd4", "Nfil3", "Tcf3", "Cd8a", "Ctla4", "Dtx1", "Ly6d", "Cd3d", "Ccr9"
)
gene2 <- c(
    # Apoptosis signature: ,
    "Bid", "Prkcd", "Fyn", "Dapk1", "Rela", "Bbc3", "Creb3l1", "Dapk3", "Bax"
)
gene3 <- c(
    # Cell adhension signature:
    "Cxcl12", "Nck2", "Cd63", "Ccr2", "Rock1", "Syk", "Fut7", "Smad7", "Ccl2", "Afdn", "Itgb3"
)
gene4 <- c(
    # Progenitors signature:
    "Hoxa9", "Spi1", "Scin", "Cd48", "Cd7", "Bcl11a", "Egr1", "Cd34", "Nfam1", "Hoxa7", "Btk", "Erg"
)

gene <- c(gene1, gene2, gene3, gene4)
deg <- read.xlsx("data/VCAB-DKO vs. VCAB.DESeq2.withTPM.xlsx")
library(tidyverse)
df <- deg %>% filter(gene %in% gene)
df <- expr[gene, c("VCAB3", "VCAB1", "VCAB2", "VCAB-DKO-2", "VCAB-DKO-1", "VCAB-DKO-3")]

pdf("out/Fig.4A VCAB-DKO vs. VCAB.heatmap.pdf", width = 4.2, height = 8)
pheatmap::pheatmap(df, scale = "row", cutree_cols = 2, cluster_cols = F, cluster_rows = F, gaps_col = 3, treeheight_row = 30, treeheight_col = 20)
dev.off()



#########################################################################################
#### Fig.S6A
#########################################################################################

select.path <- read.xlsx("data/VCAB-DKO vs. VCAB.DESeq2.lfc05&padj005.GO.p005.xlsx", sheet = "plot")
select.path$Category <- select.path$Direction

df <- select.path %>%
    group_by(Category) %>%
    arrange(qvalue) %>%
    slice_head(n = 50)
df$Description <- make.unique(df$Description, sep = "_")
df$Description <- factor(df$Description, levels = df$Description)

df[df$Category == "UP", ]$Category <- "Up-regulated"
df[df$Category == "Down", ]$Category <- "Down-regulated"


ggplot(df, aes(y = Description, x = -log10(qvalue), label = Description, fill = Category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(
        aes(x = -log10(qvalue) + 0.1), # 在柱子右侧添加微小偏移
        hjust = 0, # 左对齐文本（位于柱子右侧）
        size = 5, # 调整文本大小
        nudge_x = 0.05, # 水平偏移量（根据数据范围调整）
        check_overlap = FALSE # 可选：避免标签重叠
    ) +
    theme_classic() +
    ylab("Enriched GO terms") +
    # xlab("-Log10(q value)")+
    xlab(expression(-Log[10] ~ (italic(q) ~ value))) + # 关键修改：设置斜体q和下标10
    # ggtitle(title,subtitle = subtitle)+
    theme(
        axis.line.y = element_blank(),
        axis.title.y = element_text(size = 12), # 移除y轴标题
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(), # 移除y轴刻度标签
        axis.ticks.y = element_blank(), # 移除y轴刻度线
        plot.margin = margin(0.5, 5, 0.5, 0.5, "cm")
    ) +
    scale_fill_manual(values = c("#80a492", "salmon")) +
    scale_y_discrete(expand = c(0.05, 0.4)) # 调整y轴边距

ggsave(filename = "out/Fig.4B VCAB-DKO vs. VCAB.GO.pdf", width = 8, height = 5)




#########################################################################################
#### Figure 4C
#########################################################################################

setwd("cuttag-analysis")
library(ChIPseeker)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(patchwork)


Runx1_mm<-read.xlsx("data/Runx1.xlsx")
Bcl11a_mm1<-read.xlsx("data/B1.xlsx")
Bcl11a_mm2<-read.xlsx("data/B2.xlsx")
Bcl11a_mm3<-read.xlsx("cdata/B3.xlsx")
Bcl11a_mm<-intersect(intersect(Bcl11a_mm1$SYMBOL,Bcl11a_mm2$SYMBOL),Bcl11a_mm3$SYMBOL)
PU1_mm<-read.xlsx("data/PU.1.xlsx")
Bcl11b_mm<-read.xlsx("data/Bcl11b_mm.xlsx")
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
ggsave(plot=(p1|p2|p3)/(p1.2|p2.2|p3.2)/(p4|p4.2),filename = "out/Fig.4C Venn of Bcl11a-Runx1-PU.1-Bcl11b.pdf",width = 18,height = 18)



#########################################################################################
#### Figure 4D&E
#########################################################################################
## IGV view of bigwig files:
## 1. CUT&Tag Bcl11A (this study)
# pipline in scripts/utils/Bcl11a-CUT&Tag.sh
## 2. ChIP Runx1(GSE218147)
# GSM6735504 and GSM6735505, compared to GSM6735509 and GSM6735510
# pipline in scripts/utils/Runx1-ChIP.sh
## 3. ChIP H3K27ac(GSE283392)
# GSM8661325
# pipline in scripts/utils/H3K27ac-ChIP-enhancer.sh