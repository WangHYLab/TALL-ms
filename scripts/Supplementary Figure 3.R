#########################################################################################
#### Fig.S3A
#########################################################################################

CD8<-read.xlsx("data/L-N.CD8.DESeq2.withTPM.xlsx")
expr.cd8<-CD8  %>% distinct(gene,.keep_all = T) %>% tibble::column_to_rownames("gene")

expr.cd8<-log2(expr.cd8[,c(8:16)]+1)

expr<-expr.cd8

annocol<-data.frame(name=colnames(expr),
                    Condition=lapply(colnames(expr),function(x){str_sub(x,1,1)[[1]]}) %>%as.character()
                    # Celltype=lapply(colnames(expr),function(x){strsplit(x,split = "_")[[1]][2]}) %>%as.character()
) %>% tibble::column_to_rownames(colnames(.)[1])

library(ComplexHeatmap)
library(circlize)
library(seriation)
{
  gene.select<-c('Hes1',
                 'Ccnd1',
                 'Notch1',
                 'Dtx1',
                 'Lef1',
                 'Gata3',
                 'Phf6',
                 'Myb',
                 'Myc',
                 'Stat5a',
                 # 'Ezh2',
                 'Cxcr4',
                 'Cdkn2a',
                 'Lyl1',
                 'Lmo1',
                 'Lmo2',
                 'Cdk4',
                 'Cdk6',
                 'Bmi1',
                 'Ezh2',
                 'Pdgfrb',
                 # Down(15)
                 'Tet1',
                 'Tet2',
                 'Tet3',
                 'Dnm2',
                 'Trbv1',
                 'Trav19',
                 'Trav14-2',
                 'Bcl11b',
                 'Wt1',
                 'Cdkn1b',
                 'Spi1',
                 'Ikzf1',
                 'Dusp6',
                 'Socs1',
                 'Socs3')}
genes_to_show <- gene.select

expr.plot<-expr[c(
                  # gene.select,
                  
                  CD8 %>% filter(padj<0.001&abs(log2FoldChange)>0.5) %>% 
                    # head(500) %>% 
                    pull(gene)
) %>% unique(),]

# 构造 row_labels 向量
row_labels <- ifelse(rownames(expr.plot) %in% genes_to_show, rownames(expr.plot), "")

# 绘制热图
df <-scale(t(expr.plot)) %>% t() %>% na.omit()

ht <- Heatmap(
  df,
  name = "Scaled\nexpression",
  # row_order = get_order(o1), column_order = get_order(o2),
  # row_labels = row_labels,  # 仅显示指定基因的行名
  row_km = 2,               # 可选：行聚类分组
  column_km = 2,            # 可选：列聚类分组
  cluster_rows = T,
  cluster_columns = T,
  clustering_method_columns = "ward.D2",
  # clustering_method_rows = "ward.D2",
  show_column_names = TRUE,
  show_row_names = TRUE,    # 确保行名显示开启
  col = colorRamp2(c(-2, 0, 2), c("navy","white","firebrick3")),
  # col = colorRampPalette(c("navy","white","firebrick3"))(100),
  row_gap = unit(1.5, "mm"), column_gap = unit(1.5, "mm"), border = TRUE,
  top_annotation = columnAnnotation(df = annocol,
                                    col=list(Condition=c("L"="#9BC900","N"="#F69389"))),
)
genes <- c(genes_to_show
           # CD8 %>% filter(abs(log2FoldChange)>3)%>% head(20) %>% pull(gene)
           ) %>% unique()
genes <- as.data.frame(genes)
ht<-ht + rowAnnotation(link = anno_mark(at = which(rownames(expr.plot) %in% genes$genes), 
                                        labels = rownames(expr.plot)[which(rownames(expr.plot) %in% genes$genes)], 
                                        labels_gp = gpar(fontsize = 10)))
ht_opt$message = FALSE
# col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
# ht= HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
# 显示热图
pdf("Fig.S3A L-N.DEG.heatmap", width = 6.5, height = 8)
draw(ht)
dev.off()



#########################################################################################
#### Fig.S3B
#########################################################################################
# setwd("rna-analysis/analysis/")
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(corrplot)
library(openxlsx)
library(ggrepel)
library(enrichplot)
library(cowplot)

sample.re<-read.xlsx("data/sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])

### PCA 相关性 确认样本----
data<-read.table("../gene-TPM-matrix.mm10.txt",check.names = F)
colname.new<-lapply(colnames(data),function(x){strsplit(strsplit(x,"rna-analysis/mapping_mm10/")[[1]][2],".genes.results")[[1]][1]})%>%as.character()
colnames(data)<-colname.new
colnames(data)<-sample.re[colnames(data),]
tpm<-data
#去除未表达gene
data<-data[rowSums(data)!=0,]

## 样本一致性----
corr<-cor(data,method = "spearman")
# pdf("Nomal-Leukemia.corrplot.pdf",width = 6.4,height = 6)
pheatmap::pheatmap(corr,clustering_method = "ward.D2",border_color = NA,
                   color =colorRampPalette(c("skyblue", "salmon","red"))(100) ,cutree_rows = 5,cutree_cols = 5)
# dev.off()
# 存在异质性，尤其是样本间
## PCA----
## pca plot 
# data<-log2(tpm+1)
samples<-colnames(tpm)
sample_use<-setdiff(samples,c("N_CD8_0","N_DP_0","L0_CD8_1", "L0_DP_1",  "L0_CD8_2", "L0_DP_2"))


data<-tpm[,sample_use] #排除L0样本
# data<-count
loc<-rowSums(data)!=0
data.filter<-data[loc,]
{
  pca<-prcomp(t(data.filter),scale. = T)
  head(pca$x)
  summary(pca)
  pca_result<-as.data.frame(pca$x[,1:3])
  sample<-colnames(data.filter)
  shape<-lapply(sample, function(x){strsplit(x,split = "_")[[1]][1]})%>%as.character()
  group<-lapply(sample, function(x){strsplit(x,split = "_")[[1]][2]})%>%as.character()
  pca_result<-cbind(pca_result,sample,shape,group)
  
  # 2d pca
  colnames(pca_result)<-c("PC1" ,   "PC2",    "PC3"  ,  "name", "sample" , "class")
  pca_result[pca_result$sample!="N",]$sample<-"L"

  # ONLY CD8
  library(ggforce)
  pca_result_cd8<-pca_result %>% filter(class=="CD8")
  pca_result_DP_N<-pca_result %>% filter(class=="DP"&sample=="N")
  plot<-rbind(pca_result_cd8,pca_result_DP_N)
  p<-ggplot(plot,aes(x=PC1,y=PC2,shape=class,color=sample))+
    geom_point(size=5)+
    geom_mark_ellipse(aes(color=sample),linetype="dashed")+
    geom_text_repel(aes(label=name,x=plot[,1],y=plot[,2]),data=plot,size=2)+
    theme_bw()+
    theme(
      plot.margin = unit(c(1,1,1,1),'lines'),
      panel.grid = element_blank())+
    ylim(c(-150,100))+
    xlim(c(-100,160))+
  labs(x="PCA 1",y="PCA 2")
  
  ggsave(plot=p,"out/Fig.S3B L-N.CD8-N.PCA.pdf",width = 5.5,height = 3.8)
}



#########################################################################################
#### Fig.S3C
#########################################################################################
CD8 <- read.xlsx("data/L-N.CD8.DESeq2.withTPM.xlsx")

df <- CD8 %>%
    dplyr::select(gene, log2FoldChange) %>%
    distinct(gene, .keep_all = T) %>%
    arrange(desc(log2FoldChange)) %>%
    na.omit()
geneList <- setNames(df$log2FoldChange, df$gene)


library(clusterProfiler)
library(msigdbr)
library(enrichplot)

path <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>% dplyr::select(gs_name, gene_symbol)
trans <- path$gene_symbol %>%
    unique() %>%
    human2mouse()

path <- merge(path, trans, all.x = T, by.x = "gene_symbol", by.y = "humanGene")

path <- path %>%
    dplyr::select(gs_name, mouseGene) %>%
    na.omit()


res <- clusterProfiler::GSEA(geneList, TERM2GENE = path, pvalueCutoff = 1)

gsea <- res@result
write.xlsx(res@result, file = "data/L-N.CD8.GSEA-KEGG.xlsx")

path_id <- "KEGG_CELL_CYCLE"
plot_gsea <- function(path_id) {
    title <- paste0(
        path_id, "\n",
        "NES = ", res@result %>% filter(ID == path_id) %>% pull(NES) %>% signif(., 2), "\t",
        "p-value = ", res@result %>% filter(ID == path_id) %>% pull(pvalue) %>% signif(., 2), "\t",
        "adjusted p-value = ", res@result %>% filter(ID == path_id) %>% pull(p.adjust) %>% signif(., 2)
    )
    gseaplot2(x = res, geneSetID = path_id, title = title, pvalue_table = F, color = "darkred", subplots = 1:2)
    ggsave(filename = paste0("out/Fig.2D ", path_id, ".pdf"), width = 5.5, height = 3.2)
}

lapply(c(
    "KEGG_CELL_CYCLE", # 1
    "KEGG_CHEMOKINE_SIGNALING_PATHWAY", # 3
    "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
    "HALLMARK_APOPTOSIS",
    "KEGG_NOTCH_SIGNALING_PATHWAY",
    "KEGG_HEDGEHOG_SIGNALING_PATHWAY"
), FUN = plot_gsea)

## HALLMARK ----
path_MH <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, gene_symbol)
res2 <- clusterProfiler::GSEA(geneList, TERM2GENE = path_MH, pvalueCutoff = 1)

gsea2 <- res2@result
write.xlsx(res2@result, file = "data/rna-analysis/analysis/L-N.CD8.GSEA-HALLMARK.xlsx")

plot_gsea2 <- function(path_id) {
    title <- paste0(
        path_id, "\n",
        "NES = ", res2@result %>% filter(ID == path_id) %>% pull(NES) %>% signif(., 2), "\t",
        "p-value = ", res2@result %>% filter(ID == path_id) %>% pull(pvalue) %>% signif(., 2), "\t",
        "adjusted p-value = ", res2@result %>% filter(ID == path_id) %>% pull(p.adjust) %>% signif(., 2)
    )
    gseaplot2(x = res2, geneSetID = path_id, title = title, pvalue_table = F, color = "darkred", subplots = 1:2)
    ggsave(filename = paste0("out/Fig.S3C ", path_id, ".pdf"), width = 5.5, height = 3.2)
}
lapply("HALLMARK_MTORC1_SIGNALING", FUN = plot_gsea2)
lapply(gsea2$ID[1:20], FUN = plot_gsea2)
