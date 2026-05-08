#########################################################################################
#### Figure 2A 
#########################################################################################

# setwd("WES-analysis")

library(ggplot2)
library(tidyverse)
library(vcfR)
library(openxlsx)

## 读取所有样本并整理注释合并
files<-list.files("out/anno",pattern = ".vcf")
filter=TRUE
merge.vcf<-c()
for (file in files) {
  out_vcf<-c()
  vcf<-read.vcfR( paste0("out/anno/",file),verbose = T)
  tidy_vcf <- vcfR2tidy(vcf)
  
  VCF<-tidy_vcf$fix
  if(filter){
    VCF<-tidy_vcf$fix %>% filter(FILTER=="PASS")
  }
  
  VCF<-VCF %>% select(CHROM,POS,ID,REF,ALT,FILTER,ANN)
  VCF$ID<-paste0(VCF$CHROM,"_",VCF$POS,"_",VCF$REF,"_",VCF$ALT)
  VCF$SAMPLE<-file
  pb <- txtProgressBar(style=3)
  for (line in 1:nrow(VCF)) {
    setTxtProgressBar(pb, line/nrow(VCF),title = paste0("Process sample: ",file))
    ANNO_ori<-VCF[line,]$ANN
    # anno:
    # Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
    anno_split<-strsplit(ANNO_ori,",")[[1]]
    split_strings <- lapply(anno_split, function(x){strsplit(x,split = "\\|")[[1]]})
    # max_cols <- max(sapply(split_strings, length))
    max_cols <-16
    anno_mtx <- do.call(rbind, lapply(split_strings, function(x) c(x, rep(NA, max_cols - length(x)))))
    colnames(anno_mtx)<-c('Allele',
                          'Annotation',
                          'Annotation_Impact', 
                          'Gene_Name', 
                          'Gene_ID', 
                          'Feature_Type', 
                          'Feature_ID', 
                          'Transcript_BioType', 
                          'Rank', 
                          'HGVS.c', 
                          'HGVS.p', 
                          'cDNA.pos / cDNA.length',
                          'CDS.pos / CDS.length', 
                          'AA.pos / AA.length', 
                          'Distance', 
                          'ERRORS / WARNINGS / INFO')
    rep_len<-nrow(anno_mtx)
    vcf<-c()
    for (i in 1:rep_len) {
      vcf<-rbind(vcf,VCF[line,])
    }
    out_vcf<-rbind(out_vcf,cbind(vcf,anno_mtx))
    
  }
  close(pb)
  merge.vcf<-rbind(merge.vcf,out_vcf)
}
write.xlsx(merge.vcf,file = "further_analysis/merge.filtered.vcf")

merge.vcf<-read.xlsx("further_analysis/merge.filtered.vcf.xlsx")

## 挑选driver gene mutation landscape
driver_gene<-openxlsx::read.xlsx("further_analysis/28671688-supp.xlsx",sheet = "Table S9 Driver Mutations")

library(homologene)
gene.use<-human2mouse(driver_gene$Gene)
gene.use<-gene.use$mouseGene %>% unique()

# gene_mut<-read.table("out/anno/2021.11.19VAFT12.2e2RThymus.snpEff_summary.genes.txt")

vcf<-merge.vcf %>% filter(Gene_Name %in% gene.use) # 关注基因
vcf<-merge.vcf # 所有基因 

uni_ID<-vcf$ID %>% unique()
uni_gene<-vcf$Gene_Name %>% unique()
type<-vcf$Annotation_Impact %>% unique()
order<-c("HIGH","MODERATE","LOW","MODIFIER")
df<-data.frame()
for (gene in uni_gene) {
  mini_df<-vcf %>% filter(Gene_Name==gene)
  
  mini_df$Annotation_Impact<-factor(mini_df$Annotation_Impact, levels = order, ordered = TRUE)
  mini_df<- mini_df[order(mini_df$Annotation_Impact), ]
  
  mut_num<-mini_df$ID %>% unique() %>% length()
  sample_num<-mini_df$SAMPLE %>% unique() %>% length()
  
  d<-mini_df %>% dplyr::select(ID,Gene_Name,Annotation_Impact) %>% distinct(ID,.keep_all = T)
  high_num<-d %>% filter(Annotation_Impact=="HIGH") %>% nrow()
  moderate_num<-d %>% filter(Annotation_Impact=="MODERATE") %>% nrow()
  low_num<-d %>% filter(Annotation_Impact=="LOW") %>% nrow()
  modifier_num<-d %>% filter(Annotation_Impact=="MODIFIER") %>% nrow()
  
  d2<-mini_df %>% dplyr::select(ID,SAMPLE) %>% unique.data.frame()
  mut_more1_per_sample<-sum((table(d2$ID,d2$SAMPLE) %>% colSums())>1)
  mut_1_per_sample<-sum((table(d2$ID,d2$SAMPLE) %>% colSums())==1)
  
  row<-c(gene,mut_num,sample_num,
         high_num,moderate_num,low_num,modifier_num,
         mut_1_per_sample,mut_more1_per_sample)
  df<-rbind(df,row)
}
# df<-data.frame(df)
colnames(df)<-c('gene','mut_num','sample_num',
                'high_num','moderate_num','low_num','modifier_num',
                'mut_1_per_sample','mut_more1_per_sample')
count_df<-df

df<-df %>% tibble::remove_rownames() %>% tibble::column_to_rownames("gene") 
library(purrr)
data<-map_df(df, as.numeric)%>%as.matrix()
rownames(data)<-rownames(df)
plot<-data %>% data.frame()%>% arrange(desc(mut_num))%>% arrange(desc(sample_num))

bk <- c(seq(0,max(plot[,3:6]),by=1))
pdf("out/Fig.2A-1 Mutation number.pdf",width = 4,height = 10)
pheatmap::pheatmap(plot[1:6],
                   cluster_rows = F,cluster_cols = F,scale = "none",
                   display_numbers = T,number_format = "%.0f",fontsize_number = 10,number_color = "black",
                   # color = "grey",border_color = "white",number_color = "white",
                   
                   color = c(colorRampPalette(colors = c("white","red"))(length(bk))),
                   breaks = bk,
                   gaps_col = c(2),
                   labels_col =c("Mutations","Samples",order)
                   )
dev.off()


df<-plot %>% dplyr::select(mut_1_per_sample,mut_more1_per_sample) %>% tibble::rownames_to_column(.,var = "gene")
colnames(df)<-c("gene","1 mut/sample","≥2 mut/sample")
df<-tidyr::pivot_longer(df, 
                        cols = c(`1 mut/sample`, `≥2 mut/sample`), 
                        names_to = "%Sample", 
                        values_to = "Value")
df$gene<-factor(df$gene,levels = rev(rownames(plot)))

ggplot(df,aes(x=Value,y=gene,fill=`%Sample`))+
  geom_bar(stat = "identity", position = "fill",width = 0.9, color = "white", size = 0.8)+
  theme_minimal()+
  scale_fill_manual(values = c("#3D1F64","#7F7FB8"))+
  scale_x_continuous( breaks = c(0, 0.5, 1)) 
ggsave(filename = "out/Fig.2A-2 sample ratio.pdf",width = 4,height = 10)


### 增加详细突变类型 20250506 ----

vcf<-merge.vcf %>% filter(Gene_Name %in% gene.use)
annotation_type<-vcf$Annotation %>% unique()
mapp<-read.table("further_analysis/CategoryMapping.txt",sep = "\t",header = T)
vcf<-merge(vcf,mapp,all.x=T,by.x="Annotation",by.y="VEP")

df_cate<-vcf %>% dplyr::select(Gene_Name,ID,Category,Annotation,Annotation_Impact) %>% distinct(ID,.keep_all = T)
df_cate$Gene_Name<-factor(df_cate$Gene_Name,levels = c(rownames(plot)))

uni_ID<-vcf$ID %>% unique()
uni_gene<-vcf$Gene_Name %>% unique()
type<-vcf$Annotation_Impact %>% unique()
order<-c("HIGH","MODERATE","LOW","MODIFIER")

gene_use<-rownames(plot) #[1:10] #部分基因
gene_use<-uni_gene    #所有基因

df<-data.frame()
vcf[vcf$Annotation %in% annotation_type[grep('&',annotation_type)], ]$Annotation<-"complex_type"
annotation_type<-vcf$Annotation %>% unique()
for (gene in gene_use) { 
  mini_df<-vcf %>% filter(Gene_Name==gene)
  
  mini_df$Annotation_Impact<-factor(mini_df$Annotation_Impact, levels = order, ordered = TRUE)
  mini_df<- mini_df[order(mini_df$Annotation_Impact), ]
  
  mut_num<-mini_df$ID %>% unique() %>% length()
  sample_num<-mini_df$SAMPLE %>% unique() %>% length()
  
  d<-mini_df %>% dplyr::select(ID,Gene_Name,Annotation) %>% distinct(ID,.keep_all = T)
  anno<-d$Annotation %>% unique()
  
  for (a in anno) {
    num<-d %>% filter(Annotation==a) %>% nrow()
    df<-rbind(df,c(gene,a,num))
  }
  
}
# df<-data.frame(df)
colnames(df)<-c('gene','Annotation','Num')

mut_df<-df
wide_df <- pivot_wider(df, 
                       id_cols = "gene", 
                       names_from = "Annotation", 
                       values_from = "Num") %>% tibble::column_to_rownames("gene")

out<-cbind(count_df,wide_df[count_df$gene,])
write.xlsx(out,file = "Gene.mut.stat&detail.xlsx")


df<-mut_df %>% filter(gene %in% rownames(plot))
df$gene<-factor(df$gene,levels = rev(rownames(plot)))
df$Annotation<-factor(df$Annotation,levels = c(
  "missense_variant"    , 
  "frameshift_variant",   
  "conservative_inframe_insertion"   ,  
  "3_prime_UTR_variant" ,
  "5_prime_UTR_variant",
  "upstream_gene_variant" ,
  "downstream_gene_variant",
  "stop_gained"   ,                          
  "stop_lost"                        ,        
  "intron_variant"   ,    
  "non_coding_transcript_exon_variant"  , 
  "synonymous_variant"  , 
  'complex_type'
))

df$Num<-as.numeric(df$Num)
ggplot(df,aes(x=gene,y=Num,fill=Annotation))+
  geom_bar(stat = 'identity',width = 0.8,position = "stack")+
  scale_fill_manual(values = c("#86BEEC", "#CE403F", "#BEBADA" ,"#FB8072", "#80B1A3", "#FDB462", "#B3DE69" ,"#FCCDE5", "#D9D9D9" ,"#BC80BD", "#CCEBC5","lightblue")
                                                                                 )+
  coord_flip()+
  xlab("")+ylab('Frequency')+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12))

ggsave("out/Fig.2A-3 Frequency of mut-type.pdf",width = 6,height = 12)


#########################################################################################
#### Figure 2C 
#########################################################################################
setwd("~/Project/YuYong_TALL/zhengjie")
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(seriation)
library(tidyverse)
CD8<-read.xlsx("../rna-analysis/analysis/L-N.CD8.DESeq2.withTPM.xlsx")
expr.cd8<-CD8  %>% distinct(gene,.keep_all = T) %>% tibble::column_to_rownames("gene")

expr.cd8<-log2(expr.cd8[,c(8:16)]+1)

expr<-expr.cd8

annocol<-data.frame(name=colnames(expr),
                    Condition=lapply(colnames(expr),function(x){str_sub(x,1,1)[[1]]}) %>%as.character()
                    # Celltype=lapply(colnames(expr),function(x){strsplit(x,split = "_")[[1]][2]}) %>%as.character()
) %>% tibble::column_to_rownames(colnames(.)[1])

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
pdf("out/Fig.2B Heatmap.pdf",width = 6.5,height = 8)
draw(ht)
dev.off()


#########################################################################################
#### Figure 2D 
#########################################################################################
CD8<-read.xlsx("../rna-analysis/analysis/L-N.CD8.DESeq2.withTPM.xlsx")

df<-CD8 %>% dplyr::select(gene,log2FoldChange) %>% distinct(gene,.keep_all = T) %>% 
  arrange(desc(log2FoldChange)) %>% na.omit()
geneList <- setNames(df$log2FoldChange, df$gene)  


library(clusterProfiler)
library(msigdbr)
library(enrichplot)

path<-msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>% dplyr::select(gs_name,gene_symbol)
trans<-path$gene_symbol %>% unique() %>% human2mouse()

path<-merge(path,trans,all.x=T,by.x="gene_symbol",by.y="humanGene")

path<-path %>% dplyr::select(gs_name,mouseGene) %>% na.omit()


res<-clusterProfiler::GSEA(geneList,TERM2GENE = path,pvalueCutoff=1)

gsea<-res@result
write.xlsx(res@result,file = "../rna-analysis/analysis/L-N.CD8.GSEA-KEGG.xlsx")

path_id<-"KEGG_CELL_CYCLE"
plot_gsea<-function(path_id){
  
  title<-paste0(path_id,"\n",
                "NES = ",res@result %>% filter(ID==path_id)%>% pull(NES) %>% signif(.,2),"\t",
                "p-value = ",res@result %>% filter(ID==path_id)%>% pull(pvalue) %>% signif(.,2),"\t",
                "adjusted p-value = ",res@result %>% filter(ID==path_id)%>% pull(p.adjust) %>% signif(.,2)
                )
  gseaplot2(x = res,geneSetID = path_id,title = title,pvalue_table = F,color = "darkred",subplots = 1:2)
  ggsave(filename = paste0("out/Fig.2D ",path_id,".pdf"),width = 5.5,height = 3.2)
}

lapply(c(
'KEGG_CELL_CYCLE', #1
'KEGG_CHEMOKINE_SIGNALING_PATHWAY', #3
'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY'),FUN = plot_gsea)

## HALLMARK ----
path_MH<-msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name,gene_symbol)
res2<-clusterProfiler::GSEA(geneList,TERM2GENE = path_MH,pvalueCutoff=1)

gsea2<-res2@result
write.xlsx(res2@result,file = "../rna-analysis/analysis/L-N.CD8.GSEA-HALLMARK.xlsx")

plot_gsea2<-function(path_id){
  
  title<-paste0(path_id,"\n",
                "NES = ",res2@result %>% filter(ID==path_id)%>% pull(NES) %>% signif(.,2),"\t",
                "p-value = ",res2@result %>% filter(ID==path_id)%>% pull(pvalue) %>% signif(.,2),"\t",
                "adjusted p-value = ",res2@result %>% filter(ID==path_id)%>% pull(p.adjust) %>% signif(.,2)
  )
  gseaplot2(x = res2,geneSetID = path_id,title = title,pvalue_table = F,color = "darkred",subplots = 1:2)
  ggsave(filename = paste0("out/Fig.2D ",path_id,".pdf"),width = 5.5,height = 3.2)
}
lapply("HALLMARK_MTORC1_SIGNALING",FUN = plot_gsea2)

#########################################################################################
#### Figure 2E 
#########################################################################################

# 
### sc analysis to clone.inte.ccr.rds
source("script/funcs.R")
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(homologene)

DN_1<-process(samplename ="20221123_DN_RNA", project ="DN_1",onlyread = TRUE )
DN_2<-process(samplename ="20221127_DN_RNA", project ="DN_2",onlyread = TRUE )
Thy_1<-process(samplename ="20221123_Thy_RNA", project ="Thy_1",onlyread = TRUE )
Thy_2<-process(samplename ="20221127_Thy_RNA", project ="Thy_2",onlyread = TRUE )
DN_5<-process(samplename ="20240126_DN_RNA", project ="DN_5" ,onlyread = TRUE )
Thy_5<-process(samplename ="20240126_Thy_RNA", project ="Thy_5",onlyread = TRUE )

data<-merge(DN_1,c(DN_2,Thy_1,Thy_2,DN_5,Thy_5))
p<-VlnPlot(data,features = c('nCount_RNA','nFeature_RNA','percent.mt'),pt.size = 0,group.by = "orig.ident")
# ggsave(filename = "unfilter.vlnplot.pdf",width = 10,height = 4)

data.list<-list(DN1=DN_1,
                DN2=DN_2,
                Thy1=Thy_1,
                Thy2=Thy_2,
                DN5=DN_5,
                Thy5=Thy_5)
saveRDS(data.list,"../data.list.filtered&removeDB.rds")

seu <- merge(data.list$DN1,c(data.list$DN2,data.list$Thy1,data.list$Thy2))

DefaultAssay(seu)<-"RNA"
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
seu <- RunUMAP(seu, dims = 1:20, reduction = "integrated.cca")
seu <- FindNeighbors(seu,dims = 1:20,reduction = "integrated.cca")
seu <- FindClusters(seu,resolution = 0.8)
DimPlot(seu)
saveRDS(seu,"TALL.ccregress.integrate.rds")


### Use TCR to select tall cell
data<-readRDS("../SCIN.ccregress.integrate.rds")
data$barcode<-rownames(data@meta.data)


data.tcr<-readRDS("../data.tcr.rds")
DN_1<-data.tcr$DN_1
DN_2<-data.tcr$DN_2
Thy_1<-data.tcr$Thy_1
Thy_2<-data.tcr$Thy_2
DN_1$barcode<-paste0(DN_1$barcode,"_1")
DN_1$sample<-"DN_1"
DN_2$barcode<-paste0(DN_2$barcode,"_2")
DN_2$sample<-"DN_2"
Thy_1$barcode<-paste0(Thy_1$barcode,"_3")
Thy_1$sample<-"Thy_1"
Thy_2$barcode<-paste0(Thy_2$barcode,"_4")
Thy_2$sample<-"Thy_2"
# 所有TCR 数据
TCR<-rbind(DN_1,DN_2,Thy_1,Thy_2)

overlap<-intersect(rownames(data@meta.data),TCR$barcode)
TCR.use<-TCR%>%filter(barcode %in% overlap)

data$withTCR<-"none"
data@meta.data[data@meta.data$barcode%in%overlap,]$withTCR<-"yes"
DimPlot(data,group.by = "withTCR")+ggtitle(label = paste0("cell with TCR(",round(length(overlap)/nrow(data@meta.data),4)*100,"%)"))
# ggsave(filename = "./out/data.4.withTCR.pdf",width = 5,height = 4.6)

data.use<-subset(data,withTCR=="yes")
overlap<-intersect(rownames(data.use@meta.data),TCR$barcode)
TCR.use<-TCR%>%filter(barcode %in% overlap)


freq.tcr<-TCR.use$cdr3%>%table()%>%data.frame()
cdr3_1<-'CASSDAAGARDTQYF' #DN1 THY1
cdr3_2<-'CAWSLTGGARSQNTLYF' #DN2 THY2
bc1<-TCR.use[TCR.use$cdr3==cdr3_1,]$barcode%>%unique()
bc2<-TCR.use[TCR.use$cdr3==cdr3_2,]$barcode%>%unique()

data.use$cdr3_clone<-"other"
data.use@meta.data[data.use@meta.data$barcode%in%bc1,]$cdr3_clone<-"clone1"
data.use@meta.data[data.use@meta.data$barcode%in%bc2,]$cdr3_clone<-"clone2"
DimPlot(data.use,group.by = "cdr3_clone",split.by = "orig.ident")+ggtitle(label ="cdr3 clone",
                                                                          subtitle = paste0("clone1(",round(length(bc1)/nrow(data.use@meta.data),4)*100,"%): ",cdr3_1,
                                                                                            " (DN1/Thy1)\nclone2(",round(length(bc2)/nrow(data.use@meta.data),4)*100,"%): ",cdr3_2," (DN2/Thy2)"))

table(data.use$cdr3_clone,data.use$orig.ident) %>% heatmap()

seu <- RunUMAP(data, dims = 1:20, reduction = "integrated.cca")
seu <- FindNeighbors(seu,dims = 1:20,reduction = "integrated.cca")
seu <- FindClusters(seu,resolution = 0.8)
p<-FeaturePlot(seu,features = c("Ptprc","Cd3d","Kit","Cd4",
                             "Cd8a","Il2ra","Notch1","Bcl11a",
                             "Bcl11b","Runx1","Runx3","Mki67"),cols = c("grey","red"))
# ggsave(plot=p,filename = "Clone.features.pdf",width = 14,height = 9)

## celltype annotation ---
library(celldex)
library(SingleR)
seu<-annoBySingler(seu)
seu@active.ident<-seu$seurat_clusters

## check cell marker 
# feat<-c(
#   # "Itgam","Itgax", #Myeloid
#   # "Ms4a1","Cd79a",'Iglc2',#'Iglc3',#'Ctsh', #B cell
#   "Kit","Il2ra",
#   'Hes1','Mpzl2', #DN
#   # 'Tnfrsf4','Tnfrsf9','Foxp3',# T reg
#   "Cd3d","Ptprc","Lef1",
#   "Cd4","Cd8a", #DP
#   "Mki67", # profi
#   "Bcl11a","Bcl11b",
#   "Notch1", # T-ALL
#   'Tcrg-C2','Cpa3','Trdc'# Tgd
# )
# DotPlot(seu,cluster.idents = T,features =  feat,assay = "RNA")+coord_flip()

## rename ident
seu<-RenameIdents(seu,
                  "11"="DN2 (11)",  #Kit高, Il2ra高，略带Cd8a
                  "12"="DN4 (12)",  #Kit低, Il2ra低，略带Cd8a
                  "3"="DP (3)",
                  "0"="CD8+ (other)",
                  "1"="CD8+ (other)",
                  "2"="CD8+ (other)",
                  "4"="CD8+ (other)",
                  "5"="CD8+ (other)",
                  "6"="CD8+ (other)",
                  "7"="CD8+ (other)",
                  "8"="CD8+ (other)",
                  "9"="CD8+ (other)",
                  "10"="CD8+ (other)"
                  )
saveRDS(seu,"clone.inte.ccr.rds")

### TCR data
DN_1_contig<-read.csv("cellranger/20221123_DN_TCR/outs/filtered_contig_annotations.csv")
DN_2_contig<-read.csv("cellranger/20221127_DN_TCR/outs/filtered_contig_annotations.csv")
Thy_1_contig<-read.csv("cellranger/20221123_Thy_TCR/outs/filtered_contig_annotations.csv")
Thy_2_contig<-read.csv("cellranger/20221127_Thy_TCR/outs/filtered_contig_annotations.csv")
data.tcr<-list(DN_1=DN_1_contig,DN_2=DN_2_contig,Thy_1=Thy_1_contig,Thy_2=Thy_2_contig)
saveRDS(data.tcr,"data.tcr.rds")

### run slingshot
library(slingshot)
library(RColorBrewer)
library(paletteer) 
library(viridis)
library(Seurat)
library(tradeSeq)
library(tidyverse)
library(grDevices)
library(ggthemes)
library(Seurat)
data<-readRDS("clone.inte.ccr.rds")


sce <- as.SingleCellExperiment(data,assay = "RNA")
sce <- slingshot(sce, clusterLabels = 'ident', reducedDim = 'UMAP', start.clus= "DN2 (11)",end.clus=c("CD8+ (other)"))
saveRDS(sce,"clone.slingshot.rds")

#设置颜色
cell_pal <- function(cell_vars, colors,...) {
  categories <- sort(unique(cell_vars))
  pal <- setNames(colors, categories)
  return(pal[cell_vars])
}

cell_colors <- cell_pal(sce$ident, scales::hue_pal()(9))
pdf(file = "out/Fig.2E clone.slngshot.pdf",width = 5,height = 5)
plot(reducedDims(sce)$UMAP, col = cell_colors,pch=16, asp = 1, cex = 0.4)
lines(SlingshotDataSet(sce), lwd=2, col='black',type = "lineages")

#计算celltype坐标位置，用于图中标记
celltype_label <- data@reductions$umap@cell.embeddings%>% as.data.frame() %>%cbind(celltype = data@meta.data$ident) %>%
  group_by(celltype) %>%summarise(UMAP1 = median(umap_1),UMAP2 = median(umap_2))
index<-celltype_label
for (i in 1:nrow(index)) {text(index$celltype[i], x=index$UMAP1[i]-1, y=index$UMAP2[i])}
dev.off()



#########################################################################################
#### Figure 2F
#########################################################################################


setwd("~/Project/YuYong_TALL/out/SCENIC_mergeCTL_Ellen")

## 1.extract out_SCENIC.loom data
loom <- open_loom('aucell_output.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

## 2. load data
{
  seu<-readRDS("clone.inte.ccr.rds")
  DimPlot(seu)
  ctl<-readRDS("data.ctl.assigned.rds")  # which is finished in utils/Normal-lineage.R
  DimPlot(ctl,group.by = "orig.ident")
  ctl_use<-subset(ctl,orig.ident %in%c('DN','THY'))
  DimPlot(ctl_use,label = T)
  DimPlot(ctl_use,group.by = "orig.ident",label = T)
  main_cell<-table(ctl_use$ident) %>% sort() %>% data.frame() %>% filter(Freq>50) %>% pull(Var1) %>% as.character()

  ctl_use_main<-subset(ctl_use,ident %in% main_cell)
  DotPlot(ctl_use_main,features = c('Ptprc',"Cd3d","Cd3g","Cd3e","Cd4",
                                    "Cd8a","Ms4a1","Cd68","Cd34","Kit","Hes1","Il2ra"),cluster.idents = T)

  normal<-ctl_use_main
  normal<-RunUMAP(normal,dims = 1:30)
  DimPlot(normal,label = T,group.by = "orig.ident")
  normal$orig.ident<-paste0("N_",as.character(normal$orig.ident))

  seu@assays$SCT<-NULL
  DefaultAssay(seu)<-"RNA"
  normal@assays$SCT<-NULL
  DefaultAssay(normal)<-"RNA"
  seu@reductions
  normal@reductions

  data<-merge(seu,normal)
  data<-JoinLayers(data,assay = "RNA",layers = "counts")
  saveRDS(data,"merge_ctl&clone.rds")
} # code to get merge_ctl&clone.rds
data<-readRDS("merge_ctl&clone.rds")
{
  data<-readRDS("merge_ctl&clone.rds")

  DimPlot(data)
  normal<-subset(data,orig.ident %in% c("N_DN","N_THY"))
  seu<-NormalizeData(normal)
  seu<-FindVariableFeatures(seu,selection.method = "vst")
  seu<-ScaleData(seu)
  seu<-RunPCA(seu)
  seu<-FindNeighbors(seu,dims = 1:20)
  seu<-RunUMAP(seu,dims = 1:20)
  seu<-FindClusters(seu,resolution = 1.0)
  seu@active.ident

  DotPlot(seu,features = c("Ptprc","Cd3d","Cd3g","Cd3e","Cd4","Cd8a","Cd8b1",
                          "Kit",#"Hlf","Hoxa9",
                          "Il2ra","Cd44","Cd24a",
                          "Lef1","Bcl11a","Bcl11b","Notch1",
                          "Nkg7","Klra1","Klrb1a","Ncam1"),
          group.by = "seurat_clusters",cluster.idents = T)+
    theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0))


  seu@active.ident<-seu$seurat_clusters

  seu<-RenameIdents(seu,
                    "16"="ETP",
                    "10"="DN2",
                    "6"="DN2",
                    "21"="DN2",
                    
                    "20"="DN3",
                    "15"="DN3",
                    "1"="DN3",
                    "3"="DN3",
                    "18"="DN3",
                    "17"="DN3",
                    
                    
                    "13"="DN4",
                    "9"="DN4",
                    
                    "12"="DP",
                    "14"="DP",
                    "19"="DP",
                    "4"="DP",
                    "22"="DP",
                    "11"="DP",
                    "0"="DP",
                    "2"="DP",
                    # "9"="DP",
                    
                    "8"="CD8+",
                    "5"="CD8+",
                    "7"="CD8+"
                    )

  seu$celltype_normal<-seu@active.ident
  normal<-seu

  DotPlot(normal,features = c("Ptprc","Cd3d",#"Cd3g","Cd3e",
                          "Kit",#"Hlf","Hoxa9",
                          
                          "Il2ra","Bcl11a","Cd28",
                          # "Cd44",
                          "Cd24a",
                          "Cd4","Cd8a","Cd8b1",
                          
                          "Lef1","Bcl11b","Notch1"),cols = c("lightblue","red"),
        cluster.idents = F)+
    coord_flip()+xlab("")+ylab("")+theme_bw()+
    theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0,size = 12,color="black"),
          axis.text.y = element_text(size = 12,color="black")
          )

  # ggsave("Normal.celltype.marker.2.pdf",width = 5,height = 4.5)

  saveRDS(seu,"normal.scenic.RDS")
} # code to get normal.scenic.RDS, which is to run SCENIC in normal samples(DN & Thy)
normal<-readRDS("normal.scenic.RDS")
Ellen<-readRDS("ref_data/GSE183026_ScienceImmunology2022/out/for_SCIN.rds") # # which is finished in utils/Ellen-data.R


data$celltype<-as.character(data$ident)
data@meta.data[normal@meta.data %>% rownames(),]$celltype<-paste0("N-",normal@active.ident %>% as.character())

data$celltype %>% table()

Thy<-subset(data2,orig.ident=="Thy")
Thy$celltype<-paste0("Ellen-",Thy@active.ident %>% as.character())
Thy@meta.data[Thy@meta.data$orig.ident=="Thy",]$orig.ident<-"Ellen_Thy"
Thy$celltype %>% table()

seu<-merge(data,Thy)

seu$orig.ident %>% table()
seu$celltype %>% table()

seu$celltype



## check expression --
# seurat.data<-seu
# # DimPlot(seu)
# seurat.data<-JoinLayers(seurat.data,assay = "RNA")
# seurat.data<-NormalizeData(seurat.data)
# 
# p1<-DotPlot(Thy,features = c("Ptprc","Kit","Bcl11a","Il2ra","Bcl11b","Lef1","Cd3g","Cd4","Cd8a"))+
#   coord_flip()+xlab("")+ylab("")+theme_bw()+
#   theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0,size = 12,color="black"),
#         axis.text.y = element_text(size = 12,color="black")
#   )+ggtitle("Ellen-Thy")


seu@meta.data[seu@meta.data$celltype=="Ellen-ETP",]$celltype<-"N-ETP"
seu@meta.data[seu@meta.data$celltype=="Ellen-DN2",]$celltype<-"N-DN2"
seu@meta.data[seu@meta.data$celltype=="Ellen-DN3",]$celltype<-"N-DN3"

seu$celltype<-factor(seu$celltype,levels = c("N-ETP","N-DN2" ,"N-DN3" ,"N-DN4","N-DP" , "N-CD8+" ,     
                                             "DN2 (11)","DN4 (12)","DP (3)","CD8+ (other)"))

seu<-JoinLayers(seu,assay = "RNA")
seu<-NormalizeData(seu)



DotPlot(seu,features = c("Ptprc","Cd3d",#"Cd3g","Cd3e",
                            "Kit",#"Hlf","Hoxa9",

                            "Il2ra","Bcl11a","Cd28",
                            # "Cd44",
                            "Cd24a",
                            "Cd4","Cd8a","Cd8b1",

                            "Lef1","Bcl11b","Notch1"),cols = c("lightblue","red"),
        cluster.idents = F,group.by = "celltype")+
  coord_flip()+xlab("")+ylab("")+theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0,size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black")
  )

saveRDS(seu,"scenic_merge3data.rds")
seurat.data<-seu

## 3. view

sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
# dim(sub_regulonAUC)
# seurat.data
# identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$celltype)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

# save
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')
load("for_rss_and_visual.Rdata")

## plot ----
# TF activity mean

selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

# remove extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
out<-data.frame(regulonActivity_byGroup_Scaled,check.names = F,check.rows = F) %>% tibble::rownames_to_column("TF")

# write.xlsx(out,file = "regulonActivity.byGroup.Scaled.0609.xlsx")

dim(regulonActivity_byGroup_Scaled)
#[1] 209   9
plot=na.omit(regulonActivity_byGroup_Scaled)
Heatmap(
  plot,
  name                         = "z-score",
  
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  # row_title_rot                = 90,
  column_names_rot = -90,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


### 4.2. rss查看特异TF
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)

rss_high_tf<-rssPlot$df %>% group_by(cellType) %>% arrange(desc(RSS)) %>% slice_head(n = 5) %>% pull(Topic) %>% as.character() %>% unique()

p1<-plotly::ggplotly(rssPlot$plot) 
# ggsave(filename = "SCENIC.TF.RSS.celltype-dotplot.0609.pdf",width = 3.5,height = 15)

rss_high_tf
rssPlot <- plotRSS(rss[rss_high_tf,])
p2<-plotly::ggplotly(rssPlot$plot)
# ggsave(filename = "SCENIC.TF.RSS.top10.celltype-dotplot.0609.pdf",width = 3.5,height = 6)


score=regulonActivity_byGroup_Scaled
head(score)
df = do.call(rbind,
             lapply(1:ncol(score), function(i){
               dat= data.frame(
                 path  = rownames(score),
                 cluster =   colnames(score)[i],
                 sd.1 = score[,i],
                 sd.2 = apply(score[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)

n = score[c(top5$path) %>% unique(),] 

## show selected TFs
regulonsToPlot = c(
  "Twist1(+)","Meis1(+)","Spi1(+)","Gfi1b(+)","Erg(+)",'Tcf7l2(+)',"Nfe2(+)","Hmga2(+)","Jun(+)","Pbx1(+)","Tfap4(+)",
  "Otx2(+)","Myc(+)","Runx1(+)",
  "Stat1(+)","Relb(+)","Maf(+)",'Runx3(+)',
  "Ets2(+)","Bcl6(+)","Klf13(+)",
  "E2f1(+)","Tfdp1(+)","Hdac2(+)",
  'E2f8(+)','Vdr(+)',
  
  'Tal1(+)',"Runx1(+)",'Bcl11b(+)','Myb(+)',"Myc(+)","Dnmt1(+)","E2f7(+)","Tcf12(+)","Tcf3(+)","Ikzf1(+)",'Tcf7l1(+)','Hoxa10(+)',"Ets1(+)",
  "Smad4(+)","Tef(+)","Gata2(+)",
  "Emx2(+)","Ikzf3(+)","Dpf1(+)","Cebpa(+)",
  "Atf5(+)","Tbx21(+)"
                   ) %>% unique()
regulonsToPlot = regulonsToPlot[regulonsToPlot %in% row.names(sub_regulonAUC)]

top2 <- df %>% group_by(cluster) %>% top_n(3, fc)
rank <- df %>% group_by(cluster) %>% top_n(100, fc)
n = score[c(regulonsToPlot) %>% unique(),] 

genes <- unique(c(regulonsToPlot, top2$path))
genes <- unique(c(regulonsToPlot))
existing_genes <- genes[genes %in% rownames(n)]

at_pos <- which(rownames(n) %in% existing_genes)
labels_vec <- rownames(n)[at_pos]  # 按热图实际行顺序

pheatmap::pheatmap(n,cluster_rows = F,cluster_cols = F)

# plot heatmap
ht <- Heatmap(
  n,
  name = "Z-score",
  cluster_rows = F,
  cluster_columns = F,
  clustering_method_columns = "ward.D",
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("#1A4d73", "white", "#c12c1f")),
  row_gap = unit(1.5, "mm"), 
  column_gap = unit(1.5, "mm"), 
  border = TRUE,
  top_annotation = columnAnnotation(df = data.frame(Group=c(rep("Normal",6),rep("T-ALL",4))),
                                    col=list(Group=c("T-ALL"="#F69389","Normal"="lightblue"))),
)
  

ht_opt$message = FALSE
draw(ht)

pdf("out/Fig.2F TF.SCENIC withCTL.pdf",width = 5.2,height = 8)
draw(ht)
dev.off()



