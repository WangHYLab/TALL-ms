#########################################################################################
#### Fig.S2A
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
write.xlsx(merge.vcf,file = "data/merge.filtered.vcf")

merge.vcf<-read.xlsx("data/merge.filtered.vcf.xlsx")

## 挑选driver gene mutation landscape
driver_gene<-openxlsx::read.xlsx("data/28671688-supp.xlsx",sheet = "Table S9 Driver Mutations")

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
pdf("out/Fig.S2A-1 Mutation number.pdf",width = 4,height = 10)
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
ggsave(filename = "out/Fig.S2A-2 sample ratio.pdf",width = 4,height = 10)


### 增加详细突变类型 20250506 ----

vcf<-merge.vcf %>% filter(Gene_Name %in% gene.use)
annotation_type<-vcf$Annotation %>% unique()
mapp<-read.table("data/CategoryMapping.txt",sep = "\t",header = T)
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
write.xlsx(out,file = "data/Gene.mut.stat&detail.xlsx")


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

ggsave("out/Fig.S2A-3 Frequency of mut-type.pdf",width = 6,height = 12)



#########################################################################################
#### Fig.S2B
#########################################################################################


library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(vcfR)
library(openxlsx)
setwd("WES-analysis")


### read vcf and merge
files <- list.files("out/anno", pattern = ".vcf")
filter <- TRUE
merge.vcf <- c()
for (file in files) {
  out_vcf <- c()
  vcf <- read.vcfR(paste0("out/anno/", file), verbose = T)
  tidy_vcf <- vcfR2tidy(vcf)

  VCF <- tidy_vcf$fix
  if (filter) {
    VCF <- tidy_vcf$fix %>% filter(FILTER == "PASS")
  }

  VCF <- VCF %>% select(CHROM, POS, ID, REF, ALT, FILTER, ANN)
  VCF$ID <- paste0(VCF$CHROM, "_", VCF$POS, "_", VCF$REF, "_", VCF$ALT)
  VCF$SAMPLE <- file
  pb <- txtProgressBar(style = 3)
  for (line in 1:nrow(VCF)) {
    setTxtProgressBar(pb, line / nrow(VCF), title = paste0("Process sample: ", file))
    ANNO_ori <- VCF[line, ]$ANN
    # anno:
    # Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
    anno_split <- strsplit(ANNO_ori, ",")[[1]]
    split_strings <- lapply(anno_split, function(x) {
      strsplit(x, split = "\\|")[[1]]
    })
    # max_cols <- max(sapply(split_strings, length))
    max_cols <- 16
    anno_mtx <- do.call(rbind, lapply(split_strings, function(x) c(x, rep(NA, max_cols - length(x)))))
    colnames(anno_mtx) <- c(
      "Allele",
      "Annotation",
      "Annotation_Impact",
      "Gene_Name",
      "Gene_ID",
      "Feature_Type",
      "Feature_ID",
      "Transcript_BioType",
      "Rank",
      "HGVS.c",
      "HGVS.p",
      "cDNA.pos / cDNA.length",
      "CDS.pos / CDS.length",
      "AA.pos / AA.length",
      "Distance",
      "ERRORS / WARNINGS / INFO"
    )
    rep_len <- nrow(anno_mtx)
    vcf <- c()
    for (i in 1:rep_len) {
      vcf <- rbind(vcf, VCF[line, ])
    }
    out_vcf <- rbind(out_vcf, cbind(vcf, anno_mtx))
  }
  close(pb)
  merge.vcf <- rbind(merge.vcf, out_vcf)
}
write.xlsx(merge.vcf, file = "data/merge.filtered.vcf.xlsx")

### mut count file extract

merge.vcf <- read.xlsx("data/merge.filtered.vcf.xlsx")

vcf <- merge.vcf %>% filter(Gene_Name %in% gene.use)
annotation_type <- vcf$Annotation %>% unique()
mapp <- read.table("data/CategoryMapping.txt", sep = "\t", header = T)
vcf <- merge(vcf, mapp, all.x = T, by.x = "Annotation", by.y = "VEP")

df_cate <- vcf %>%
  dplyr::select(Gene_Name, ID, Category, Annotation, Annotation_Impact) %>%
  distinct(ID, .keep_all = T)
df_cate$Gene_Name <- factor(df_cate$Gene_Name, levels = c(rownames(plot)))

uni_ID <- vcf$ID %>% unique()
uni_gene <- vcf$Gene_Name %>% unique()
type <- vcf$Annotation_Impact %>% unique()
order <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

gene_use <- rownames(plot) # [1:10] #部分基因
gene_use <- uni_gene # 所有基因

df <- data.frame()
vcf[vcf$Annotation %in% annotation_type[grep("&", annotation_type)], ]$Annotation <- "complex_type"
annotation_type <- vcf$Annotation %>% unique()
for (gene in gene_use) {
  mini_df <- vcf %>% filter(Gene_Name == gene)

  mini_df$Annotation_Impact <- factor(mini_df$Annotation_Impact, levels = order, ordered = TRUE)
  mini_df <- mini_df[order(mini_df$Annotation_Impact), ]

  mut_num <- mini_df$ID %>%
    unique() %>%
    length()
  sample_num <- mini_df$SAMPLE %>%
    unique() %>%
    length()

  d <- mini_df %>%
    dplyr::select(ID, Gene_Name, Annotation) %>%
    distinct(ID, .keep_all = T)
  anno <- d$Annotation %>% unique()

  for (a in anno) {
    num <- d %>%
      filter(Annotation == a) %>%
      nrow()
    df <- rbind(df, c(gene, a, num))
  }
}
# df<-data.frame(df)
colnames(df) <- c("gene", "Annotation", "Num")

mut_df <- df
wide_df <- pivot_wider(df,
  id_cols = "gene",
  names_from = "Annotation",
  values_from = "Num"
) %>% tibble::column_to_rownames("gene")

out <- cbind(count_df, wide_df[count_df$gene, ])
write.xlsx(out, file = "data/Gene.mut.stat&detail.xlsx")

### plot 
path_MH<-msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name,gene_symbol)
path_GO<-msigdbr(species = "Mus musculus", category = "C5","BP") %>% dplyr::select(gs_name,gene_symbol)
path_BIOCARTA<-msigdbr(species = "Mus musculus", category = "C2","BIOCARTA") %>% dplyr::select(gs_name,gene_symbol)
path_set<-rbind(path_MH,path_GO,path_BIOCARTA)
stat<-read.xlsx("data/Gene.mut.stat&detail.xlsx",rowNames = T)

path <- c(
  "HALLMARK_NOTCH_SIGNALING",
  "GOBP_THYMOCYTE_APOPTOTIC_PROCESS",
  "GOBP_REGULATION_OF_LEUKOCYTE_APOPTOTIC_PROCESS",
  "GOBP_REGULATION_OF_TRANSCRIPTION_REGULATORY_REGION_DNA_BINDING",
  "GOBP_REGULATION_OF_HISTONE_MODIFICATION",
  "GOBP_NEGATIVE_REGULATION_OF_GENE_EXPRESSION_EPIGENETIC"
)



plot_mut_path<-function(path_select){
  gene_in_path<-path_set %>% filter(gs_name == path_select) %>% pull(gene_symbol)
  
  df<-stat[gene_in_path,]
  
  ## top mut gene
  tmp<-apply(df[,1:2],2,as.numeric) %>% data.frame()
  rownames(tmp)<-rownames(df)
  gene_order<-tmp %>% arrange(desc(mut_num)) %>% rownames() %>% head(10)
  
  ## gene 2 plot
  df<-stat[gene_order,] %>% 
    arrange(desc(mut_num))
  
  plot_df<-apply(df,2,as.numeric) %>% data.frame()
  
  rownames(plot_df)<-rownames(df)
  # pheatmap::pheatmap(plot_df[,3:6],cluster_rows = F,cluster_cols = F,display_numbers = T)
  
  gg_df<-gather(plot_df[,3:6],key="class",value = "value")
  gg_df$gene<-rep(rownames(df),4)
  
  gg_df$gene<-factor(gg_df$gene,levels = gene_order)
  gg_df[gg_df$class=="high_num",]$class<-"HIGH"
  gg_df[gg_df$class=="moderate_num",]$class<-"MODERATE"
  gg_df[gg_df$class=="low_num",]$class<-"LOW"
  gg_df[gg_df$class=="modifier_num",]$class<-"MODIFIER"
  
  
  gg_df$class<-factor(gg_df$class,levels = c("HIGH","MODERATE","LOW","MODIFIER"))
  
  p<-ggplot(gg_df,aes(x=gene,y=value,fill=class))+
    geom_bar(stat="identity",position = "stack",width = 0.7)+
    scale_fill_manual(values = c("darkred","salmon","pink","grey"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0,size = 10),
          panel.grid = element_blank())+
    ylab("Altered number")+
    xlab("Top altered gene")+
    ggtitle(path_select)
  
  return(p)
}

res<-lapply(path,FUN = plot_mut_path)

ggsave(plot = (res[[1]]|res[[2]]|res[[3]])/(res[[4]]|res[[5]]|res[[6]])
,filename = "out/Fig.S2B mutFreq pathway gene.pdf",width = 15,height = 6
)



