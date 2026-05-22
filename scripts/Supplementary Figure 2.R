#########################################################################################
#### Fig.S2A
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
write.xlsx(merge.vcf, file = "further_analysis/merge.filtered.vcf.xlsx")

### mut count file extract

merge.vcf <- read.xlsx("further_analysis/merge.filtered.vcf.xlsx")

vcf <- merge.vcf %>% filter(Gene_Name %in% gene.use)
annotation_type <- vcf$Annotation %>% unique()
mapp <- read.table("further_analysis/CategoryMapping.txt", sep = "\t", header = T)
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
write.xlsx(out, file = "Gene.mut.stat&detail.xlsx")

### plot 
path_MH<-msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name,gene_symbol)
path_GO<-msigdbr(species = "Mus musculus", category = "C5","BP") %>% dplyr::select(gs_name,gene_symbol)
path_BIOCARTA<-msigdbr(species = "Mus musculus", category = "C2","BIOCARTA") %>% dplyr::select(gs_name,gene_symbol)
path_set<-rbind(path_MH,path_GO,path_BIOCARTA)
stat<-read.xlsx("Gene.mut.stat&detail.xlsx",rowNames = T)

path<-c("HALLMARK_NOTCH_SIGNALING",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_KRAS_SIGNALING_UP",
        "GOBP_T_CELL_DIFFERENTIATION",
        "GOBP_MYELOID_LEUKOCYTE_ACTIVATION")

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
,filename = "out/Fig.S2A mutFreq pathway gene.pdf",width = 15,height = 6
)


#########################################################################################
#### Fig.S2B
#########################################################################################

# Filtered Notch1 mutations from file Gene.mut.stat&detail.xlsx
# Mutation domains added manually