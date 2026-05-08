#########################################################################################
#### Figure 7F
#########################################################################################
setwd("ov-analysis/analysis/")

library(tidyverse)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(corrplot)
library(openxlsx)
sample.re<-read.xlsx("sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])

#### DESeq2 analysis
run_DESeq2<-function(key1,key2,out_postfix,tpm,raw_count,write=TRUE){
  # id2gene<-read.table("/home/zhengjie/DATA/ensmusg2symbol.GRCm39.txt",header = T,sep = "\t",row.names = 1)
  id2gene<-read.table("mm10.id2symbol.txt",header = T,sep = "\t",row.names = 1)
  id2description<-read.table("ensmusg2description.GRCm39.txt",header = T,sep = "\t",row.names = 1,check.names = F,quote = "")
  sample.re<-read.xlsx("sample.xlsx",sheet = "bulk")%>%tibble::column_to_rownames(names(.)[1])
  
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
  gmt.mh<-read.gmt('/home/zhengjie/DATA/msigdb/mh.all.v2023.2.Mm.symbols.gmt')
  egmt.mh <- GSEA(geneList, TERM2GENE=gmt.mh,eps = 0,seed = 1,nPermSimple=10000,pvalueCutoff = 1)
  out.mh<-egmt.mh@result
  message(paste0("Sig path: ",sum(out.mh$pvalue<0.05%>%as.numeric())))
  
  message("Run GSEA in Reactome...")
  gmt.reactome<-read.gmt('/home/zhengjie/DATA/msigdb/m2.cp.reactome.v2023.2.Mm.symbols.gmt')
  egmt.reactome <- GSEA(geneList, TERM2GENE=gmt.reactome,eps = 0,seed = 1,nPermSimple=10000,pvalueCutoff = 1)
  out.reactome<-egmt.reactome@result
  message(paste0("Sig path: ",sum(out.reactome$pvalue<0.05%>%as.numeric())))
  
  message("Run GSEA in GO...")
  gmt.go<-read.gmt('/home/zhengjie/DATA/msigdb/m5.go.v2023.2.Mm.symbols.gmt')
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
## run DESeq2 and GSEA
{
    raw_count<-read.table("../gene-count-matrix.txt",check.names = F)
    colname.new<-lapply(colnames(raw_count),function(x){strsplit(strsplit(x,"/home/zhengjie/Project/YuYong_TALL/ov-analysis/mapping/")[[1]][2],"_")[[1]][1]})%>%as.character()
    colnames(raw_count)<-colname.new
    colnames(raw_count)<-sample.re[colnames(raw_count),]

    # raw_count<-raw_count[rowSums(raw_count)!=0,]
    raw_count<-round(as.matrix(raw_count))

    data<-read.table("../gene-TPM-matrix.txt",check.names = F)
    colname.new<-lapply(colnames(data),function(x){strsplit(strsplit(x,"/home/zhengjie/Project/YuYong_TALL/ov-analysis/mapping/")[[1]][2],"_")[[1]][1]})%>%as.character()
    colnames(data)<-colname.new
    colnames(data)<-sample.re[colnames(data),]
    tpm<-data

    sample<-colnames(raw_count)

    # OV_EGFP vs CTL_EGFP----
    dis<-sample[grep("OV_EGFP[1-3]",sample)]
    ctl<-sample[grep("CTL_EGFP[1-3]",sample)]
    ov_ctl<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "EGFP.OV-CTL",tpm = tpm,raw_count = raw_count,write = TRUE)
    gsea.ov_ctl<-run_gsea(deseq2_res = ov_ctl,file_postfix = "EGFP.OV-CTL")

    # OV_td vs CTL_td----
    dis<-sample[grep("OV_TD[1-3]",sample)]
    ctl<-sample[grep("CTL_TD[1-3]",sample)]
    ov_ctl<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "TD.OV-CTL",tpm = tpm,raw_count = raw_count,write = TRUE)
    gsea.ov_ctl<-run_gsea(deseq2_res = ov_ctl,file_postfix = "TD.OV-CTL")


    # OV_TD vs OV_EGFP 20241120----
    sample<-colnames(raw_count)
    dis<-sample[grep("OV_EGFP[1-3]",sample)]
    ctl<-sample[grep("OV_TD[1-3]",sample)]
    EGFP_TD<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "OV.EGFP-TD",tpm = tpm,raw_count = raw_count,write = TRUE)
    gsea.EGFP_TD<-run_gsea(deseq2_res = EGFP_TD,file_postfix = "OV.EGFP-TD",merge_save = TRUE)


    # CTL EGFP-TD ----
    sample<-colnames(raw_count)
    dis<-sample[grep("CTL_EGFP[1-3]",sample)]
    ctl<-sample[grep("CTL_TD[1-3]",sample)]
    EGFP_TD<-run_DESeq2(key1 = dis,key2 = ctl,out_postfix = "CTL.EGFP-TD",tpm = tpm,raw_count = raw_count,write = TRUE)
    gsea.EGFP_TD<-run_gsea(deseq2_res = EGFP_TD,file_postfix = "CTL.EGFP-TD",merge_save = TRUE)

}

## heatmap
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

data<-read.xlsx("OV.EGFP-TD.DESeq2.withTPM.xlsx")
deg<-data %>% filter(abs(log2FoldChange)>1&padj<0.05) %>% arrange(desc(log2FoldChange)) %>% pull(gene)
df<-df[c(select_gene,deg) %>% unique(),c( "CTL_EGFP1", "CTL_EGFP2", "CTL_EGFP3",
              "CTL_TD1" ,  "CTL_TD2" ,  "CTL_TD3" ,
              
              "OV_TD1" ,   "OV_TD2"  ,  "OV_TD3",
              "OV_EGFP1",  "OV_EGFP2" , "OV_EGFP3" )]

pdf(file = "out/Fig.7f OV.EGFP-TD.deg.heatmap.withCTL.pdf",width = 4.5,height = 9)
bk <- c(seq(-2,2,by=0.01))
pheatmap::pheatmap(df,scale = "row",
                   color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),
                             colorRampPalette(colors = c("white","#d23918"))(length(bk)/2)),
                   cluster_cols = F,cluster_rows =T,
                   gaps_col = c(6,9),
                   cutree_rows = 2,
                   treeheight_row = 20,
                   # gaps_row = 20
                   )
dev.off()




#########################################################################################
#### Figure 7G
#########################################################################################

## heatmap
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

# data<-read.xlsx("TD.OV-CTL.DESeq2.withTPM.xlsx")
# deg1<-data %>% filter(log2FoldChange>1&padj<0.05) %>% arrange(desc(log2FoldChange)) %>% head(30)%>% pull(gene) 
# deg2<-data %>% filter(log2FoldChange< -1&padj<0.05) %>% arrange(desc(log2FoldChange)) %>% head(30)%>% pull(gene) 
# deg<-c(deg1,deg2)

select_gene<-read.xlsx("TD.OV-CTL.DESeq2.withTPM.xlsx",sheet = "plot") %>% pull(gene)


df.plot<-df[c(select_gene) %>% unique(),
  # setdiff(c(select_gene,deg) %>% unique(),deg[grep("Rik$",deg)]),
       c( "CTL_EGFP1", "CTL_EGFP2", "CTL_EGFP3",
          "OV_EGFP1",  "OV_EGFP2" , "OV_EGFP3" ,
          
          "CTL_TD1" ,  "CTL_TD2" ,  "CTL_TD3" ,
          "OV_TD1" ,   "OV_TD2"  ,  "OV_TD3"
          )] 

pdf(file = "out/Fig.7G TD.OV-CTL.deg.heatmap.withCTL.pdf",width = 4.2,height = 7)
bk <- c(seq(-2,2,by=0.01))
pheatmap::pheatmap(df.plot,scale = "row",
                   color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),
                             colorRampPalette(colors = c("white","#d23918"))(length(bk)/2)),
                    clustering_method = "ward.D2",
                   cluster_cols = F,cluster_rows =T,
                   gaps_col = c(6,9),
                   cutree_rows = 2,
                   treeheight_row = 20,
                   # gaps_row = 20
)
dev.off()