#########################################################################################
#### Fig.S2A
#########################################################################################

library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(patchwork)

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