#########################################################################################
#### Fig.S1G
#########################################################################################

# setwd("out/")
library(ggplot2)
library(dplyr)
library(tidyr)

tcr_dn1<-read.csv("cellranger/20221123_DN_TCR/outs/clonotypes.csv")
tcr_dn2<-read.csv("cellranger/20221127_DN_TCR/outs/clonotypes.csv")
tcr_dn5<-read.csv("cellranger/20240126_DN_TCR/clonotypes.csv")

tcr_thy1<-read.csv("cellranger/20221123_Thy_TCR/outs/clonotypes.csv")
tcr_thy2<-read.csv("cellranger/20221127_Thy_TCR/outs/clonotypes.csv")
tcr_thy5<-read.csv("cellranger/20240126_Thy_TCR/clonotypes.csv")

## stack bar
df_dn1<-tcr_dn1 %>% mutate(Proportion=100*frequency/sum(tcr_dn1$frequency)) %>% head(12)
df_dn2<-tcr_dn2 %>% mutate(Proportion=100*frequency/sum(tcr_dn2$frequency)) %>% head(12)
df_dn5<-tcr_dn5 %>% mutate(Proportion=100*frequency/sum(tcr_dn5$frequency)) %>% head(12)

df_thy1<-tcr_thy1 %>% mutate(Proportion=100*frequency/sum(tcr_thy1$frequency)) %>% head(12)
df_thy2<-tcr_thy2 %>% mutate(Proportion=100*frequency/sum(tcr_thy2$frequency)) %>% head(12)
df_thy5<-tcr_thy5 %>% mutate(Proportion=100*frequency/sum(tcr_thy5$frequency)) %>% head(12)

cdr3_comm<-intersect(df_dn1$cdr3s_aa,df_thy1$cdr3s_aa)
intersect(df_dn5$cdr3s_aa,c(df_dn1$cdr3s_aa,df_thy1$cdr3s_aa))

df<-data.frame(Group=c(rep(c("DN_Normal","DN_TALL","Thy_TALL"),each=length(cdr3_comm))),
               Clonotype=c(cdr3_comm,
                           df_dn1 %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(cdr3s_aa),
                           df_thy1 %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(cdr3s_aa)),
               Proportion=c(rep(0.001,length(cdr3_comm)),
                            df_dn1 %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(Proportion),
                            df_thy1 %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(Proportion)))


p1<-ggplot(df, aes(x =Group, y= Proportion, fill = Clonotype,
                  stratum=Clonotype, alluvium=Clonotype)) +
  
  geom_flow(width=0.5,alpha=0.6, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  geom_col(width = 0.5, color='black',linewidth=0.7)+
  theme_classic() +
  labs(x='',y = 'Clonetype proportion (%)')+
  ggtitle("Sample1")+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 12,color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5,size = 12,color = "black"),
    axis.title.x  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
  )

head_n<-155
df_dn2<-tcr_dn2 %>% mutate(Proportion=100*frequency/sum(tcr_dn2$frequency)) %>% head(head_n) %>% distinct(cdr3s_aa,.keep_all = T)
df_dn5<-tcr_dn5 %>% mutate(Proportion=100*frequency/sum(tcr_dn5$frequency)) %>% head(head_n) %>% distinct(cdr3s_aa,.keep_all = T)

df_thy2<-tcr_thy2 %>% mutate(Proportion=100*frequency/sum(tcr_thy2$frequency)) %>% head(head_n) %>% distinct(cdr3s_aa,.keep_all = T)
df_thy5<-tcr_thy5 %>% mutate(Proportion=100*frequency/sum(tcr_thy5$frequency)) %>% head(head_n) %>% distinct(cdr3s_aa,.keep_all = T)
cdr3_comm<-intersect(df_dn2$cdr3s_aa %>% unique(),df_thy2$cdr3s_aa %>% unique())
cdr3_comm
intersect(df_dn5$cdr3s_aa,c(df_dn2$cdr3s_aa,df_thy2$cdr3s_aa))

df<-data.frame(Group=c(rep(c("DN_Normal","DN_TALL","Thy_TALL"),each=length(cdr3_comm))),
               Clonotype=c(cdr3_comm,
                           df_dn2 %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(cdr3s_aa),
                           df_thy2 %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(cdr3s_aa)),
               Proportion=c(rep(0.001,length(cdr3_comm)),
                            df_dn2 %>% head(head_n)%>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(Proportion),
                            df_thy2 %>% head(head_n) %>% filter(cdr3s_aa %in% cdr3_comm) %>% pull(Proportion)))


p2<-ggplot(df, aes(x =Group, y= Proportion, fill = Clonotype,
                   stratum=Clonotype, alluvium=Clonotype)) +
  
  geom_flow(width=0.5,alpha=0.6, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  geom_col(width = 0.5, color='black',linewidth=0.7)+
  theme_classic() +
  labs(x='',y = 'Clonetype proportion (%)')+
  ggtitle("Sample2")+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 12,color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5,size = 12,color = "black"),
    axis.title.x  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
  )
p2

library(patchwork)
p1/p2
ggsave(plot=p1/p2,filename = "out/Fig.S1G TCR top10 FlowBar.pdf",width = 10,height = 8)


#########################################################################################
#### Fig.S1H
#########################################################################################


# setwd("~/Project/YuYong_TALL/out/")
library(ggplot2)
library(dplyr)
library(tidyr)

tcr_dn1<-read.csv("../cellranger/20221123_DN_TCR/outs/clonotypes.csv")
tcr_dn2<-read.csv("../cellranger/20221127_DN_TCR/outs/clonotypes.csv")
tcr_dn5<-read.csv("../cellranger/20240126_DN_TCR/clonotypes.csv")

tcr_thy1<-read.csv("../cellranger/20221123_Thy_TCR/outs/clonotypes.csv")
tcr_thy2<-read.csv("../cellranger/20221127_Thy_TCR/outs/clonotypes.csv")
tcr_thy5<-read.csv("../cellranger/20240126_Thy_TCR/clonotypes.csv")

table<-read.csv("../cellranger/20221127_DN_TCR/outs/consensus_annotations.csv")
table2<-read.csv("../cellranger/20221123_Thy_TCR/outs/filtered_contig_annotations.csv")

# 计算D50的函数
calculate_D50 <- function(abundance) {
  # 按降序排序
  r_sorted <- sort(abundance, decreasing = TRUE)
  J <- sum(r_sorted)
  if (J == 0) stop("The sum of abundances is zero, cannot calculate D50.")
  half_J <- J / 2
  cum_r <- cumsum(r_sorted)
  # 找到第一个满足累积和 ≥ J/2 的索引
  C <- which(cum_r >= half_J)[1]
  S <- length(r_sorted)
  D50 <- (C / S) * 100
  return(D50)
}

# sample_data <- c(50, 30, 15, 5)  # D50应为25
# calculate_D50(sample_data) 

d50.dn1<-calculate_D50(tcr_dn1$frequency) 
d50.dn2<-calculate_D50(tcr_dn2$frequency) 
d50.dn5<-calculate_D50(tcr_dn5$frequency) 
d50.thy1<-calculate_D50(tcr_thy1$frequency) 
d50.thy2<-calculate_D50(tcr_thy2$frequency) 
d50.thy5<-calculate_D50(tcr_thy5$frequency) 

df<-data.frame(group=c('DN_Normal',"DN_TALL","Thy_Normal","Thy_TALL"),
               type=c("Normal","TALL","Normal","TALL"),
               D50=c(d50.dn5,mean(c(d50.dn1,d50.dn2)),d50.thy5,mean(c(d50.thy1,d50.thy2))))

ggplot(data = df,aes(x=group,y=D50,fill=type))+
  geom_bar(stat = "identity",width = 0.6)+
  theme_classic()+
  labs(
    x = "Sample",
    y = "TCR diversity (D50)",
    fill = "Type"
  ) +
  scale_fill_manual(values = c("Normal" = "#4E79A7", "TALL" = "#E15759")) +
  theme(
    axis.text.x = element_text(angle = -30, hjust = 0,size = 12,color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5,size = 12,color = "black"),
    axis.title.x  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
  )
  

# PLOT
df_plot <- df %>%
  separate(group, into = c("main_group", "sub_group"), sep = "_", remove = FALSE)

ggplot(df_plot, aes(x = main_group, y = D50, fill = sub_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.65), width = 0.6) +
  scale_fill_manual(values = c("Normal" = "#4E79A7", "TALL" = "#E15759")) +
  labs(
    x = "Group",
    y = "TCR diversity (D50)",
    fill = "Type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 14,color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5,size = 12,color = "black"),
    axis.title.x  = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
  )
ggsave(filename = "out/Fig.S1H TCR.diversity.pdf",width = 4,height = 5)

