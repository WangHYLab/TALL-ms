#########################################################################################
#### Fig.S9D
#########################################################################################
setwd("cuttag-analysis")
library(tidyverse)
library(openxlsx)
library(eulerr)

Runx1_mm <- read.xlsx("ref_data/GSE218147/peakCalling_mm10/Runx1.xlsx")
Bcl11a_mm1 <- read.xlsx("cuttag-analysis/peakCalling/anno_bed/B1.xlsx")
Bcl11a_mm2 <- read.xlsx("cuttag-analysis/peakCalling/anno_bed/B2.xlsx")
Bcl11a_mm3 <- read.xlsx("cuttag-analysis/peakCalling/anno_bed/B3.xlsx")
Bcl11a_mm <- intersect(intersect(Bcl11a_mm1$SYMBOL, Bcl11a_mm2$SYMBOL), Bcl11a_mm3$SYMBOL)

## 0722 ----
bcl11a.deg <- read.xlsx("data/VCAB-DKO vs. VCAB.DESeq2.withTPM.xlsx")
bcl11a.deg$sig <- "ns"
bcl11a.up <- bcl11a.deg %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange > 0) %>%
  pull(gene) %>%
  unique()
bcl11a.dn <- bcl11a.deg %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange < 0) %>%
  pull(gene) %>%
  unique()

gene2des <- bcl11a.deg %>%
  dplyr::select(gene, description) %>%
  distinct(gene, .keep_all = T) %>%
  tibble::column_to_rownames("gene")

## Runx1 OE data come from GSE218147
# Runx factors launch T cell and innate lymphoid programs via direct and gene network-based mechanisms [https://www.nature.com/articles/s41590-023-01585-z]
# SUPPLEMENTARY TABLE 3: Total_list_of_differentially_expressed_genes_in_OE-diverging-UMAP_2_window_from_Figure_4
Runx1.deg <- read.xlsx("data/Venn.xlsx", sheet = "Runx1", startRow = 2, rowNames = T)
runx1.up <- Runx1.deg %>%
  filter(Note == "OE.UP") %>%
  rownames()
runx1.dn <- Runx1.deg %>%
  filter(Note == "OE.DN") %>%
  rownames()
gene <- intersect(Bcl11a_mm, Runx1_mm$SYMBOL)

df <- data.frame(
  Gene = gene,
  Description = gene2des[gene, ]
  # Bcl11a.UP=gene %in% bcl11a.up,
  # Bcl11a.DN=gene %in% bcl11a.dn,
  # Runx1.UP=gene %in% runx1.up,
  # Runx1.DN=gene %in% runx1.dn
)
df$Bcl11a <- ""
df[df$Gene %in% bcl11a.up, ]$Bcl11a <- "Up"
df[df$Gene %in% bcl11a.dn, ]$Bcl11a <- "Down"
df$Runx1 <- ""
df[df$Gene %in% runx1.up, ]$Runx1 <- "Up"
df[df$Gene %in% runx1.dn, ]$Runx1 <- "Down"

write.xlsx(df, file = "data/Bcl11a&Runx1.overlap-1450.gene.withDEGinfo.xlsx")

# 读取上面生成的基因表
df <- read.xlsx("data/Bcl11a&Runx1.overlap-1450.gene.withDEGinfo.xlsx")

Bcl11a_KO_up <- df %>%
  filter(Bcl11a == "Up") %>%
  pull(Gene)
Bcl11a_KO_dn <- df %>%
  filter(Bcl11a == "Down") %>%
  pull(Gene)

Runx_OE_up <- df %>%
  filter(Runx1 == "Up") %>%
  pull(Gene)
Runx_OE_dn <- df %>%
  filter(Runx1 == "Down") %>%
  pull(Gene)



gene_lists <- list(
  Bcl11a_KO_up = Bcl11a_KO_up,
  Bcl11a_KO_dn = Bcl11a_KO_dn,
  Runx_OE_up = Runx_OE_up,
  Runx_OE_dn = Runx_OE_dn
)
a <- euler(gene_lists)

pdf("out/Fig.S9D Bcl11a&Runx1.venn4.elipse.pdf", width = 5, height = 4)
plot(a,
  labels = TRUE, # 显示标签
  counts = TRUE, # 显示数量
  fill = RColorBrewer::brewer.pal(n = 8, name = "Paired")[c(1, 3, 5, 7)],
  alpha = 0.5
)
dev.off()


intersect(Bcl11a_KO_up, Runx_OE_dn)
intersect(Bcl11a_KO_up, Runx_OE_up)
intersect(Bcl11a_KO_dn, Runx_OE_dn)
