#########################################################################################
#### Fig.S7A
#########################################################################################

data.tlineage <- readRDS("~/Project/YuYong_TALL/data.Tlineage.assign.rds")

data.tlineage$celltype <- factor(gsub("T\\.", "", data.tlineage@active.ident), levels = c("HSC", "MLP", "MPP", "ETP", "DN2", "DN3", "DN4", "ISP", "DP", "NKT", "Tgd"))
data.tlineage@meta.data[data.tlineage@meta.data$celltype_singleR_ImmGenData == "Stem cells (SC.MPP34F)", ]$celltype <- "MPP"
DimPlot(data.tlineage, label = T, group.by = "celltype")
ggsave(filename = "Tlineage.celltype.UMAP.pdf", width = 5, height = 4.5)


seu <- data.tlineage
seu@active.assay <- "RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu)


### 将normal的etp迁移到细胞上

normal.celltype <- normal@meta.data %>%
  dplyr::select(celltype_normal, orig.ident) %>%
  tibble::rownames_to_column("cell")


normal.tl <- seu@meta.data %>%
  filter(orig.ident %in% c("DN", "THY")) %>%
  dplyr::select(celltype, orig.ident) %>%
  tibble::rownames_to_column("cell")

normal.tl$cell_map <- paste0(normal.tl$cell, "_2")
intersect(normal.tl$cell_map, normal.celltype$cell)
# ok的 进行合并
df <- merge(normal.tl, normal.celltype, all.x = T, by.x = "cell_map", by.y = "cell")

rownames(df) <- df$cell


table(df$celltype, df$celltype_normal) %>% pheatmap::pheatmap(., scale = "column")


## 迁移回tlineage
seu$celltype_new <- as.character(seu$celltype)
seu$celltype_new[df$cell] <- as.character(df$celltype_normal)
seu$celltype_new <- factor(seu$celltype_new, levels = c("HSC", "MLP", "MPP", "ETP", "DN2", "DN3", "DN4", "ISP", "DP", "CD8+", "NKT", "Tgd"))
seu$celltype_new %>% table()
# only 1 cell for Tgd and ISP
# so remove it
seu.new <- subset(seu, celltype_new %in% c("HSC", "MLP", "MPP", "ETP", "DN2", "DN3", "DN4", "DP", "CD8+", "NKT"))

DimPlot(seu.new, group.by = "celltype_new", label = T)
ggsave(filename = "Tlineage.normal-reAssign.celltype.UMAP.pdf", width = 5, height = 4.5)

hsc <- c(
  "Runx1", "Mllt3", "Hoxa9", "Mecom", "Hlf",
  # "Spink2",

  # "Epor",
  "Mpl", "Kit", "Il3ra"
)
MPP1 <- c(
  "Cd34", # cellmarker2 来源MPP marker
  "Ly6a",
  "Kit",
  "Fcgr3", "Fcgr4", #' Cd16',
  "Fcgr2b", #' Cd32',
  "Flt3", # CD135
  "Gcnt2",
  "Hlf",
  "Il7r"
)
MPP2 <- c(
  "Cd34", "Ly6a", "Kit", "Flt3", "Gcnt2", "Hlf",
  "Gata2"
)

DotPlot(seu.new,
  features = c(
    "Cd34", "Hoxa9", "Mecom", "Hlf",
    "Fcgr3", "Fcgr4", #' Cd16',
    "Fcgr2b", #' Cd32',
    "Kit", "Il7r", "Il2ra", "Cd4", "Cd8a", "Cd8b1",
    "Bcl11a", "Bcl11b", "Notch1",
    "Rag1", "Rag2", "Gata3", "Dntt", "Cd247", "Ptprc",
    # "Cd52",
    "Cd3d",
    # "Cd3e",
    "Tcf7", "Klrb1c", "Ccl5"
  ) %>% unique(),
  cluster.idents = F, group.by = "celltype_new"
) +
  xlab("") + ylab("") + theme_bw() +
  theme(
    axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )

ggsave(filename = "Lineage.normal-reAssign.celltype.marker.pdf", width = 8, height = 3.5)




MTX <- LayerData(seu.new, assay = "RNA", layer = "data")
seu.new$Scin_expression <- MTX["Scin", ]

# VlnPlot(seu,features = "Scin_expression",group.by = "celltype")

df <- seu.new@meta.data %>% dplyr::select(celltype_new, Scin_expression)
df$celltype <- df$celltype_new
df <- df %>%
  group_by(celltype) %>%
  mutate(Scin = mean(Scin_expression)) %>%
  dplyr::select(celltype, Scin) %>%
  unique.data.frame()
df$celltype <- factor(df$celltype, levels = levels(seu.new$celltype_new))

# plot
ggplot(seu.new@meta.data %>% as.data.frame(), aes(x = celltype_new, y = Scin_expression)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1, color = "grey20") +
  geom_point(data = df, aes(x = celltype, y = Scin), color = "salmon") +
  geom_line(data = df, aes(x = celltype, y = Scin), group = 1, color = "salmon") +
  theme_classic()

ggsave(filename = "Fig.S7A Scin.expression.pdf", height = 3, width = 6)
