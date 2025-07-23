library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(here)
library(RColorBrewer)
library(monocle3)
library(SeuratWrappers)

# Fig.4A-----------------------------
SO.all <- readRDS(here("ALL_DCT.rds"))
Idents(SO.all) <- "subclass.l1"
mycols.DCT <- c("#66c2a5", # DCT1
                "#fc8d62", # DCT2
                "#A91554"  # Prolif
)

cells_with_Top2a <- WhichCells(SO.all, expression = Top2a > 0)
cells_with_Mki67 <- WhichCells(SO.all, expression = Mki67 > 0)
SO.all$Top2a_expr <- "noTop2a"
SO.all$Top2a_expr[cells_with_Top2a] <- "Top2a"
SO.all$Mki67_expr <- "noMki67"
SO.all$Mki67_expr[cells_with_Mki67] <- "Mki67"
meta.SO.all <- head(SO.all@meta.data)

SO.all@meta.data <- SO.all@meta.data %>%
  mutate(subset = ifelse(subclass.l1 == "Prolif" | 
                           Top2a_expr == "Top2a" | 
                           Mki67_expr == "Mki67", 
                         "proliferating", "no"))

Idents(SO.all) <- "subset"
DimPlot(SO.all, reduction = "umap", 
        cells.highlight = WhichCells(SO.all, idents = c("proliferating")),  # Highlight cells in selected clusters
        cols.highlight = "magenta",                                       # Color for the highlighted cluster
        cols = "lightgrey",                                                      # Color for other clusters
        pt.size = 1)+
  ggtitle("Proliferating cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  scale_x_reverse()                                                    

ggsave(
  "F6A_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 8,
  height = 7,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

FeaturePlot(SO.all, 
            features = c("Top2a"), 
            pt.size = 0.1,
            cols = c("lightgrey", "royalblue"),
            order = TRUE) & scale_x_reverse()
ggsave(
  "F6A_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

FeaturePlot(SO.all, 
            features = c("Mki67"), 
            pt.size = 0.1,
            cols = c("lightgrey", "royalblue"),
            order = TRUE) & scale_x_reverse()
ggsave(
  "F6A_3.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6B-----------------------------
SO <- subset(SO.all, subset = subset == "proliferating")
DefaultAssay(SO) <- "integrated"
SO_1 <- RunPCA(SO, verbose = FALSE)
ElbowPlot(SO_1, n=50)
SO_2 <- RunUMAP(SO_1, dims = 1:50)
SO_2 <- FindNeighbors(SO_2, reduction = "pca", dims = 1:50)
SO3 <- FindClusters(SO_2, resolution = 0.1)

SO3@meta.data <- SO3@meta.data %>% 
  mutate(Proliferation = dplyr::case_when(
    seurat_clusters == 0  ~ "DCT1",
    seurat_clusters == 1  ~ "DCT2",
    seurat_clusters == 2  ~ "Prolif",
    seurat_clusters == 3  ~ "pro-DCT1"
  ))
Idents(SO3) <- "Proliferation"
my_levels <- c("DCT1", "DCT2", "pro-DCT1", "Prolif")
Idents(SO3) <- factor(x = Idents(SO3), levels = my_levels)

mycols.prolif <- c("#66c2a5", # DCT1
                   "#fc8d62", # DCT2
                   "#1E88E5", # pro-DCT1
                   "#A91554") # Prolif

DimPlot(object = SO3, 
        reduction = "umap", 
        pt.size = 0.5,
        label = TRUE,
        label.size = 4,
        repel = TRUE,
        cols = mycols.prolif)

ggsave(
  "F6B_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

DimPlot(object = SO3, 
        reduction = "umap", 
        pt.size = 0.5,
        split.by = "Diet",
        cols = mycols.prolif,
        ncol = 1)

ggsave(
  "F6B_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 8,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)


# Fig.6C-----------------------------
Idents(SO3) <- "Proliferation"
my_levels <- c("DCT1", "DCT2", "pro-DCT1", "Prolif")
Idents(SO3) <- factor(x = Idents(SO3), levels = my_levels)
mycols.prolif <- c("#66c2a5", # DCT1
                   "#fc8d62", # DCT2
                   "#1E88E5", # pro-DCT1
                   "#A91554") # Prolif

t1 <- table(Idents(SO3), SO3$Rep)[, c("NK1", "NK2", "NK3", "KD1", "KD2", "KD3")]
t1 <- as.data.frame(t1)
colnames(t1) <- c('Cell_type', 'Rep', 'Frequency')

t2 <- t1 %>%
  dplyr::mutate(Diet = ifelse(grepl("NK", Rep, ignore.case = TRUE), "NK", "KD"))
t3 <- t2 %>%
  group_by(Cell_type, Diet) %>%
  summarise(mean_freq = mean(Frequency, na.rm=T),
            sd_freq = sd(Frequency, na.rm = TRUE)
  )
t4.norm <- t3 %>%
  group_by(Diet) %>%
  mutate(mean_freq_norm = mean_freq / sum(mean_freq)) %>%
  mutate(sd_freq_norm = sd_freq / sum(mean_freq)) %>%
  ungroup()

t4.norm$Diet <- factor(t4.norm$Diet, levels = c("KD", "NK"))
ggplot(t4.norm,
                    aes(fill=Cell_type, x=mean_freq_norm, y=Diet)) + 
  geom_bar(position="fill", stat = "identity", colour="black") +
  theme_classic() +
  labs(x = "Proportion") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = mycols.prolif)

ggsave(
  "F6C.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 8,
  height = 2,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6D-----------------------------
cds <- as.cell_data_set(SO3)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
all_partitions <- unique(cds@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions %in% all_partitions] <- "1"
cds <- learn_graph(cds, use_partition = F)

get_earliest_principal_node <- function(cds, time_bin="DCT"){
  cell_ids <- which(colData(cds)[, "class"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds = cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1.5,
           label_groups_by_cluster=FALSE,
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           graph_label_size=5,
           cell_size = 0.5) + scale_color_viridis_c()

ggsave(
  "F6D.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 5,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.6E-----------------------------
SO3 <- PrepSCTFindMarkers(SO3, verbose = TRUE)
SO3.prolif.markers <- FindAllMarkers(SO3, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)

top5 <- SO3.prolif.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5
DoHeatmap(SO3,
          features = top5$gene,
          label = TRUE,
          size = 3,
          hjust = 0.5,
          angle = 0,
          group.colors = mycols.prolif) + 
  NoLegend() +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 3, name = "RdBu"))) +
  theme(axis.text.y = element_text(size = 12))
ggsave(
  "F6E.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 12,
  height = 4,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)


