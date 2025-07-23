
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)

# Fig.4B ------------------------------------------------------------------------------------------------
SO <- readRDS(here("ALL_DCT.rds"))
Idents(SO) <- "subclass.l1"
my_levels <- c("DCT1", "DCT2", "Prolif")
Idents(SO) <- factor(x = Idents(SO), levels = my_levels)

mycols.DCT <- c("#66c2a5", # DCT1
                "#fc8d62", # DCT2
                "#A91554") # Prolif

DimPlot(SO,
        reduction = "umap",
        pt.size = 0,
        cols = mycols.DCT) +
  scale_x_reverse()

ggsave(
  "F4B.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)


# Fig.4C ------------------------------------------------------------------------------------------------

VlnPlot(SO,
        features = c("DCT1_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F4C_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

VlnPlot(SO,
        features = c("DCT2_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F4C_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4D ------------------------------------------------------------------------------------------------

DCT1_score_gene_list <- (c("Erbb4", "Egf", "Trpm7", "Fgf13", "Col5a2", "Umod", "Ptgfr", "Stk32b", "Rtl4", "Abca13"))
DCT2_score_gene_list <- (c("Slc8a1", "Arl15", "Calb1", "Slc2a9", "Phactr1", "Gls", "S100g", "Kl", "Klk1", "Egfem1"))

DotPlot(SO,
        features = DCT1_score_gene_list,
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        scale.max = 100,
        scale.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  RotatedAxis() +
  coord_flip() + 
  scale_y_discrete(limits = c("DCT1", "DCT2"))

ggsave(
  "F4D_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 3.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

DotPlot(SO,
        features = DCT2_score_gene_list,
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        scale.max = 100,
        scale.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  RotatedAxis() +
  coord_flip() + 
  scale_y_discrete(limits = c("DCT1", "DCT2"))

ggsave(
  "F4D_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 3.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4E ------------------------------------------------------------------------------------------------

VlnPlot(SO,
        features = c("ENaC_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F4E_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

VlnPlot(SO,
        features = c("Ca_score1"),
        split = TRUE,
        split.by = "Diet",
        cols = c("white", "lightskyblue"),
        flip = TRUE,
        pt.size = 0) & 
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = .8),
        plot.title = element_text(hjust = 0.5)) &
  stat_summary(fun = mean,
               geom = "crossbar",
               width = 0.5,
               size = 0.3,
               position = position_dodge(width = 0.5)) &
  scale_x_discrete(limits=c("DCT1", "DCT2"))

ggsave(
  "F4E_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 2.2,
  height = 3.5,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4F ------------------------------------------------------------------------------------------------

ENaC_score_gene_list <- (c("Nr3c2", "Klk1", "Scnn1a", "Scnn1b", "Scnn1g"))
Ca_score_gene_list <- (c("Slc8a1", "Vdr", "Trpv5", "Calb1", "S100g", "Ryr2"))

DotPlot(SO,
        features = ENaC_score_gene_list,
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        scale.max = 100,
        scale.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  RotatedAxis() +
  coord_flip() + 
  scale_y_discrete(limits = c("DCT1", "DCT2"))

ggsave(
  "F4F_1.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 3.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

DotPlot(SO,
        features = Ca_score_gene_list,
        cols = c("light grey", "royal blue"),
        dot.scale = 8,
        dot.min = 0,
        scale.max = 100,
        scale.min = 0,
        col.min = -2.5,
        col.max = 2.5) + 
  theme(axis.text.x = element_text(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  RotatedAxis() +
  coord_flip() + 
  scale_y_discrete(limits = c("DCT1", "DCT2"))

ggsave(
  "F4F_2.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 3.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

