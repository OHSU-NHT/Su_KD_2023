
library(Seurat)
library(dplyr)
library(here)

# Fig. 5A ------------------------------------------------------------------------------------------------
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
        cols = mycols.DCT,
        split.by = "Diet") &
  scale_x_reverse()

ggsave(
  "F5A.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 6,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

