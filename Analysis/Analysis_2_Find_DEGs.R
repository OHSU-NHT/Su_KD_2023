
library(Seurat)
library(here)
library(DESeq2)
library(dplyr)
library(openxlsx)

# 1. Computing potassium-deficiency-induced DEGs ---------------------------------
SO <- readRDS(here("ALL_DCT.rds"))
Idents(SO) <- "analysis"

# DCT1
DCT1.kd <- FindMarkers(SO,
                       ident.1 = "DCT1_KD",
                       ident.2 = "DCT1_NK",
                       test.use = "DESeq2",
                       assay = "SCT",
                       min.pct = 0.1,
                       logfc.threshold = 0.3219)

# DCT2
DCT2.kd <- FindMarkers(SO,
                       ident.1 = "DCT2_KD",
                       ident.2 = "DCT2_NK",
                       test.use = "DESeq2",
                       assay = "SCT",
                       min.pct = 0.1,
                       logfc.threshold = 0.3219)

# Proliferation population
Prolif.kd <- FindMarkers(SO,
                         ident.1 = "Prolif_KD",
                         ident.2 = "Prolif_NK",
                         test.use = "DESeq2",
                         assay = "SCT",
                         min.pct = 0.1,
                         logfc.threshold = 0.3219)

# 2. Save all DEG lists -----------------------------------------------------
# data.frame
save(DCT1.kd, DCT2.kd, Prolif.kd, file = here("KD DEGs.RData"))
save(DCT1.kd, DCT2.kd, Prolif.kd, file = here("KD DEGs.rda"))

# excel file
dataset_names <- list('DCT1.kd' = DCT1.kd, 'DCT2.kd' = DCT2.kd, 'Prolif.kd' = Prolif.kd,)
write.xlsx(dataset_names, file = 'KD DEGs.xlsx', rowNames=TRUE)
