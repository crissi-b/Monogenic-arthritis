library(Seurat)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(dplyr)
library(magrittr)
library(grid)
library(gridExtra)
library(ggExtra)
library(stringr)
library(data.table)
library(RColorBrewer)
library(splitstackshape)
library(harmony)
library("scProportionTest")
options(bitmapType='cairo')

# Load amp2 dataset (Zhang et al 2023)- called amp2_new_f
load(file= "amp2_full.RData")

# Load JIA tissue dataset (Bolton et al 2025)- called tissue
load(file = "sc-tissue.RData")

# Load adult arthritis dataset (Edalat et al 2024)- called synov
load(file="final-subclusters.RData")

#######################################################
#######################################################
#######################################################

# Subset amp2 data for biopsies only

Idents(amp2_new_f) <- "sample"
biopsy <- c('BRI-401', 'BRI-405', 'BRI-407', 'BRI-409', 'BRI-411', 'BRI-413', 'BRI-415', 'BRI-417', 'BRI-419',
'BRI-421', 'BRI-425', 'BRI-427', 'BRI-429', 'BRI-431', 'BRI-436', 'BRI-440', 
'BRI-458', 'BRI-462', 'BRI-475', 'BRI-481', 'BRI-483', 'BRI-485', 'BRI-503', 
'BRI-505', 'BRI-507', 'BRI-509', 'BRI-511', 'BRI-515', 'BRI-521', 'BRI-527', 'BRI-529', 
'BRI-534', 'BRI-536', 'BRI-538', 'BRI-540', 'BRI-542', 'BRI-544', 'BRI-546', 'BRI-548', 
'BRI-552', 'BRI-554', 'BRI-558', 'BRI-560', 'BRI-562', 'BRI-564', 'BRI-566', 'BRI-570', 
'BRI-581', 'BRI-583', 'BRI-585', 'BRI-589', 'BRI-601', 'BRI-603', 'BRI-605', 'BRI-607', 'BRI-611', 'BRI-613', 'BRI-615', 'BRI-617', 'BRI-619', 'BRI-621', 'BRI-623', 'BRI-625', 'BRI-627', 'BRI-629', 'BRI-631', 'BRI-633R', 'BRI-635', 'BRI-637R')

amp2 <- subset(amp2_new_f, sample %in% biopsy)
rm(amp2_new_f)

#######################################################
#######################################################
#######################################################

# Label plasma cells in AMP2 data and rename to give a common naming approach
Idents(amp2) <- "cluster_name"
rest <- subset(amp2, idents=c("B-2: IgG1+IgG3+ plasma", "B-6: IgM+ plasma"), invert=T)
plasma <- subset(amp2, idents=c("B-2: IgG1+IgG3+ plasma", "B-6: IgM+ plasma"))
plasma$cell_type <- "Plasma cells"
new <- rbind(plasma[["cell_type"]], rest[["cell_type"]])
amp2 <- AddMetaData(amp2, new, col.name="cell_type2")
Idents(amp2) <- "cell_type2"
amp2 <- RenameIdents(amp2, "T cell"="T cells", "NK"="NK cells / ILCs", "B cell/plasma cell"="B cells", "Myeloid cell"= "Myeloid cells", 
                     "Stromal cell"="Stromal cells", "Endothelial cell"="Endothelial cells")
amp2$cell_type2 <- amp2@active.ident

#######################################################
#######################################################
#######################################################

# Transfer cell labels to JIA dataset for finer subclusters
transfer.anchors <- FindTransferAnchors(reference = amp2, query = tissue, features = VariableFeatures(object = amp2), reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = amp2$cluster_name, weight.reduction = tissue[["pca"]], dims = 1:30)
tissue <- AddMetaData(tissue, metadata = celltype.predictions, col.name="predicted.id")
tissue$predicted.id -> tissue$RA_clusters

# Transfer cell labels to JIA dataset for major cell types 
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = amp2$cell_type2, weight.reduction = tissue[["pca"]], dims = 1:30)
tissue <- AddMetaData(tissue, metadata = celltype.predictions, col.name="predicted.id")
tissue$predicted.id -> tissue$main_type_RA

#######################################################
#######################################################
#######################################################

# Repeat for Edalat adult arthritis dataset
DefaultAssay(synov) <- "RNA"

synov <- synov %>% ScaleData(features = rownames(synov)) %>% FindVariableFeatures()
synov <- RunPCA(synov, features = VariableFeatures(synov), reduction.name = "pca", verbose = FALSE)

amp2 <- amp2 %>% ScaleData(features = rownames(amp2)) %>% FindVariableFeatures() 
amp2 <- RunPCA(amp2, features = VariableFeatures(amp2), reduction.name = "pca", verbose = FALSE)

# Confirm match
all(colnames(synov[["pca"]]) %in% colnames(synov)) 

#####################

transfer.anchors <- FindTransferAnchors(reference = amp2, query = synov, 
                                      features = VariableFeatures(object = amp2),
                                      reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = amp2$cluster_name,
                                     weight.reduction = synov[["pca"]], dims = 1:30)

synov <- AddMetaData(synov, metadata = celltype.predictions)
