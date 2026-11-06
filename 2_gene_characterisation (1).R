
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
library(RColorBrewer)
library(pheatmap)
library(patchwork)

options(bitmapType='cairo')

#######################################################
#######################################################
#######################################################

# Plot tissue transcription of genes 
# Import Gene Summary data from Data Object sheet 1 and reformat 
syndromes <- read.csv(file="CSV/monogenic_sheet_1.csv", header=T, blank.lines.skip = TRUE)

excluded_genes <- syndromes %>% filter(exclude == "Y") %>% pull(Gene) # Remove genes with low penetrance
gen <- syndromes %>% filter(exclude == "N") %>% pull(Gene)

#######################################################
#######################################################
#######################################################

syndromes_tb <- syndromes %>% filter(exclude=="N")

# Features
syn_anno <- subset(syndromes_tb, select=rev(c("Gene_D", 
                                              "Autoinflammatory", "Uveitis", "Psoriasis",
                                              "SLE", "IBD","HLH...MAS", "Vasculitis_Panniculitis", "Endocrinopathy",
                                              "Syndrome")))

# Replace blank with none
syn_anno <- syn_anno %>% mutate(across(c(Syndrome), ~ replace_na(str_trim(na_if(.x, "")), "None")))
rownames(syn_anno) <-syn_anno$Gene_D
syn_anno$Gene_D <- NULL
syn_anno[,c(2:9) ][is.na(syn_anno[,c(2:9) ])] <- 0 # Convert NA values to 0

# Prepare annotation data frame for pheatmap
syn_anno$Syndrome <- as.factor(syn_anno$Syndrome)

syn_anno$Vasculitis_Panniculitis <- as.numeric(syn_anno$Vasculitis_Panniculitis)
syn_anno$HLH...MAS <- as.numeric(syn_anno$HLH...MAS)
syn_anno$Autoinflammatory <- as.numeric(syn_anno$Autoinflammatory)
syn_anno$IBD <- as.numeric(syn_anno$IBD)
syn_anno$SLE <- as.numeric(syn_anno$SLE)
syn_anno$Uveitis <- as.numeric(as.character(syn_anno$Uveitis))
syn_anno$Psoriasis <- as.numeric(syn_anno$Psoriasis)

syndrome_colours <- c("Hyper IgM" = "purple3", "Wiskott-Aldrich /-like"  = "darkolivegreen2","LSD" = "grey45", "None"="white",
                      "Hyper IgD"  = "deeppink3",  "IPEX /-like"  = "cyan", "Hypomorphic SCID"  = "deepskyblue3",  
                      "AGS /-like"  = "cornsilk3", "FCAS" = "lightcoral",  "EDA-ID /-like" = "forestgreen", "CGD"= "goldenrod4")

# Combine all annotation colours into a list for pheatmap
annotation_colours_list <- list(
  Syndrome = syndrome_colours)

anno_gene <- syn_anno[,c(1), drop=F]

ph <- pheatmap(t(syn_anno[,c(2:9)]), color = c("white", "orange"), na_col = "white", 
               fontsize=7, border_color = "lightgrey", legend = F, 
               clustering_distance_cols = "binary",
               annotation_col = anno_gene, annotation_colors = annotation_colours_list)

#######################################################
#######################################################
#######################################################

arth <- read.csv(file="CSV/arthritis_csv.csv", header=T, row.names = 1)
arth[arth == ""] <- 0
arth[arth == "NR"] <- NA

arth <- arth %>%
  mutate(across(everything(), as.numeric))

arth$RF_status <- with(arth,
                       ifelse(RF_pos == 1 & RF_neg == 1, "Both",
                              ifelse(RF_pos == 1 & RF_neg == 0, "RF_pos",
                                     ifelse(RF_pos == 0 & RF_neg == 1, "RF_neg",
                                            "Neither"))))

arth$ANA_status <- with(arth,
                        ifelse(ANA_pos == 1 & ANA_neg == 1, "Both",
                               ifelse(ANA_pos == 1 & ANA_neg == 0, "ANA_pos",
                                      ifelse(ANA_pos == 0 & ANA_neg == 1, "ANA_neg",
                                             "Neither"))))

arth$eros_status <- with(arth,
                         ifelse(Erosions == 1 & Non.erosive == 1, "Both",
                                ifelse(Erosions == 1 & Non.erosive == 0, "Erosions",
                                       ifelse(Erosions == 0 & Non.erosive == 1, "Non-erosive",
                                              "Neither"))))

arth$size <- with(arth,
                  ifelse(Large.joints == 1 & Small.joints == 1, "Both",
                         ifelse(Large.joints == 1 & Small.joints == 0, "Large_joints",
                                ifelse(Large.joints == 0 & Small.joints == 1, "Small_joints",
                                       "Neither"))))


# Plot number of genes
plot_data <- arth %>%
  dplyr::select(RF_pos, RF_neg, ANA_pos, ANA_neg) %>%
  summarise(across(everything(), ~sum(. == 1, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = c("Marker", "Status"), names_sep = "_")

plot_data$Marker <- factor(plot_data$Marker, levels=rev(c("ANA", "RF")))

p1 <- ggplot(plot_data, aes(y = Marker, x = value, fill = Status)) +
  geom_col(width = 0.7) + geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("pos" = "#00A087", "neg" = "gold2")) +
  labs(y = NULL, x = "Number of genes", fill = "Status") +
  theme_classic(base_size = 14) + theme(axis.text = element_text(color = "black"))

arth.size <- arth[,rev(c("Polyarthritis", "Oligoarthritis", "Monoarthritis")), drop=F]

gene_counts <- data.frame(Arthritis = colnames(arth.size), Count = colSums(arth.size == 1, na.rm = TRUE))

p2 <- ggplot(gene_counts, aes(y = Arthritis, x = Count, fill = Arthritis)) +
  geom_col(width = 0.7) + 
  scale_fill_manual(values = c("#4DBBD5", "#E64B35", "#00A087")) +
  labs(y = NULL, x = "Number of genes") +
  theme_classic(base_size = 14) + theme(legend.position = "none", axis.text = element_text(color = "black")) +
  expand_limits(x = max(gene_counts$Count) * 1.1)

plot_data <- arth %>%
  dplyr::select(Large.joints, Small.joints, Axial...spondyloarthritis) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "Status") %>%
  count(Group, Status) %>% na.omit()

plot_data <- data.frame(Group = c("Large.joints", "Small.joints", "Axial...spondyloarthritis"),
                        Count = colSums(arth[, c("Large.joints", "Small.joints", "Axial...spondyloarthritis")] == 1, na.rm = TRUE))

p3 <- ggplot(plot_data, aes(y = Group, x = Count)) +
  geom_col(width = 0.7) +
  labs(y = NULL, x = "Number of genes", fill = "Status") +
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(color = "black"))


# Combine horizontally
combined_plot <- (p1 | p2 | p3) + plot_layout(widths = c(1, 1, 1))
combined_plot

#######################################################
#######################################################
#######################################################

arth2 <- arth[,rev(c("size","Axial...spondyloarthritis", "Monoarthritis", "eros_status",  "RF_status", "ANA_status")), drop=F]

arth2$Monoarthritis <- as.character(arth2$Monoarthritis)
arth2$Axial...spondyloarthritis <- as.character(arth2$Axial...spondyloarthritis)

syn_anno2 <- merge(anno_gene, arth2, by="row.names")
syn_anno2 <- syn_anno2 %>% column_to_rownames("Row.names")

anno_col <- c("RF_neg"="grey90", "RF_pos"="red", "Both"="orange", Neither="white")
anno_col2 <- c("ANA_neg"="grey90", "ANA_pos"="goldenrod3", "Both"="gold", Neither="white")
anno_col3 <- c("Erosions"="blue", "Non-erosive"="grey90", "Both"="purple", Neither = "white")  
anno_col6 <- c("0"="grey90", "1"="hotpink2", "NR"="white")
anno_col7 <- c("0"="grey90", "1"="green4", "NR"="white")
anno_col9 <- c("Large_joints"="turquoise3", "Small_joints"="pink2", "Both"="lavender", "Neither"="white")

annotation_colours_list <- list(
  eros_status = anno_col3,
  RF_status = anno_col,
  ANA_status = anno_col2,
  Axial...spondyloarthritis=anno_col7,
  Monoarthritis = anno_col6, 
  size=anno_col9,
  Syndrome=syndrome_colours)

ph <- pheatmap(t(syn_anno[,c(2:9)]), color = c("white", "grey50"), na_col = "white", 
               fontsize=7, border_color = "lightgrey", legend = F, 
               clustering_distance_cols = "binary", clustering_distance_rows = "binary",
               annotation_col = syn_anno2, annotation_colors = annotation_colours_list)

#######################################################
#######################################################
#######################################################


