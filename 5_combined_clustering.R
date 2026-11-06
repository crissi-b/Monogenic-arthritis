
library(SingleCellExperiment)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(scales)
library(ggpubr)
library(viridis)
library(vegan)

col_fun <- c("turquoise4", "turquoise3", "#AAEBEE", "#FFFFFF", "#D8C1B0", "#8B4513", "brown4")
heatcol2 <- colorRampPalette(c("blue", "dodgerblue", "lightblue1","white","white", "darkorange","red"))(10)


options(bitmapType='cairo')

#######################################################

# Load scRNA-seq expression data for adults and children, pathway outputs and syndromic features from monogenic_sheet_1
# Combine together

#######################################################

load("2601_pathway-data.RData")

load(file="2609_expr_clustering_objects.RData")

#######################################################
#######################################################
#######################################################

syndromes <- read.csv(file="CSV/monogenic_sheet_1.csv", header=T, blank.lines.skip = TRUE)
syndromes <- syndromes %>% filter(exclude == "N")
included_genes <- syndromes %>% pull(Gene_D)

# Pull top pathways
filtered <- all_enrichment %>% mutate(BgNum = as.numeric(str_extract(BgRatio, "^[0-9]+"))) %>% filter(BgNum < 300)
uni_gene <- filtered %>% distinct(geneID, .keep_all = TRUE)
top200 <- uni_gene %>% arrange(qvalue) %>% slice_head(n = 300)
top200 -> pw_list

# Separate geneID column into individual genes per row
combi_long <- pw_list %>% dplyr::select(Description, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>% unnest(geneID)

# Create binary matrix
gene_pathway_matrix <- combi_long %>%
  mutate(value = 1) %>% pivot_wider(names_from = geneID, values_from = value, values_fill = 0)

# Set pathway descriptions as row names (optional)
gene_matrix <- gene_pathway_matrix %>% column_to_rownames("Description") %>% as.matrix()

colnames(gene_matrix)[colnames(gene_matrix) == "STING1"] <- "TMEM173"
colnames(gene_matrix)[colnames(gene_matrix) == "RIGI"] <- "DDX58"

# Remove infection and unrelated disease terms
remove <- rownames(gene_matrix)[grepl("Measles|Tumor|Diabet|Yersinia|Hepatitis|Alcoholic|leukaemia|Toxoplasmosis|Herpes|Chagas|Influenza|Sars|Pertussis|Leishmania|covid|cancer|tuberculosis", rownames(gene_matrix), ignore.case = TRUE)]
gene_matrix <- gene_matrix[!rownames(gene_matrix) %in% remove,]

gene_matrix %>% as.data.frame() -> matrix_order 

matrix_order <- matrix_order %>%
  mutate(`PIK3R1_LOF` = PIK3R1, `PIK3R1_GOF` = PIK3R1, 
         `PIK3CD_LOF` = PIK3CD, `PIK3CD_GOF` = PIK3CD) %>% dplyr::select(-PIK3R1, -PIK3CD)

# Create full matrix
included_pathways <- rownames(matrix_order)

# Initialize a pathways x genes matrix
full_mat <- matrix(0, nrow = length(included_pathways), ncol = length(included_genes),
  dimnames = list(included_pathways, included_genes))

# Copy existing values in
common_genes <- intersect(colnames(matrix_order), included_genes)
if (length(common_genes) > 0) {
  full_mat[, common_genes] <- as.matrix(matrix_order[, common_genes])
}

# Transpose to have genes as rows like your final object and convert to data.frame
full_matrix <- as.data.frame(t(full_mat))

pheatmap(t(full_matrix),  color=heatcolx, scale="row", fontsize = 6,
         cutree_cols=16)

#######################################################
#######################################################
#######################################################

# Associated phenotypes
syndromes_tb <- subset(syndromes, select=-c(Variant.type, Gene,  Direction, Immunodeficiency, 
                                            IUIS_neut, Metabolic_major, Metabolic_minor, 
                                               All_subgroups, Onset_median, exclude, Category))
syndromes_tb[, ][is.na(syndromes_tb[, ])] <- 0 

# Create matrix of column values
syndromes_tb <- syndromes_tb %>%
  pivot_wider(names_from = Syndrome, names_prefix = "Syndrome_", values_from = Syndrome, values_fn = ~ 4, values_fill = 0)  %>%
  pivot_wider(names_from = Subgroups_n2, names_prefix = "Subgroup_", values_from = Subgroups_n2, values_fn = ~1, values_fill = 0) 

syn2 <- syndromes_tb %>% column_to_rownames(var = "Gene_D")
syn2$Syndrome_ <- NULL
syn2$Subgroup_ <- NULL
syn2$Category_ <- NULL
syn2$`Subgroup_Non-inflammasome` <- NULL
syn2$`Subgroup_Other combined immunodeficiency` <- NULL

syn2 <- syn2 %>%  mutate(across(all_of(c("HLH...MAS", "Vasculitis_Panniculitis", "SLE")), 
                                ~ ifelse(. == 1, 1, .)))
syn2 <- syn2 %>%  mutate(across(all_of(c("IBD", "Autoinflammatory")), 
                                ~ ifelse(. == 1, 0.5, .)))
syn2 <- syn2 %>%  mutate(across(all_of(c("Uveitis", "Endocrinopathy", "Psoriasis")),
                                ~ ifelse(. == 1, 0.05, .)))

#######################################################
#######################################################

# Use objects generated in 3.2 expression code

# Create a combined dataset that takes an average if present in both datasets or just the value for one dataset if not present in the other 
# This assumes technical aspects may limit retrieval of some cell types for example neutrophils / mast cells, as oppose to them not being present in tissue 
# Identify which rows (genes) and columns (cell types) have any data in each matrix
top_both2 <- unique(c(top_jia2,top_ed2))

mat1 <- as.matrix(t(df_select_gen))
mat2 <- as.matrix(t(df_select_gen_ed))

# Union of row and column names
all_genes <- unique(c(rownames(mat1), c(rownames(mat2))))

# Function to expand matrix and fill missing entries with 0
expand_with_zeros <- function(mat, rows, cols) {
  out <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
  r <- intersect(rownames(mat), rows); c <- intersect(colnames(mat), cols)
  out[r, c] <- mat[r, c, drop = FALSE]
  out
}

mat1_full <- expand_with_zeros(mat1, all_genes, top_both2)
mat2_full <- expand_with_zeros(mat2, all_genes, top_both2)

gene_in_mat1 <- rowSums(mat1_full != 0) > 0
gene_in_mat2 <- rowSums(mat2_full != 0) > 0

cell_in_mat1 <- colSums(mat1_full != 0) > 0
cell_in_mat2 <- colSums(mat2_full != 0) > 0

# A cell (i, j) is represented only if both that Gene and that Cell Type exist in the dataset
rep1 <- outer(gene_in_mat1, cell_in_mat1, "&")
rep2 <- outer(gene_in_mat2, cell_in_mat2, "&")

# Create the divisor matrix for determining how many datasets to divide the value by
# This will be 2 if the gene/cell combo exists in both, and 1 if only in one
divisor <- rep1 + rep2
divisor[divisor == 0] <- 1  # Prevent division by zero for genes/cells missing in both

# 4. Calculate the conditional average
mat_avg <- (mat1_full + mat2_full) / divisor

# Remove arthropathy genes
arthrop <- c("MMP2", "PRG4", "IDUA", "CCN6")
mat_avg <- mat_avg[!rownames(mat_avg) %in% arthrop,]

#######################################################
#######################################################

# Take average values for both datasets
pheno_s <- scale(syn2)
path_s  <- scale(full_matrix) 
mat_s <- scale(mat_avg)

#######################################################
#######################################################
#######################################################

# Take average values for both datasets
combi <- merge(mat_s, path_s, by=0, all=TRUE)
combi <- combi %>% column_to_rownames("Row.names")
combi <- merge(combi, pheno_s, by=0, all=TRUE)
combi <- combi %>% column_to_rownames("Row.names")

ph <- pheatmap(t(combi), color = heatcol,fontsize = 7, border_color = "grey91", scale = "column",
               cutree_rows = 16, cutree_cols = 20,
               clustering_distance_cols = "correlation", clustering_distance_rows="correlation")

gene_order = rownames(combi)[ph$tree_col$order]
col_order <- colnames(combi)[ ph$tree_row$order ]

mat <- t(combi) # Matrix supplied to pheatmap
mat_ordered <- mat[ph$tree_row$order, ph$tree_col$order] # Reorder according to the heatmap clustering
write.csv(mat_ordered, "heatmap_matrix_ordered.csv")

save.image(file="260913_clustering_objects.RData")

#######################################################
#######################################################
#######################################################

# Meta data
arth <- read.csv(file="CSV/arthritis_csv.csv", header=T, row.names=1)
arth[arth == ""] <- 0
arth2 <- arth[,c(1:3), drop=F]

arth$RF_status <- with(arth, ifelse(RF_pos == 1 & RF_neg == 1, "Both", ifelse(RF_pos == 1 & RF_neg == 0, "RF_pos", ifelse(RF_pos == 0 & RF_neg == 1, "RF_neg", "Neither"))))
arth$ANA_status <- with(arth, ifelse(ANA_pos == 1 & ANA_neg == 1, "Both", ifelse(ANA_pos == 1 & ANA_neg == 0, "ANA_pos", ifelse(ANA_pos == 0 & ANA_neg == 1, "ANA_neg", "Neither"))))
arth$eros_status <- with(arth, ifelse(Erosions == 1 & Non.erosive == 1, "Both", ifelse(Erosions == 1 & Non.erosive == 0, "Erosions", ifelse(Erosions == 0 & Non.erosive == 1, "Non-erosive", "Neither"))))
arth$size <- with(arth, ifelse(Large.joints == 1 & Small.joints == 1, "Both", ifelse(Large.joints == 1 & Small.joints == 0, "Large_joints", ifelse(Large.joints == 0 & Small.joints == 1, "Small_joints", "Neither"))))
arth$Large_joints_only <- ifelse(arth$size == "Large_joints", 1, 0)

arth2 <- arth[,c("ANA_status", "RF_status", "eros_status",  "Axial...spondyloarthritis", "size", "Monoarthritis", "Oligoarthritis", "Polyarthritis"), drop=F]

#######################################################
#######################################################
#######################################################

syndromes_tb <- subset(syndromes, select=-c(Variant.type,  Gene, All_subgroups, Onset_median, exclude, Immunodeficiency, Category, 
                                            IUIS_neut, Metabolic_minor, Metabolic_major, Subgroups_n2))

syn_anno <- syndromes_tb %>% column_to_rownames("Gene_D")
syn_anno <- syn_anno %>% mutate(Syndrome = replace_na(str_trim(na_if(Syndrome, "")), "None"))

# Prepare annotation data frame for pheatmap, convert categorical columns to factors for proper annotation coloring
syn_anno$Syndrome <- as.factor(syn_anno$Syndrome)
syn_anno$Direction <- as.factor(syn_anno$Direction)

# Generate distinct colours for each unique level of the annotation factors
syndrome_colours <- c("Hyper IgM" = "purple3", "Wiskott-Aldrich /-like"  = "darkolivegreen2","LSD" = "grey45", "None"="white", "Hyper IgD"  = "deeppink3",  "IPEX /-like"  = "cyan", "Hypomorphic SCID"  = "deepskyblue3",  
                      "AGS /-like"  = "cornsilk3", "FCAS" = "lightcoral",  "EDA-ID /-like" = "forestgreen", "CGD"= "goldenrod3")

direction_colours <- setNames(c("red", "white")[1:length(levels(syn_anno$Direction))], levels(syn_anno$Direction))

syn_anno$Vasculitis_Panniculitis <- as.numeric(syn_anno$Vasculitis_Panniculitis)
syn_anno$Autoinflammatory <- as.numeric(syn_anno$Autoinflammatory)
syn_anno$IBD <- as.numeric(syn_anno$IBD)

syn_anno2 <- merge(syn_anno, arth2, by="row.names")
syn_anno2 <- syn_anno2 %>% column_to_rownames("Row.names")

anno_col <- c("RF_neg"="grey95", "RF_pos"="purple1", "Both"="purple", Neither="white")
anno_col2 <- c("ANA_neg"="grey95", "ANA_pos"="turquoise4", "Both"="turquoise2", Neither="white")
anno_col3 <- c("Erosions"="navy", "Non-erosive"="lightblue", "Both"="blue", Neither = "white")  
anno_col4 <- c("0"="grey95", "1"="purple3", "NR"="white")
anno_col5 <- c("0"="grey95", "1"="orange", "NR"="white")
anno_col6 <- c("0"="grey95", "1"="skyblue", "NR"="white")
anno_col7 <- c("0"="grey95", "1"="lightgreen", "NR"="white")
anno_col8 <- c("0"="grey95", "1"="red", "NR"="white")

annotation_colours_list <- list(
  Direction = direction_colours,
  Syndrome = syndrome_colours,
  Neutrophils = c("white", "grey"),
  Autoinflammatory = c("white", "plum3"), 
  SLE = c("white", "salmon"), 
  Endocrinopathy = c("white", "forestgreen"), 
  HLH...MAS = c("white", "lightpink2"), 
  Psoriasis = c("white", "cadetblue2"), 
  Uveitis = c("white", "steelblue"),
  IBD = c("white", "tan3"),
  Vasculitis_Panniculitis = c("white", "gold"),
  RF_status = anno_col,
  ANA_status = anno_col2,
  eros_status = anno_col3,
  Axial...spondyloarthritis=anno_col4,
  Monoarthritis = anno_col5,
  Polyarthritis = anno_col6,
  Oligoarthritis = anno_col7)

ph <- pheatmap(t(combi), color = heatcol,fontsize = 7, border_color = "grey91", scale = "column",
               cutree_rows = 16, cutree_cols = 18,
               clustering_distance_cols = "correlation", clustering_distance_rows="correlation",
               annotation_col = arth2, annotation_colors = annotation_colours_list)

#######################################################
#######################################################
#######################################################

both_combi <- mat_avg[gene_order, colnames(mat_avg) %in% top_both2]
cell_order <- gsub(".x", "", col_order)
cell_order <- cell_order[cell_order %in% colnames(mat1_full)]
expr_b <- both_combi[gene_order, cell_order]

column <- colnames(mat_avg)

col_label <- read.csv(file="CSV/monogenic_sheet_7.csv", header=T, stringsAsFactors = FALSE) %>%
      mutate(cell_types = str_squish(as.character(cell_types)))

col_label <- col_label %>%
  mutate(main = str_squish(main), main = case_when(
    main %in% c("Treg", "NK_ILC") ~ "T cell",
    main %in% c("LAMP3_DC","Macrophage_DC") ~ "Macrophage", 
    main %in% c("pDC", "B cell") ~ "pDC/B cell", TRUE ~ main)) 

unique(col_label$main)

# Create colour annotation for main cell types
avg_major <- both_combi %>% as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_types", values_to = "expr") %>%
  mutate(cell_types = str_squish(as.character(cell_types))) %>%
  left_join(col_label, by = "cell_types")

# Summarise expression for the main cell groups
avg_major <- avg_major  %>% group_by(gene, main) %>%                  
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = main, values_from = mean_expr) %>%
  column_to_rownames("gene")

cell_ord <- c("Macrophage","Granulocyte",  "T cell", "pDC/B cell", "Fibroblast","Vascular")

col_label$main <- str_squish(col_label$main)
colnames(avg_major) <- str_squish(colnames(avg_major))
cell_ord <- str_replace_all(cell_ord, "\u00A0", " ")

rev(gene_order)

# rotate dendrogram for clearer visualisation
gene_dend <- c(genes <- c(
  "LACC1", "UNC93B1", "C1QC", "C1QA", "C1QB", "CARD14", "ADA2", "TYMP",
  "NCF1", "NCF2", "ASAH1", "GLA", "NFIL3", "SLC29A3", "ACP5", "ELF4",
  "DNASE2", "SAMHD1", "TREX1", "DOCK11", "IL10RA", "IL10RB", "LYN", "SPI1",
  "PLCG2", "BTK", "SYK", "NLRP1", "NLRP3", "NLRP12", "MEFV", "NLRC4",
  "ALPK1", "NOD2", "ARPC1B", "WAS", "RMRP", "IRF2BP2", "TRNT1", "AIRE",
  "ICOS", "PDCD1", "PRF1", "RAB27A", "UNC13D", "IKZF1", "PSTPIP1", "DEF6",
  "LRBA", "CARMIL2", "CD27", "IFITM5", "IL2RG", "SOCS1", "CTLA4", "FOXP3",
  "DCLRE1C", "PRKDC", "RAG1", "RAG2", "AICDA", "CD40LG", "ATM", "TNFRSF13B",
  "OTULIN", "TNFAIP3", "TNFRSF1A", "SHARPIN", "XIAP", "RIPK1", "RBCK1",
  "RNF31", "DNASE1L3", "EPCAM", "GATA2",
  "HYAL1", "C1S", "C2", "C4A", "C4B", "MVK", "PMVK", "C2orf69", "RNASEH2B",
  "COPA", "SHOC2", "G6PC3", "MAN2B2", "RFXANK", "DDX58", "IFIH1", "TBK1", "TMEM173", 
  "PIK3R1_GOF", "PIK3R1_LOF",
  "PIK3CD_GOF", "PIK3CD_LOF", "STAT1", "STAT3", "NFKB2", "NFKB1", "RELA",
  "IKBKG", "NFKBIA", "PSMB8", "RELB"))

gene_dend
avg_major <- avg_major[gene_dend, cell_ord]


pheatmap(avg_major, main = "", color=heatcol, scale="row", fontsize=7,
         cluster_rows=F, cluster_cols = F, border_color = "grey91", 
         annotation_row = syn_anno, annotation_colors = annotation_colours_list)

#######################################################
#######################################################
#######################################################

col_label <- read.csv(file="CSV/monogenic_sheet_7.csv", header=T, stringsAsFactors = FALSE) %>%
  mutate(cell_types = str_squish(as.character(cell_types)))

unique(col_label$main)

# Create colour annotation for main cell types
avg_major <- both_combi %>% as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_types", values_to = "expr") %>%
  mutate(cell_types = str_squish(as.character(cell_types))) %>%
  left_join(col_label, by = "cell_types")

# Summarise expression for the main cell groups
avg_major <- avg_major  %>% group_by(gene, main) %>%                  
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = main, values_from = mean_expr) %>%
  column_to_rownames("gene")

col_label$main <- str_squish(col_label$main)
colnames(avg_major) <- str_squish(colnames(avg_major))
cell_ord <- str_replace_all(cell_ord, "\u00A0", " ")

cell_ord <- str_squish(c("Granulocyte", "Macrophage_DC", "LAMP3_DC", "T cell", "Treg", "NK_ILC", "pDC", "B cell", "Vascular", "Fibroblast"))

avg_major <- avg_major[gene_order, cell_ord]

pheatmap(avg_major, main = "", color=heatcol, scale="row", fontsize=7,
         cluster_rows=F, cluster_cols=F, annotation_row = syn_anno2, annotation_colors = annotation_colours_list)

#######################################################
#######################################################

#######################################################

path_order <- col_order[col_order %in% colnames(path_s)]
plot_matrix <- full_matrix[gene_dend, path_order]

cl <- cutree(ph$tree_col, k = 20)
check <- as.data.frame(cl)
check

overlap <- sapply(1:20, \(x) colSums(plot_matrix[names(cl)[cl == x], , drop = FALSE] > 0))
overlap <- as.data.frame(overlap)
colnames(overlap) <- paste0("Cluster_", 1:20)

top_pathways <- do.call(rbind, lapply(names(overlap), \(x)
                                      head(data.frame(cluster=x, pathway=rownames(overlap)[order(overlap[[x]], decreasing=TRUE)],
                                                      overlap=sort(overlap[[x]], decreasing=TRUE)), 4)))

#######################################################
#######################################################
#######################################################

genes_dend <- gene_order 

# Treatment response
rx <- read.csv(file="CSV/monogenic_sheet_8.csv", header=T)
rx <- rx[rx$Gene %in% genes_dend,]
rx <- rx %>% dplyr::select(-c(Ref))

treatments <- c("MTX","TNF","JAK","IL1","IL6","CTLA4","CD20","HSCT")

rx_long <- rx %>% as.data.frame() %>% 
  pivot_longer(cols = -Gene,
    names_to = c("treatment", "kind"),
    names_sep = "\\.",
    values_to = "value") %>%
  filter(treatment %in% treatments) %>%           
  pivot_wider(names_from = kind, values_from = value)

#Proportions of response
rx_summary <- rx_long %>% mutate(total = as.numeric(total),
    response = as.numeric(response),
    prop = if_else(is.na(total) | total == 0, NA_real_, response / total))

# Gene order
ordered_genes <- c(intersect(genes_dend, unique(rx_summary$Gene)),
                   setdiff(unique(rx_summary$Gene), genes_dend))

# Categorise proportions
rx_summary <- rx_summary %>%
  mutate(prop_cat = case_when(
      prop <= 0.2 ~ "0-20%",
      prop <= 0.4 ~ "21-40%",
      prop <= 0.6 ~ "41-60%",
      prop <= 1   ~ "61-100%",
      TRUE        ~ NA_character_))

# Specify colours
prop_cols <- c("0-20%"   = "#1F78B4",  "21-40%"  = "lightblue",  "41-60%" = "darkorange", "61-100%" = "#E31A1C")

# Specify order
rx_summary$Gene <- factor(rx_summary$Gene, levels = rev(gene_dend))
treat_order <- c("IL1", "TNF", "IL6", "MTX", "JAK", "CTLA4", "CD20","HSCT")
rx_summary$treatment <- factor(rx_summary$treatment, levels=treat_order)

rx_summary <- rx_summary[rx_summary$treatment %in% c("IL1", "TNF", "IL6", "MTX","HSCT"),]

# Plot
ggplot(rx_summary, aes(x = treatment, y = Gene)) + geom_point(aes(size = total, color = prop_cat), alpha = 0.9) +
  scale_color_manual(name = "Proportion responding", values = prop_cols, drop = FALSE) +
  scale_size_continuous(name = "Total treated", range = c(1.5, 8)) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=7),
    panel.grid.major.y = element_blank())

#######################################################
#######################################################
#######################################################

# Calculate contributions to clustering
genes <-  rownames(both_combi)
expr2 <- mat_s[genes,,drop=FALSE]
path2 <- path_s[genes,,drop=FALSE]
pheno2 <- pheno_s[genes,,drop=FALSE]
combi <- combi[genes,,drop=FALSE]

d_all <- as.dist(1-cor(t(combi), use="pairwise.complete.obs"))
d_expr <- as.dist(1-cor(t(expr2), use="pairwise.complete.obs"))
d_path <- as.dist(1-cor(t(path2), use="pairwise.complete.obs"))
d_pheno <- as.dist(1-cor(t(pheno2), use="pairwise.complete.obs"))

cor(d_all, d_expr, method="spearman")
cor(d_all, d_path, method="spearman")
cor(d_all, d_pheno, method="spearman")

combi_no_expr <- cbind(path2, pheno2)
combi_no_path <- cbind(expr2, pheno2)
combi_no_pheno <- cbind(expr2, path2)

d_no_expr <- as.dist(1-cor(t(combi_no_expr), use="pairwise.complete.obs"))
d_no_path <- as.dist(1-cor(t(combi_no_path), use="pairwise.complete.obs"))
d_no_pheno <- as.dist(1-cor(t(combi_no_pheno), use="pairwise.complete.obs"))

cor(d_all, d_no_expr, method="spearman")
cor(d_all, d_no_path, method="spearman")
cor(d_all, d_no_pheno, method="spearman")

#######################################################
#######################################################
#######################################################

library(mclust)

k <- 20

# Full clustering
d_full <- as.dist(1 - cor(t(combi), use="pairwise.complete.obs"))
cl_full <- cutree(hclust(d_full, method="complete"), k=k)

# Leave one feature out
loo_results <- lapply(seq_len(ncol(combi)), function(i) {
  x <- combi[, -i, drop=FALSE]
  d <- as.dist(1 - cor(t(x), use="pairwise.complete.obs"))
  cl <- cutree(hclust(d, method="complete"), k=k)
  data.frame(feature=colnames(combi)[i],
             ARI=adjustedRandIndex(cl_full, cl),
             distance_cor=cor(d_full, d, method="spearman"))
})

loo_results <- do.call(rbind, loo_results)

summary(loo_results$ARI)
summary(loo_results$distance_cor)
