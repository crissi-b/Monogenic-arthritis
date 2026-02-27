
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


col_fun <- c("turquoise4", "turquoise3", "#AAEBEE", "#FFFFFF", "#D8C1B0", "#8B4513", "brown4")
heatcol2 <- colorRampPalette(c("blue", "dodgerblue", "lightblue1","white","white", "darkorange","red"))(10)
heatcol3 <- rev(brewer.pal(8, "RdYlBu"))
heatcol <- rev(brewer.pal(8, "RdBu"))

#######################################################

# Load scRNA-seq expression data for adults and children, pathway outputs and syndromic features from monogenic_sheet_1


#######################################################

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
included_genes <- colnames(both_combi)
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

#######################################################

# Associated phenotypes
syndromes <- read.csv(file="monogenic_sheet_1.csv", header=F, blank.lines.skip = TRUE)
syndromes <- syndromes[3:nrow(syndromes),]
colnames(syndromes) <- syndromes[1, ]
syndromes <- syndromes[-1, ]
syndromes <- syndromes %>% filter(exclude == "N")

syndromes_tb <- subset(syndromes, select=-c(Variant.type, age_group, Gene,  Direction, age_group, 
                                               X, All_subgroups, Onset_median, exclude, Category))
syndromes_tb[, ][is.na(syndromes_tb[, ])] <- 0 

# Create matrix of column values
syndromes_tb <- syndromes_tb %>%
  pivot_wider(names_from = Syndrome, names_prefix = "Syndrome_", values_from = Syndrome, values_fn = ~ 3, values_fill = 0)  %>%
  pivot_wider(names_from = Subgroups_n2, names_prefix = "Subgroup_", values_from = Subgroups_n2, values_fn = ~0.5, values_fill = 0) 

syn_long <- syndromes_tb %>% column_to_rownames(var = "Gene_D")

syn2 <- syn_long[colnames(both_combi),]
syn2$Syndrome_ <- NULL
syn2$Subgroup_ <- NULL

syn2$`Subgroup_Non-inflammasome` <- NULL
syn2$`Subgroup_Other combined immunodeficiency` <- NULL

syn2 <- syn2 %>%  mutate(across(all_of("Neutrophils"), ~ ifelse(. == 1, 2, .)))
syn2 <- syn2 %>%  mutate(across(all_of(c("Immunodeficiency", "Autoinflammatory", "SLE", 
                                         "Endocrinopathy", "HLH...MAS", "Vasculitis_Panniculitis")), ~ ifelse(. == 1, 0.5, .)))
syn2 <- syn2 %>%  mutate(across(all_of(c("Psoriasis",
                                         "Uveitis", "IBD")), ~ ifelse(. == 1, 0.1, .)))

#######################################################

# Create a combined dataset that takes an average if present in both datasets or just the value for one dataset if not present in the other 
# This assumes technical aspects may limit retrieval of some cell types for example neutrophils / mast cells, as oppose to them not being present in tissue 
# Use objects generated in 3.2 expression code

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

# Filter for highest expressing cell types
both_combi <- mat_avg[top_both, ]

pheatmap(both_combi,  color=heatcolx, scale="column")

data <- t(both_combi)
max_col_indices <- apply(data, 1, which.max)
max_col <- colnames(data)[max_col_indices]
top_cells <- data.frame(Gene_D=colnames(both_combi), Cells = max_col)

#######################################################
#######################################################
#######################################################

expr_j <- scale(t(mat1_full), center=TRUE, scale=TRUE)
expr_r <- scale(t(mat2_full), center=TRUE, scale=TRUE)

pheno_s <- scale(syn2)
path_s  <- scale(full_matrix) 

#######################################################
#######################################################
#######################################################

combi <- merge(expr_j, expr_r, by=0, all=TRUE)
combi <- combi %>% column_to_rownames("Row.names")

combi <- merge(expr_s, path_s, by=0, all=TRUE)
combi <- combi %>% column_to_rownames("Row.names")

combi <- merge(combi, pheno_s, by=0, all=TRUE)
combi <- combi %>% column_to_rownames("Row.names")

ph <- pheatmap(t(combi), main = "", color=heatcol, 
               fontsize = 7, border_color = "grey91", scale="row",
               cutree_rows = 15, cutree_cols = 16, clustering_distance_rows = "correlation",
               clustering_distance_cols="correlation")

ph

# which are genes (rownames(combi)). So gene order is:
gene_order = rownames(combi)[ ph$tree_col$order ]
gene_order

col_order <- colnames(combi)[ ph$tree_row$order ]
combi <- combi[gene_order, col_order]

# Rotate around dendrogram for ease of interpretation
genes_dend <- gene_order[c(2:6,1,20:11,7:10,21:26,33:29,27:28,34:57,66:72,58:65,96:99,100:106,81:84,85:95,73:80)]
genes_dend

#######################################################

cell_order <- gsub(".x", "", col_order)
cell_order <- cell_order[cell_order %in% rownames(mat1_full)]
expr_b <- both_combi[cell_order, gene_order] %>% t()

column <- colnames(expr_s)

col_label <- read.csv(file="monogenic_sheet_7.csv", header=T, stringsAsFactors = FALSE) %>%
      mutate(cell_types = str_squish(as.character(cell_types)))

# Create colour annotation for main cell types
avg_major <- t(both_combi) %>% as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_types", values_to = "expr") %>%
  mutate(cell_types = str_squish(as.character(cell_types))) %>%
  left_join(col_label, by = "cell_types")

# Summarise expression for the main cell groups
avg_major <- avg_major  %>% group_by(gene, main) %>%                  
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = main, values_from = mean_expr) %>%
  column_to_rownames("gene")

pheatmap(avg_major, main = "", color=heatcol, scale="column", cutree_cols = 8, cutree_rows = 12,
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")

# Order manually by previous order
avg_major <- avg_major[genes_dend2,]

pheatmap(avg_major, main = "", 
         color=heatcol, scale="column", cutree_cols = 8, cutree_rows = 12)

avg_major <- avg_major[genes_dend,]

pheatmap(avg_major, main = "", 
         color=heatcol, scale="row", cutree_cols = 8, cutree_rows = 12,
         cluster_rows=F, 
         clustering_distance_cols = "correlation")


#######################################################
#######################################################

syndromes_tb <- subset(syndromes, select=-c(Variant.type,  Gene, 
                                               Neutrophils, All_subgroups, Onset_median, exclude))

syn_anno <- syndromes_tb %>% column_to_rownames("Gene_D")
syn_anno <- syn_anno %>% mutate(Syndrome = replace_na(str_trim(na_if(Syndrome, "")), "None"))

# Prepare annotation data frame for pheatmap, convert categorical columns to factors for proper annotation coloring
syn_anno$Syndrome <- as.factor(syn_anno$Syndrome)
syn_anno$Direction <- as.factor(syn_anno$Direction)
syn_anno$age_group <- as.factor(syn_anno$age_group)

# Generate distinct colours for each unique level of the annotation factors
syndrome_colours <- c("Hyper IgM" = "purple3", "Wiskott-Aldrich /-like"  = "darkolivegreen2","LSD" = "grey45", "None"="white",
                      "Hyper IgD"  = "deeppink3",  "IPEX /-like"  = "cyan", "Hypomorphic SCID"  = "deepskyblue3",  
                      "AGS /-like"  = "cornsilk3", "FCAS" = "lightcoral",  "EDA-ID /-like" = "forestgreen", "CGD"= "goldenrod3")


direction_colours <- setNames(c("red", "white")[1:length(levels(syn_anno$Direction))], levels(syn_anno$Direction))

syn_anno$Vasculitis_Panniculitis <- as.numeric(syn_anno$Vasculitis_Panniculitis)
syn_anno$Autoinflammatory <- as.numeric(syn_anno$Autoinflammatory)
syn_anno$IBD <- as.numeric(syn_anno$IBD)

# Combine all annotation colours into a list for pheatmap
annotation_colours_list <- list(
  Direction = direction_colours,
  Syndrome = syndrome_colours,
  Autoinflammatory = c("white", "plum3"), 
  IBD = c("white", "tan3"),
  Vasculitis_Panniculitis = c("white", "gold"))

syn2 <- syn_anno[,rev(c("Syndrome", "Direction", 
                        "Autoinflammatory",  "IBD", "Vasculitis_Panniculitis"))]

pheatmap(avg_major, color=heatcol, fontsize = 7, border_color = "grey91", 
         scale="row", cluster_rows = F,
         annotation_row = syn2, annotation_colors = annotation_colours_list) 

colnames(avg_major)

avg_matrix <- avg_major[genes_dend,c(8,7,6,10,3,5,4,1,2,9)]
pheatmap(avg_matrix, color=heatcol, fontsize = 6, border_color = "grey91", 
         scale="row", cluster_cols=F, cluster_rows = F,
         annotation_row = syn2, annotation_colors = annotation_colours_list) 

#######################################################

path_order <- col_order[col_order %in% colnames(path_s)]
plot_matrix <- full_matrix[genes_dend, path_order]

pheatmap(plot_matrix, color=c( "white", "orange"), fontsize = 5, border_color = "grey91", 
         cluster_rows=F, clustering_distance_cols = "correlation") 

pheatmap(t(plot_matrix), color=c( "white", "orange"), fontsize = 8, border_color = "grey91", 
         cluster_rows=F, cluster_cols=F) 

#######################################################

# Treatment response
rx <- read.csv(file="monogenic_sheet_8.csv", header=T)
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
rx_summary$Gene <- factor(rx_summary$Gene, levels = rev(genes_dend))
treat_order <- c("CTLA4", "MTX","TNF","IL1","IL6","JAK","CD20","HSCT")
rx_summary$treatment <- factor(rx_summary$treatment, levels=treat_order)

# Plot
ggplot(rx_summary, aes(x = treatment, y = Gene)) + geom_point(aes(size = total, color = prop_cat), alpha = 0.9) +
  scale_color_manual(name = "Proportion responding", values = prop_cols, drop = FALSE) +
  scale_size_continuous(name = "Total treated", range = c(1.5, 8)) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=7),
    panel.grid.major.y = element_blank())


