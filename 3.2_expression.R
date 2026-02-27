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
library(SeuratObject)
library(pheatmap)

heatcol <- rev(brewer.pal(8, "RdBu"))
heatcolx <- c("#2166AC", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7","#F4A582", "#D6604D" ,"#B2182B","#B2182B")

#######################################################
#######################################################
#######################################################

# Plot PBMC expression of monogenic genes
syndromes <- read.csv(file="monogenic_sheet_1.csv", header=F, blank.lines.skip = TRUE)
syndromes <- syndromes[3:nrow(syndromes),]
colnames(syndromes) <- syndromes[1, ]
syndromes <- syndromes[-1, ]

excluded_genes <- syndromes %>% filter(exclude == "Y") %>% pull(Gene)
gen <- syndromes %>% filter(exclude == "N") %>% pull(Gene)
  
# Check expression in PBMC dataset from Human Protein Atlas, available here: https://www.proteinatlas.org/download/tsv/rna_immune_cell_monaco.tsv.zip
immune_full <- read.delim(file="Datasets/rna_immune_cell_monaco.tsv")

# Find genes in dataset
immune <- immune_full[immune_full$Gene.name %in% gen,]

# Find missing genes
setdiff(gen, unique(immune$Gene.name))

# Add gene aliases
gen <- c(gen, "STING1", "RIGI", "GBA1")
setdiff(unique(immune$Gene.name),gen)

# Reformat
immune <- immune_full[immune_full$Gene.name %in% gen,]
immune1 <- reshape2::dcast(immune, formula = Gene.name ~ Immune.cell, value.var= "TPM")
rownames(immune1) <- immune1$Gene.name
immune1 <- immune1 %>% as.data.frame() %>% dplyr::select(-c(Gene.name, "total PBMC"))

# Amend gene names to align with other datasets
rownames(immune1)[rownames(immune1) == "STING1"] <- "TMEM173"
rownames(immune1)[rownames(immune1) == "GBA1"] <- "GBA"
rownames(immune1)[rownames(immune1) == "RIGI"] <- "DDX58"

immune1 <- immune1 %>% as.matrix() %>% t()
absent <- immune1[, colSums(immune1)==0] # Check zero columns
immune1 <- immune1[, colSums(immune1 != 0) > 0]

# Reformat names for plotting
rownames(immune1) <- str_to_sentence(rownames(immune1))
rownames(immune1) <- gsub("cd4", "CD4+", rownames(immune1))
rownames(immune1) <- gsub("cd8", "CD8+", rownames(immune1))
rownames(immune1) <- gsub("t-cell", "T cell", rownames(immune1))
rownames(immune1) <- gsub("b-cell", "B cell", rownames(immune1))
rownames(immune1) <- gsub("dc", "DC", rownames(immune1))
rownames(immune1) <- gsub("th1", "Th1", rownames(immune1))

pheatmap(immune1, scale="column", col= heatcol, 
         main="Immune Cell Expression (Human Protein Atlas)", 
         fontsize=8, border_color = "grey", cutree_cols = 5, clustering_distance_cols = "correlation")

# Select cell types which a gene is maximally expressed in for more condensed visualisation
data <- immune1
max_row_indices <- apply(data, 2, which.max)
max_rows <- unique(rownames(data)[max_row_indices])
immune2 <- data[max_rows, ]
top_immune <- rownames(immune2)
max_immune <- as.data.frame(rownames(data)[max_row_indices])
table(max_immune)

cell.type <- as.data.frame(rownames(data)[max_row_indices])
table(cell.type$`rownames(data)[max_row_indices]`)

immune2 <- immune2 %>% as.data.frame() %>%
  mutate(`PIK3R1_LOF` = PIK3R1, `PIK3R1_GOF` = PIK3R1, 
         `PIK3CD_LOF` = PIK3CD, `PIK3CD_GOF` = PIK3CD) %>% dplyr::select(-PIK3R1, -PIK3CD)

###############################################################

# Add annotation of immune category, and IUIS neutrophil/T cell listed effects
syn_anno <- syndromes %>% filter(exclude == "N")
syn_anno <- subset(syn_anno, select=c(Gene_D, Neutrophils))
syn_anno[, ][is.na(syn_anno[, ])] <- 0 # Convert NA values to 0

# Replace missing values with None
syn_anno <- syn_anno %>% column_to_rownames("Gene_D")
syn_anno$Neutrophils <- as.numeric(syn_anno$Neutrophils)

annotation_colours_list <- list("Neutrophils"=c("white", "forestgreen"))
  
syn2 <- syn_anno[colnames(immune2), , drop = FALSE]

pheatmap(as.matrix(immune2), scale="column", col= heatcol, 
         main="PBMC Expression", fontsize=7, border_color = "grey",  cutree_cols = 6,
         annotation_col = syn2, annotation_colors = annotation_colours_list, na_col = "white",
         clustering_distance_cols = "correlation")

#######################################################
#######################################################
#######################################################

# Extract expression of IUIS neutrophil labelled cells to examine gene expression in neutrophils vs other cells
neut <- syn_anno %>% dplyr::select("Neutrophils") %>% rownames_to_column("Gene_D")
immune3 <- t(immune2) %>% as.data.frame() %>% rownames_to_column("Gene_D")
immune3 <- merge(immune3, neut, by="Gene_D") 

immune_long <- immune3 %>%
  pivot_longer(cols = -c(Gene_D, Neutrophils), names_to = "CellType", values_to = "Expression")

# Compare expression in 'Neutrophil' column vs all other cell types
immune_combine <- immune_long %>%
  mutate(IsNeutrophil = ifelse(CellType == "Neutrophil", "Neutrophil", "Other")) %>%
  group_by(Neutrophils, IsNeutrophil) %>% summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

ggplot(immune_combine, aes(x = factor(Neutrophils), y = mean_expr, fill = IsNeutrophil)) +
  geom_col(position = "dodge") +
  labs(x = "", y = "Mean expression", title = "Average expression in Neutrophil vs other cell types") +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank(),
                                   panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

# Are Neutrophil-expressed genes higher in the Neutrophil column?
# Compare expression in 'Neutrophil' column vs all other cell types

expr_test <- immune3 %>%
  dplyr::select(Gene_D, Neutrophil, Neutrophils)

# Wilcoxon rank-sum test comparing expression in the Neutrophil column 
# between genes with Neutrophils == 1 vs == 0
wilcox.test(Neutrophil ~ Neutrophils, data = expr_test)

#######################################################
#######################################################
#######################################################

# Visualising Synovial expression of genes- additional code provided of initial processing
# Visualise aggregate gene counts in JIA dataset, coloured by max tissue expression

# Load JIA tissue dataset
load(file="JIA-tissue-annotated.rds")

# First format the tissue dataset to keep cell types uniform 
# Remove cell types with low counts in tissue
tissue.small <- subset(tissue, percent.mt < 6 & nCount_RNA > 500)

# Remove cell types with < 120 cells in tissue
check <- as.data.frame(table(tissue.small$RA_clusters))
low.counts <- check[check$Freq<120,] %>% pull(Var1)
low.counts
Idents(tissue.small) <- "RA_clusters"
tissue.small <- subset(tissue.small, idents=low.counts, invert=T)

#######################################################

gen <- syndromes %>% filter(exclude == "N") %>% pull(Gene)

# Pseudobulk counts- using cell names from label tranGene# Pseudobulk counts- using cell names from label transfer of AMP2 RA dataset (Zhang et al 2023)
jia_counts <- AggregateExpression(tissue.small, group.by = "RA_clusters",
                                  assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA %>% as.matrix()

# Format row names
colnames(jia_counts) <- gsub("g ", "", colnames(jia_counts))
jia_total <- jia_counts[rownames(jia_counts) %in% gen,]
jia_total <- as.data.frame(rowSums(jia_total))

#######################################################

# Repeat for adult RA dataset (Edalat et al 2024)

load(file="synov_annotated.RData")

# Create unique ID of each patient
synov$patient <- sub("\\..*", "", colnames(synov))
synov$patient <- sub("Syn_Bio_", "", synov$patient)
unique(synov$patient)

# Correct names for analysis
Idents(synov) <- "RA_clusters"
synov <- RenameIdents(synov, "E-5 SERPINF1+ capillary"="E-5: SERPINF1+ capillary", "Neutrophils"="M-15: Neutrophils")
synov$RA_clusters <- synov@active.ident

# QC filtering
synov.low <- subset(synov, percent.mt < 25 & nCount_RNA > 500)

# Remove cell types with low counts in tissue
check <- as.data.frame(table(synov.low$RA_clusters))
low.counts <- check[check$Freq<120,] %>% pull(Var1)
low.counts
Idents(synov.low) <- "RA_clusters"
synov.low <- subset(synov.low, idents=low.counts, invert=T)

# Calculate aggregate counts of genes
ra_counts <- AggregateExpression(synov.low, group.by = "RA_clusters",
                                 assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA %>% as.matrix()

colnames(ra_counts) <- gsub("g ", "", colnames(ra_counts))
ra_total <- ra_counts[rownames(ra_counts) %in% gen,]
ra_total <- as.data.frame(rowSums(ra_total))

setdiff(gen,rownames(ra_total))

# Combine JIA and RA datasets for overall total
total_synov <- merge(jia_total, ra_total, by=0, all=T)
colnames(total_synov) <- c("gene", "JIA", "RA")
total_synov$RA[is.na(total_synov$RA)] <- 0 
total_synov <- total_synov %>% mutate(total_count = JIA + RA)

#######################################################

# Extract which tissue has highest expression of each gene for plotting
# Gene list amended for different gene names
hpa <- read.delim(file="Datasets/rna_tissue_consensus.tsv")

# Extract rows with monogenic gene names
gen_hpa <- c(gen, "GBA1", "RIGI", "STING1")
hpa_gen <- hpa[hpa$Gene.name %in% gen_hpa,]
setdiff(gen, unique(hpa_gen$Gene.name)) # Find missing genes # Not found: "DDX58"   "IFIH1 "  "RMRP"  

# Ensure excluded genes not included
hpa_gen <- hpa_gen[!(hpa_gen$Gene.name %in% excluded_genes),]

# Reformat
hpa_gen1 <- reshape2::dcast(hpa_gen, formula = Gene.name ~ Tissue, value.var= "nTPM")
rownames(hpa_gen1) <- hpa_gen1$Gene.name
hpa_gen1 <- hpa_gen1 %>% dplyr::select(-Gene.name)
hpa_gen1 <- hpa_gen1 %>% as.matrix() %>% t()

# Create dataframe with maximum value for any gene
data <- hpa_gen1
max_tissue <- apply(hpa_gen1, 2, function(x) {
  rownames(hpa_gen1)[which.max(x)] })

max.tissue <- data.frame(gene = names(max_tissue), max_tissue = unname(max_tissue))

# Correct gene names
max.tissue$gene[max.tissue$gene=="RIGI"] <- "DDX58"
max.tissue$gene[max.tissue$gene=="STING1"] <- "TMEM173"
max.tissue$gene[max.tissue$gene=="GBA1"] <- "GBA"

#######################################################

# Merge expression with max_tissue info
all_total <- merge(total_synov, max.tissue, by = "gene")

# Amalgamate tissue types for plotting
test <- all_total %>%
  mutate(max_tissue = case_when(
    max_tissue %in% c("colon", "stomach", "duodenum", "small intestine") ~ "gut",
    max_tissue %in% c("pituitary gland", "midbrain", "choroid plexus", "parathyroid gland", "retina") ~ "brain",
    max_tissue %in% c("epididymis", "testis", "urinary bladder", "seminal vesicle", "placenta") ~ "genitourinary",
    max_tissue %in% c("heart muscle", "skeletal muscle", "tongue") ~ "muscle",
    max_tissue %in% c("liver", "gallbladder", "pancreas") ~ "biliary system",
    max_tissue %in% c("tonsil", "appendix", "lymph node") ~ "secondary lymphoid organs",
    max_tissue %in% c("adipose tissue", "breast") ~ "adipose/breast",
    TRUE ~ max_tissue))

palette <- c("#1B9E77", "#E5C494" , "#A6CEE3" ,"#E41A1C" ,"#FFFFB3","#984EA3", "#D9D9D9" ,
             "#CCEBC5", "cornsilk" ,  "#B3DE69" ,"#FFD92F" , "#FCCDE5","#1F78B4")

test$max_tissue <- str_to_sentence(test$max_tissue)
test <- test[!(test$gene %in% excluded_genes),]

# Plot
ggplot(test, aes(x = reorder(gene,total_count), y = total_count, fill = max_tissue)) + geom_bar(stat = "identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "Total counts per gene in synovium",
       x = "",y = "Total counts (log10)", fill = "Max tissue") + 
  scale_fill_manual(values=palette) + scale_y_continuous(trans = 'log10') 

#######################################################
#######################################################
#######################################################

# Plot specific expression for genes in JIA synovial tissue using dreamlet package
# Create SCE object for analysis and pseudobulk
sce <- as.SingleCellExperiment(tissue.small)
pb <- aggregateToPseudoBulk(sce, assay = "counts", cluster_id = "RA_clusters", sample_id = "patient")

# Compute cell type specificity of each gene
df <- cellTypeSpecificity(pb)
genes_present <- intersect(gen, rownames(tissue))
genes_present <- genes_present[!(genes_present %in% excluded_genes)]

df_select <- df[genes_present,] 
df_select <- df_select %>% as.data.frame()
df_select$totalCPM <- NULL

# Reformat names using cluster names in single cell dataset
df_prefix <- colnames(df_select) %>% as.data.frame()
df_prefix$prefix <- df_prefix$. %>%
  str_extract("^.*?(?=\\.\\.)") %>% 
  str_replace("\\.", "-")     
clean.names <- unique(tissue$RA_clusters) %>% as.data.frame()
clean.names$prefix <- clean.names$. %>% str_extract("^.*?(?=:)")
format.names <- merge(df_prefix, clean.names, by="prefix")
rename_vec <- setNames(format.names$..y, format.names$..x)
colnames(df_select) <- rename_vec[colnames(df_select)]

setdiff(gen, rownames(df_select))
pheatmap(t(df_select), main = "JIA pseudo-bulk", color=heatcolx, 
         fontsize = 7, scale="column", border_color = "grey91")

#######################################################

# Extract highest expressing cell type per gene and only visualise these cell types for condensed visualisation
data <- df_select
max_col_indices <- apply(data, 1, which.max)
max_col <- unique(colnames(data)[max_col_indices])
filtered_jia <- data[, max_col]
top_jia <- colnames(filtered_jia)

pheatmap(t(filtered_jia), main = "JIA specificity scores", color=heatcol, 
         fontsize = 7, border_color = "grey91", scale="column")

#######################################################
#######################################################
#######################################################

# Repeat steps for adult RA (Edalat et al 2024) dataset to come up with an order for both

# Create more specific global annotation
Idents(synov.low) <- "global"
synov.low <- RenameIdents(synov.low, "Cycling cells"="T cells", "Cycling myeloid cells"="Myeloid cells","Cycling stromal cells"= "Stromal cells",
                          "Mast cells"= "Granulocytes","Neutrophils"= "Granulocytes",
                          "Pericytes"="Stromal cells")
synov.low@active.ident -> synov.low$global_detail
table(synov.low$global_detail)

# Create SCE object
sce.ed <- as.SingleCellExperiment(synov.low)

pb.ed <- aggregateToPseudoBulk(sce.ed, assay = "counts", cluster_id = "RA_clusters", sample_id = "patient")

# Compute cell type specificity of each gene
df.ed <- cellTypeSpecificity(pb.ed)

genes_present <- intersect(gen, rownames(synov))
genes_present <- genes_present[!(genes_present %in% c(excluded_genes))]

df_select.ed <- df.ed[genes_present,] 
df_select.ed <- df_select.ed %>% as.data.frame()
df_select.ed$totalCPM <- NULL

# Reformatting names
df_prefix.ed <- colnames(df_select.ed) %>% as.data.frame()
df_prefix.ed$prefix <- df_prefix.ed$. %>%
  str_extract("^.*?(?=\\.\\.)") %>% 
  str_replace("\\.", "-")     
clean.names.ed <- unique(synov$RA_clusters) %>% as.data.frame()
clean.names.ed$prefix <- clean.names.ed$. %>% str_extract("^.*?(?=:)")
format.names.ed <- merge(df_prefix.ed, clean.names.ed, by="prefix")
rename_veced <- setNames(format.names.ed$..y, format.names.ed$..x)
colnames(df_select.ed) <- rename_veced[colnames(df_select.ed)]

setdiff(gen, rownames(df_select.ed))
pheatmap(t(df_select.ed), main = "RA specificity scores", color=heatcolx, 
         fontsize = 7, scale="column", border_color = "grey91")

#######################################################

# Top cell types only
data <- df_select.ed
max_col_indices <- apply(data, 1, which.max)
max_col <- unique(colnames(data)[max_col_indices])
filtered_ed <- data[, max_col]
top_ed <- colnames(filtered_ed)

#######################################################

# Add GOF / LOF columns
df_select <- t(df_select) %>% as.data.frame() %>%
  mutate(`PIK3R1_LOF` = PIK3R1, `PIK3R1_GOF` = PIK3R1, 
         `PIK3CD_LOF` = PIK3CD, `PIK3CD_GOF` = PIK3CD) %>% 
  dplyr::select(-PIK3R1, -PIK3CD) 

df_select.ed <- t(df_select.ed) %>% as.data.frame() %>%
  mutate(`PIK3R1_LOF` = PIK3R1, `PIK3R1_GOF` = PIK3R1, 
         `PIK3CD_LOF` = PIK3CD, `PIK3CD_GOF` = PIK3CD) %>% 
  dplyr::select(-PIK3R1, -PIK3CD) 

# Find a combined order
# Top cell types across both datasets
top_both <- unique(c(top_jia,top_ed))

mat1 <- as.matrix(df_select) # Set as matrices
mat2 <- as.matrix(df_select.ed)

# Union of row and column names
all_rows <- union(rownames(mat1), rownames(mat2))
all_cols <- union(colnames(mat1), colnames(mat2))

# Function to expand matrix with zeros
expand_with_zeros <- function(mat, rows, cols) {
  out <- matrix(0, nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

# Apply to both
mat1_full <- expand_with_zeros(mat1, all_rows, all_cols)
mat2_full <- expand_with_zeros(mat2, all_rows, all_cols)

# Take average values for both datasets
mat_avg <- (mat1_full + mat2_full) / 2

# Filter for highest expressing cell types
mat_top <- mat_avg[top_both, ]

ch <- pheatmap(mat_top, color=heatcolx,  
         clustering_distance_rows="manhattan", clustering_distance_cols="manhattan",
         fontsize = 7, scale="column", border_color = "grey91")

# cluster rows
cell_order <- rownames(mat_top)[ch$tree_row$order]
gene_order <- colnames(mat_top)[ch$tree_col$order]

cell_dend <- cell_order[c(1:10,20:21,14,11,15,13,22:29,12,30:32,34,33,35:39,40,16:19,45:50,41:44)]
gene_dend <- gene_order[c(1:65, 87,66:70,93:94,71:86,88:92,95:106)]


#######################################################
#######################################################

# Correct order
gene_dend2 <- gene_order[c(1:60, 62, 93, 94, 61, 63:65, 87, 66:70, 97, 98, 71:82, 83, 95, 96, 91, 92, 99, 84:86, 88:90, 100:106)]
cell_dend2 <- cell_order[c(1:10,20:21,14,11,15,13,22:29,12,30:32,34,33,35:39,40,16:19,41:44,45:50)]

#######################################################
#######################################################

# Repeat same process for arthropathy genes: "MMP2", "PRG4", "IDUA", "CCN6"

#######################################################
#######################################################

# Set order of JIA and RA datasets
jia2 <- mat1_full[cell_dend2, gene_dend2] %>% as.data.frame()
ed2 <- mat2_full[cell_dend2, gene_dend2] %>% as.data.frame()

# Plot with direction, neutrophils and onset age annotations
# Ensure LOF and GOF versions of gene defects are present and remove original columns
# Remove excluded genes
syn_brief <- syndromes %>% filter(exclude == "N")

# Cut ages into categories
syn_brief <- syn_brief %>% mutate(
  age_group = cut(Onset_median, breaks = c(-Inf, 2.9, 6.9, 16.9,  Inf), labels = c("0-2", "3-6", "7-16", ">17")))

syn_brief <- subset(syn_brief, select=c(Gene_D, age_group, Neutrophils, Direction))
syn_brief[, ][is.na(syn_brief[, ])] <- 0 # Convert NA values to 0
syn_anno <- syn_brief %>% column_to_rownames("Gene_D")
syn_anno <- syn_anno[colnames(jia2), , drop=FALSE] %>% as.data.frame()

# Prepare annotation data frame for pheatmap, convert categorical columns to factors for proper annotation coloring
age_colours <- c("0-2" = "blue", "3-6" = "lightblue1", "7-16" = "lemonchiffon", ">17"= "cornsilk3")

# Add cell type colours for cell types in rows of the plot
row_anno <- as.data.frame(rownames(jia2))
colnames(row_anno) <- "celltype"
row_anno$type <- substr(row_anno$celltype, 0,2)
row_anno <- row_anno %>% column_to_rownames("celltype")
row_anno["M-u1: SERPINF1+ sublining", "type"] <- "F-"
cell_type_colours <- c("M-" = "forestgreen", "F-" = "brown", "E-"="gold", "T-" = "navy",
                       "NK" = "lightblue", "B-" = "grey")

# Generate distinct colours for each unique level of the annotation factors
syn_anno$Direction <- as.factor(syn_anno$Direction)
direction_colours <- setNames(c("red", "white")[1:length(levels(syn_anno$Direction))], levels(syn_anno$Direction))

# Combine all annotation colours into a list for pheatmap
annotation_colours_list <- list(
  Direction = direction_colours,
  age_group = age_colours,
  type = cell_type_colours,
  Neutrophils=c("white", "forestgreen"))

syn2 <- syn_anno[colnames(jia2),]

heatcolx <- c("#2166AC", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7","#F4A582", "#D6604D" ,"#B2182B","#B2182B")

# Plot JIA synovium specificity scores
pheatmap(jia2, main = "JIA specificity scores", color=heatcolx, 
               annotation_col = syn2, annotation_row = row_anno,
               annotation_colors = annotation_colours_list, 
               na_col = "white", fontsize = 7, border_color = "grey91", 
               scale="column",
               cluster_rows=F, cluster_cols = F)

# Plot RA synovium specificity scores
pheatmap(ed2, main = "RA specificity scores", color=heatcolx, 
         annotation_col = syn2, annotation_row = row_anno,
         annotation_colors = annotation_colours_list, 
         na_col = "grey91", fontsize = 7, border_color = "grey91", 
         scale="column",
         cluster_rows=F, cluster_cols = F)

#######################################################
#######################################################
#######################################################

# Create a table of max cell type for LOF and GOF genes

data <- t(ed2)
max_col_indices <- apply(data, 1, which.max)
max_col <- unique(rownames(data)[max_col_indices])
direction <- syndromes[,c("Direction", "Gene_D")]
max <- data.frame(Gene_D = colnames(ed2), Cell = colnames(data)[max_col_indices])
max <- merge(max, direction, by="Gene_D")
max3 <- table(max$Cell, max$Direction) %>% as.data.frame() 
max3 <- max3 %>% pivot_wider(names_from=c(Var2), values_from=Freq)

data <- t(jia2)
max_col_indices <- apply(data, 1, which.max)
max_col <- unique(rownames(data)[max_col_indices])
direction <- syndromes[,c("Direction", "Gene_D")]
colnames(direction) <- c("Direction_JIA", "Gene_D")
maxj <- data.frame(Gene_D = colnames(jia2), Cell = colnames(data)[max_col_indices])
maxj <- merge(maxj, direction, by="Gene_D")
max2 <- table(maxj$Cell, maxj$Direction) %>% as.data.frame()
max2 <- max2 %>% pivot_wider(names_from=c(Var2), values_from=Freq)

maxb <- full_join(max2, max3, by="Var1")

#######################################################

# Identify top expressing cell by gene
data <- t(jia2)
max_col_indices <- apply(data, 1, which.max)
max_col <- colnames(data)[max_col_indices]
top_cells <- data.frame(Gene_D=colnames(jia2), Cells = max_col)
write.csv(top_cells, file="CSV/2601_top_cells_by_gene_JIA.csv")

data <- t(ed2)
max_col_indices <- apply(data, 1, which.max)
max_col <- colnames(data)[max_col_indices]
top_cells <- data.frame(Gene_D=colnames(ed2), Cells = max_col)
write.csv(top_cells, file="CSV/2602_top_cells_by_gene_RA.csv")


#######################################################
#######################################################
#######################################################
# GOF gene expression in Neutrophils
syn_brief <- subset(syndromes, select=c(Gene_D, Onset_median, Direction))

# Age analysis
sc_data <- ed2
max_celltype <- sc_data %>% t() %>% as.data.frame()
max_celltype <- apply(max_celltype, 1, function(x) {
  colnames(max_celltype)[which.max(x)]})
max_celltype <- max_celltype %>% as.data.frame() %>% rownames_to_column()
colnames(max_celltype) <- c("Gene_D", "celltype")
age.corr <- merge(max_celltype, syn_brief, by="Gene_D")
direction.types <- age.corr %>% dplyr::count(celltype, Direction)

# GOF gene set enrichment analysis
direct <- syndromes[,c("Direction", "Gene_D"),drop=FALSE]
direct$Genes <- rownames(direct)

ed_rows <- t(sc_data) %>% as.data.frame() %>% rownames_to_column("Gene_D")
gof_analysis <- merge(ed_rows, direct, by="Gene_D")

wilcox.test(`M-15: Neutrophils` ~ Direction, data = gof_analysis)

ggplot(gof_analysis, aes(x = Direction, y = `M-15: Neutrophils`, fill = Direction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) + geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1), limits=c(0,0.15)) +
  labs(title = "Neutrophil expression of GOF vs LOF genes", y = "log(Expression + 1)", x = "") +
  scale_fill_manual(values = c("LOF" = "darkgrey", "GOF" = "red3")) + theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

#######################################################
#######################################################
#######################################################

# Plot Neutrophil vs IL1B macrophage expression
# See which correlates most with neutrophils
expr <- t(df_select.ed) %>% as.data.frame()       
target_col <- "M-15: Neutrophils"

# Get target vector and other columns
target <- as.numeric(expr[, target_col])
other_cols <- setdiff(colnames(expr), target_col)

# Function to compute correlation + p-value + n
cor_stats <- function(x, y, method = c("pearson","spearman")) {
  method <- match.arg(method)
  ok <- which(!is.na(x) & !is.na(y))
  n <- length(ok)
  if (n < 3) return(list(r = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  list(r = unname(ct$estimate), p = ct$p.value, n = n)
}

# Iterate over columns
res_list <- lapply(other_cols, function(col) {
  vec <- as.numeric(expr[, col])
  p <- cor_stats(vec, target, method = "pearson")
  s <- cor_stats(vec, target, method = "spearman")
  data.frame(cell_type = col,pearson_r = p$r,pearson_p = p$p,
             spearman_rho = s$r,spearman_p = s$p,
             n = p$n)
})
res <- do.call(rbind, res_list)

# Add absolute-correlation columns and sort (by absolute Pearson by default)
res$abs_pearson <- abs(res$pearson_r)
res$abs_spearman <- abs(res$spearman_rho)
res_sorted <- res[order(-res$abs_pearson, -res$abs_spearman), ]
head(res_sorted, 10)

# correlation
cor_val <- cor(x=expr$`M-15: Neutrophils`, y=expr$`M-7: IL1B+ FCN1+`, use = "pairwise.complete.obs", method = "pearson")
cor_val

# Create labels of high values
label_df <- expr %>% filter(`M-7: IL1B+ FCN1+` > 0.051|`M-15: Neutrophils` > 0.07) %>%
  { tmp <- .; tmp$gene <- rownames(tmp); tmp }

# scatterplot with linear regression line
ggplot(expr, aes(x = `M-15: Neutrophils`, y = `M-7: IL1B+ FCN1+`)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", linetype = "dashed") +
  geom_text(data = label_df, aes(label = gene),size = 5, hjust=0.5, vjust=2) +
  scale_y_continuous(trans = 'log1p') + scale_x_continuous(trans = 'log1p') +
  annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(cor_val, 2)),
           hjust = 1.1, vjust = 1.5, size = 6, color = "black") +
  theme_minimal(base_size = 12) + theme(panel.grid=element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))

#######################################################
#######################################################
#######################################################

# Age analysis- looking at cell types with maximum specific expression vs onset age in RA dataset
pbe <- aggregateToPseudoBulk(sce.ed, assay = "counts",
                             cluster_id = "global_detail",
                             sample_id = "patient",
                             verbose = FALSE)

# Compute cell type specificity of each gene
dfe <- cellTypeSpecificity(pbe)

genes_present <- intersect(gen, rownames(synov.low))

ed3 <- dfe[genes_present,] 
ed3 <- ed3 %>% as.data.frame()
ed3$totalCPM <- NULL

# Age analysis
max_celltype <- apply(ed3, 1, function(x) {
  colnames(ed3)[which.max(x)]})

# Build a dataframe with gene names + corresponding max cell type
etop_main <- data.frame(Gene = rownames(ed3), Type = max_celltype, stringsAsFactors = FALSE)
etop_main <- as.data.frame(t(etop_main)) %>%
  mutate(`PIK3R1_LOF` = PIK3R1, `PIK3R1_GOF` = PIK3R1, 
         `PIK3CD_LOF` = PIK3CD, `PIK3CD_GOF` = PIK3CD) %>% dplyr::select(-PIK3R1, -PIK3CD)
etop_main <- etop_main %>% t() %>% as.data.frame() %>% rownames_to_column("Gene_D")

syndromes <- read.csv(file="CSV/2601_Syndromes.csv")
syndromes <- syndromes %>% filter(exclude == "N")
syn_brief <- subset(syndromes, select=c(Gene_D, Onset_median))

# Group ages
etop_main <- merge(etop_main, syn_brief, by="Gene_D")
eall_age <- etop_main %>% mutate(
  age_group = cut(Onset_median, breaks = c(-Inf, 2.9, 6.9, 16.9, Inf), labels = c("0-2", "3-6", "7-16", ">17")))

# Count number of observations per Type × age_group
eall_age <- na.omit(eall_age)
eage_counts <- eall_age %>% dplyr::count(age_group, Type)

age_colours <- c("0-2" = "blue", "3-6" = "lightblue1", "7-16" = "lemonchiffon", ">17" = "cornsilk3")
eage_counts$age_group <- factor(eage_counts$age_group, levels = rev(c("0-2", "3-6", "7-16", ">17")))
eage_counts$Type <- factor(eage_counts$Type, levels=c("Granulocytes","Myeloid.cells",  "T.cells", "pDC", "B.cells","Stromal.cells", "Endothelial.cells","Dendritic.cells"))

ggplot(eage_counts, aes(x = Type, y = n, fill = age_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=age_colours, guide = guide_legend(reverse = TRUE)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(color="black", angle =90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color="black")) + 
  labs(fill = "Age group", x="", y="")


# Reformat for stats
eage_stats <- eage_counts[eage_counts$Type %in% c("Granulocytes", "T.cells"),]
eage_stats <- eage_stats %>% mutate(early = ifelse(age_group == "0-2", "0-2", ">2"))
eage_stats <- eage_stats %>% mutate(Cells = ifelse(Type == "Granulocytes", "Granulocytes", "Other"))

eage_m <- eage_stats %>% group_by(early, Cells) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = early, values_from = n, values_fill = 0) %>%
  column_to_rownames("Cells") %>% as.matrix()
eage_m

# Fisher exact test
fisher_res <- fisher.test(eage_m)
fisher_res

#######################################################
#######################################################
#######################################################

