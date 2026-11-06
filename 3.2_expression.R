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
library(dreamlet)

heatcol <- rev(brewer.pal(8, "RdBu"))
heatcolx <- c("#2166AC", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7","#F4A582", "#D6604D" ,"#B2182B","#B2182B")

options(bitmapType='cairo')

#######################################################
#######################################################
#######################################################

# Plot PBMC expression of monogenic genes
syndromes <- read.csv(file="CSV/monogenic_sheet_1.csv", header=T)

excluded_genes <- syndromes %>% filter(exclude == "Y") %>% pull(Gene)
gen <- syndromes %>% filter(exclude == "N") %>% pull(Gene)

#######################################################
#######################################################
#######################################################

# Load Human Protein Atlas RNA Tissue consensus fil, available here: https://www.proteinatlas.org/humanproteome/tissue/data#consensus_tissues_rna
hpa <- read.delim(file="Datasets/rna_tissue_consensus.tsv")
hpa_gen <- hpa[hpa$Gene.name %in% gen,]

# Find missing genes
setdiff(gen, unique(hpa_gen$Gene.name)) # Not found: "RMRP", "TMEM173", "DDX58"

# Add aliases to list- "DDX58"="RIGI", "TMEM173"="STING1"
gen <- c(gen, "RIGI", "STING1")
hpa_gen <- hpa[hpa$Gene.name %in% gen,]

# Reformat
hpa_gen1 <- reshape2::dcast(hpa_gen, formula = Gene.name ~ Tissue, value.var= "nTPM")
rownames(hpa_gen1) <- hpa_gen1$Gene.name

# Correct gene names for consistency
rownames(hpa_gen1)[rownames(hpa_gen1)=="RIGI"] <- "DDX58"
rownames(hpa_gen1)[rownames(hpa_gen1)=="STING1"] <- "TMEM173"

# Ensure excluded genes are not there
hpa_gen1 <- hpa_gen1[!(rownames(hpa_gen1) %in% excluded_genes),]
hpa_gen1 <- hpa_gen1 %>% dplyr::select(-Gene.name)
hpa_gen1 <- hpa_gen1 %>% as.matrix() %>% t()

heatcol <- rev(brewer.pal(8, "RdBu"))
pheatmap(hpa_gen1, scale="column", col= heatcol, main="Tissue Expression", fontsize=8)

# Visualise top expressing tissues for clearer visualisation
# Calculate row means (average expression across genes for each tissue) and take top 10
data <- hpa_gen1 %>% as.data.frame()
data$avg_expression <- rowMeans(data, na.rm = TRUE)
top_tissues <- head(data[order(-data$avg_expression), ], 15)
rownames(top_tissues)
top_tissues$avg_expression <- NULL

tissue_order <- c("thymus","lymph node", "tonsil", "appendix", "bone marrow", "spleen", 
                  "lung", "adipose tissue", "urinary bladder","choroid plexus", "colon", "duodenum", "small intestine", "liver")
top_order <- top_tissues[tissue_order,]
rownames(top_order) <- str_to_sentence(rownames(top_order))

pheatmap(as.matrix(top_order), scale="column", col= heatcol, 
         main="Tissue Expression", cluster_rows=FALSE,
         fontsize=7, border_color = "grey", clustering_distance_cols = "correlation") 

# Create PCA of expression data
pca_result <- prcomp(hpa_gen1, scale. = TRUE)
pca_scores <- data.frame(pca_result$x[, 1:2])
pca_scores$Labels <- str_to_sentence(row.names(hpa_gen1))

# Create PCA plot
ggplot(data = pca_scores, aes(x = PC1, y = PC2, label = Labels)) +
  geom_point(color = "navy", size = 2) +
  geom_text(aes(label = Labels), hjust = 0.5, vjust = -1, size = 3.5) +
  theme_minimal() +
  coord_cartesian(clip = "off") + # Prevent clipping of labels
  theme(plot.margin = unit(c(1, 1, 1, 3), "cm"), panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))


#######################################################
#######################################################
#######################################################

# Visualise aggregate gene counts in JIA dataset, coloured by max tissue expression

# Load JIA tissue dataset
load(file="/rds/projects/c/croftap-mapjagbatch2/Chrissy_monogenic_objects/JIA-tissue-annotated.rds")

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

load(file="/rds/projects/c/croftap-mapjagdata/MAPJAGv2/2306/Global/Monogenic/Run-relabel/synov_annotated.RData")

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
gen_hpa <- c(gen, "RIGI", "STING1")
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

# Add in arthropathy genes
arthrop <- c("MMP2", "PRG4", "IDUA", "CCN6")

gen2 <- unique(c(gen, arthrop))

# Create SCE object for analysis and pseudobulk
sce <- as.SingleCellExperiment(tissue.small)
pb <- aggregateToPseudoBulk(sce, assay = "counts", cluster_id = "RA_clusters", sample_id = "patient")

# Compute cell type specificity of each gene
df <- cellTypeSpecificity(pb)
genes_present <- intersect(gen2, rownames(tissue))
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

genes_present <- intersect(gen2, rownames(synov))
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

# Create df without arthropathy genes
df_select_gen_ed <- df_select.ed[, !colnames(df_select.ed) %in% arthrop]
data <- df_select_gen_ed %>% t()
max_col_indices <- apply(data, 1, which.max)
max_col <- unique(colnames(data)[max_col_indices])
filtered_ed3 <- data[, max_col]
top_ed2 <- colnames(filtered_ed3)

# Create df without arthropathy genes
df_select_gen <- df_select[, !colnames(df_select) %in% arthrop]
data <- df_select_gen %>% t()
max_col_indices <- apply(data, 1, which.max)
max_col <- unique(colnames(data)[max_col_indices])
filtered_jia3 <- data[, max_col]
top_jia2 <- colnames(filtered_jia3)

save.image(file="2609_expr_clustering_objects.RData")

#######################################################

# Combine matrices to support finding an order for visualisation of expression heatmaps
# Top cell types across both datasets
top_both <- unique(c(top_jia,top_ed))

mat1 <- as.matrix(t(df_select))
mat2 <- as.matrix(t(df_select.ed))

# Union of row and column names
all_genes <- unique(c(rownames(mat1), c(rownames(mat2))))

# Function to expand matrix and fill missing entries with 0
expand_with_zeros <- function(mat, rows, cols) {
  out <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
  r <- intersect(rownames(mat), rows); c <- intersect(colnames(mat), cols)
  out[r, c] <- mat[r, c, drop = FALSE]
  out
}

mat1_full <- expand_with_zeros(mat1, all_genes, top_both)
mat2_full <- expand_with_zeros(mat2, all_genes, top_both)

#######################################################

cell_dend2 <- c("M-3: MERTK+ HBEGF+", "M-4: SPP1+", "M-2: MERTK+ S100A8+", "M-0: MERTK+ SELENOP+ LYVE1+", "M-1: MERTK+ SELENOP+ LYVE1-", 
                "M-5: C1QA+", "M-10: DC2", "M-9: DC3", "M-6: STAT1+ CXCL10+", "M-11: DC4", "M-7: IL1B+ FCN1+", "M-15: Neutrophils", 
                "M-12: DC1", "M-14: LAMP3+", "M-8: PLCG2+", "M-13: pDC",
                "B-3: CD24++CD27+IgM+ unswitched memory", "B-0: CD24+CD27+CD11b+ switched memory", "B-1: IgM+IgD+TCL1A+ naive", "B-2: IgG1+IgG3+ plasma",
                "NK-7: CD56bright CD16- GNLY+CD69+", "NK-5: CD56bright CD16- GZMA+CD69+", 
                "NK-1: CD56dim CD16+ IFNG+CD160+", "NK-0: CD56dim CD16+ IFNG-", "T-12: CD4+ GNLY+", "T-13: CD8+ GZMK/B+ memory", 
                "T-5: CD4+ GZMK+ memory", "T-14: CD8+ GZMK+ memory", "T-6: CD4+ memory", "T-1: CD4+ CD161+ memory", "T-2: CD4+ IL7R+CCR5+ memory", 
                "T-21: Innate-like", "T-0: CD4+ IL7R+ memory", "T-4: CD4+ naive", "T-9: CD4+ CD25-low Treg", "T-10: CD4+ OX40+NR3C1+", 
                "T-3: CD4+ IFNG- Tph/Tfh", "T-7: CD4+ IFNG+ Tph/Tfh", "T-8: CD4+ CD25-high Treg", "T-18: Proliferating", 
                "E-6: Lymphatic", "E-3: Arterial", "E-2: ICAM1+ venous", "E-0: SELP-hi venous", "M-u1: SERPINF1+ sublining", 
                "F-1: SFRP1+ sublining", "F-6: Proliferating sublining", "F-5: CD82+ sublining", "F-4: POSTN+ sublining", 
                "F-8: DKK3+ sublining", "F-2: CD34+ sublining", "F-3: CLIC5+ lining", "F-0: VCAM1+ lining")

gene_dend <- c("UNC93B1", "ADA2", "SLC29A3", "LACC1", "C1QC", "C1QA", "C1QB", 
               "C2", "MAN2B2", "ASAH1", "ACP5", "IL10RB", "DNASE2", "TYMP", "SPI1", "CARD14", "NCF1",
               "XIAP", "IFIH1", "STAT1", "DOCK11", "NOD2", "ALPK1", "GLA", "BTK",
               "LYN", "SYK", "COPA", "NFIL3", "MEFV", "NLRP3", "TBK1", "NCF2", "NLRC4", 
               "WAS", "ELF4", "TREX1", "NLRP12", "DDX58", "ARPC1B", "SAMHD1", "C2orf69", "DCLRE1C", 
               "RNASEH2B", "EPCAM", "DNASE1L3","RNF31", "RELB", "NFKB2", "NFKBIA", "NFKB1", 
               "OTULIN","RIPK1", "RMRP","AIRE", "IRF2BP2",
               "RAG1", "MVK", "AICDA", "RAG2", "SHOC2", "TNFRSF13B", "PLCG2", "TRNT1", "ATM", 
               "PRF1", "PIK3R1_LOF", "PIK3R1_GOF", "RELA", "IKBKG",  "IL10RA", "NLRP1", 
               "TNFAIP3", "PSTPIP1", "UNC13D", "SOCS1", "PIK3CD_LOF", "PIK3CD_GOF", "IKZF1", "CARMIL2", "DEF6", 
               "IL2RG", "RAB27A", "LRBA", "CD27", "IFITM5", "STAT3", "CTLA4", "ICOS", "CD40LG", "PDCD1", "FOXP3", 
               "RBCK1", "PSMB8", "RFXANK","PMVK", "TMEM173", "GATA2", "HYAL1", "PRKDC", "TNFRSF1A", 
               "G6PC3", "SHARPIN",  "C4A", "C4B", "C1S","IDUA", "MMP2", "CCN6", "PRG4")

# Set order of JIA and RA datasets
jia2 <- mat1_full %>% t() %>% as.data.frame()
ed2 <- mat2_full %>% t() %>% as.data.frame()
jia2 <- jia2[cell_dend2, gene_dend]
ed2 <- ed2[cell_dend2, gene_dend]
jia2[is.na(jia2)] <- 0
ed2[is.na(ed2)] <- 0
jia2 <- jia2[rowSums(jia2) > 0, ]
ed2 <- ed2[rowSums(ed2) > 0, ]

#######################################################
#######################################################

# Plot with direction, neutrophils and onset age annotations
# Ensure LOF and GOF versions of gene defects are present and remove original columns
# Remove excluded genes
syn_brief <- syndromes %>% filter(exclude == "N")

# Cut ages into categories
syn_brief <- syn_brief %>% mutate(
  age_group = cut(Onset_median, breaks = c(-Inf, 2.9, 6.9, 16.9,  Inf), labels = c("0-2", "3-6", "7-16", ">17")))

syn_brief <- subset(syn_brief, select=c(Gene_D, age_group, Direction, Neutrophils))
syn_brief[, ][is.na(syn_brief[, ])] <- 0 # Convert NA values to 0
syn_anno <- syn_brief %>% column_to_rownames("Gene_D")

# Prepare annotation data frame for pheatmap, convert categorical columns to factors for proper annotation coloring
age_colours <- c("0-2" = "blue", "3-6" = "lightblue1", "7-16" = "lemonchiffon", ">17"= "cornsilk3")

# Add cell type colours for cell types in rows of the plot
row_anno <- as.data.frame(rownames(jia2))
colnames(row_anno) <- "celltype"
row_anno$type <- substr(row_anno$celltype, 0,2)
row_anno <- row_anno %>% column_to_rownames("celltype")
row_anno["M-u1: SERPINF1+ sublining", "type"] <- "F-"
row_anno["M-13: pDC", "type"] <- "P-"
cell_type_colours <- c("M-" = "forestgreen", "F-" = "brown", "E-"="gold", "T-" = "navy",
                       "NK" = "lightblue", "B-" = "hotpink", "P-"="lavender")

# Generate distinct colours for each unique level of the annotation factors
syn_anno$Direction <- as.factor(syn_anno$Direction)
direction_colours <- setNames(c("red", "white")[1:length(levels(syn_anno$Direction))], levels(syn_anno$Direction))

# Combine all annotation colours into a list for pheatmap
annotation_colours_list <- list(
  Direction = direction_colours,
  age_group = age_colours,
  type = cell_type_colours)

heatcolx <- c("#2166AC", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7","#F4A582", "#D6604D" ,"#B2182B","#B2182B")

#######################################################
#######################################################

# Plot JIA synovium specificity scores
pheatmap(jia2, main = "JIA specificity scores", color=heatcolx, 
         annotation_col = syn_anno, annotation_row = row_anno,
         annotation_colors = annotation_colours_list, 
         na_col = "white", fontsize = 7, border_color = "grey91", 
         scale="column",
         cluster_rows=F, cluster_cols = F)

pheatmap(ed2, main = "RA specificity scores", color=heatcolx, 
         annotation_col = syn_anno, annotation_row = row_anno,
         annotation_colors = annotation_colours_list, 
         na_col = "white", fontsize = 7, border_color = "grey91", 
         scale="column",
         cluster_rows=F, cluster_cols = F)

#######################################################
#######################################################

# Check expression in PBMC dataset from Human Protein Atlas, available here: https://www.proteinatlas.org/download/tsv/rna_immune_cell_monaco.tsv.zip
immune_full <- read.delim(file="Datasets/rna_immune_cell_monaco.tsv")

# Find genes in dataset
immune <- immune_full[immune_full$Gene.name %in% gen,]

# Find missing genes
setdiff(gen, unique(immune$Gene.name))

# Add gene aliases
gen <- c(gen, "STING1", "RIGI")

# Reformat
immune <- immune_full[immune_full$Gene.name %in% gen,]
setdiff(gen, unique(immune$Gene.name))

immune1 <- reshape2::dcast(immune, formula = Gene.name ~ Immune.cell, value.var= "TPM")
rownames(immune1) <- immune1$Gene.name
immune1 <- immune1 %>% as.data.frame() %>% dplyr::select(-c(Gene.name, "total PBMC"))

# Amend gene names to align with other datasets
rownames(immune1)[rownames(immune1) == "STING1"] <- "TMEM173"
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

pheatmap(as.matrix(immune2), scale="column", col= heatcol, 
         main="PBMC Expression", fontsize=7, border_color = "grey",  cutree_cols = 6,
         na_col = "white",
         clustering_distance_cols = "correlation")

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

arth <- read.csv(file="CSV/arthritis_csv.csv", header = T, row.names=1)

arth[arth == ""] <- 0
arth[arth == "NR"] <- 0

top_cells <- read.csv("CSV/2602_top_cells_by_gene_RA.csv")

top_cells <- top_cells %>%
  mutate(macrophage = str_detect(Cells, "MERTK|SPP1"))

arth$Gene_D <- rownames(arth)

ana_dat <- top_cells %>%
  inner_join(arth %>% dplyr::select(Gene_D, ANA_pos),
             by=c("Gene_D"="Gene_D"))

tab <- table(ANA_positive = ana_dat$ANA_pos == 1,
             Macrophage = ana_dat$macrophage)

tab
fisher.test(tab)

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
eage_counts$Type <- gsub("\\.", " ", eage_counts$Type)
eage_counts$Type <- factor(eage_counts$Type, levels=rev(c("Granulocytes","Myeloid cells",  "T cells", "pDC", "B cells","Stromal cells", "Endothelial cells","Dendritic cells")))

ggplot(eage_counts, aes(y = Type, x = n, fill = age_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=age_colours, guide = guide_legend(reverse = TRUE)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(color="black", angle =0, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color="black")) + 
  labs(fill = "Age group", x="", y="")


# Reformat for stats
eage_stats <- eage_counts[eage_counts$Type %in% c("Granulocytes", "T cells"),]
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

# Autoimmunity disorders by age
medians2 <- syndromes %>% filter(exclude == "N") %>% dplyr::select(c(Gene_D, Autoinflammatory, Onset_median))

medians2 <- medians2 %>% mutate(
  age_group = cut(Onset_median, breaks = c(-Inf, 2.9, 6.9, 16.9, Inf), labels = c("0-2", "3-6", "7-16", ">17")))

medians2$Autoinflammatory <- ifelse(medians2$Autoinflammatory == 1, 1, 0)
medians2$Autoinflammatory[is.na(medians2$Autoinflammatory)] <- 0

med_tab <- table(medians2$age_group, medians2$Autoinflammatory) %>% as.data.frame()

plot_age <- medians2 %>% filter(!is.na(age_group)) %>%
  group_by(age_group, Autoinflammatory) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(age_group = factor(age_group, levels = c("0-2", "3-6", "7-16", ">17")),
    Autoinflammatory = factor(Autoinflammatory, levels = c(0, 1), labels = c("No", "Yes")))

ggplot(plot_age, aes(fill = age_group,y = n,x = Autoinflammatory)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) + 
  scale_fill_manual(values=age_colours, guide = guide_legend(reverse = TRUE)) +
  theme_bw(base_size = 14) +
  labs( x = "Autoinflammatory features reported in cases",
    y = "Number of gene defects",
    fill = "Median age of onset") + theme_classic() 

#######################################################
#######################################################
#######################################################




