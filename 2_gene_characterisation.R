
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


#######################################################
#######################################################
#######################################################

# Plot tissue transcription of genes 
# Import Gene Summary data from Data Object sheet 1 and reformat 
syndromes <- read.csv(file="monogenic_sheet_1.csv", header=F, blank.lines.skip = TRUE)
syndromes <- syndromes[3:nrow(syndromes),]
colnames(syndromes) <- syndromes[1, ]
syndromes <- syndromes[-1, ]

excluded_genes <- syndromes %>% filter(exclude == "Y") %>% pull(Gene) # Remove genes with low penetrance
gen <- syndromes %>% filter(exclude == "N") %>% pull(Gene)

# Load Human Protein Atlas RNA Tissue consensus fil, available here: https://www.proteinatlas.org/humanproteome/tissue/data#consensus_tissues_rna
hpa <- read.delim(file="Datasets/rna_tissue_consensus.tsv")
hpa_gen <- hpa[hpa$Gene.name %in% gen,]

# Find missing genes
setdiff(gen, unique(hpa_gen$Gene.name)) # Not found: "RMRP", "TMEM173", "DDX58"

# Add aliases to list- "DDX58"="RIGI", "GBA" ="GBA1", "TMEM173"="STING1"
gen <- c(gen, "RIGI", "GBA1", "STING1")
hpa_gen <- hpa[hpa$Gene.name %in% gen,]

# Reformat
hpa_gen1 <- reshape2::dcast(hpa_gen, formula = Gene.name ~ Tissue, value.var= "nTPM")
rownames(hpa_gen1) <- hpa_gen1$Gene.name

# Correct gene names for consistency
rownames(hpa_gen1)[rownames(hpa_gen1)=="RIGI"] <- "DDX58"
rownames(hpa_gen1)[rownames(hpa_gen1)=="STING1"] <- "TMEM173"
rownames(hpa_gen1)[rownames(hpa_gen1)=="GBA1"] <- "GBA"

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

# To determine representation in immune categories of IUIS Inborn Errors of Immunity, first find genes represented in IUIS tables
# Available from monogenic sheet 3 in Data_Objects
iuis <- read.csv(file="monogenic_sheet_3.csv", header=T)

# Add gene alias
iuis_gen <- iuis[grepl(paste(gen, collapse ="|"), iuis$Genetic.defect),]
write.csv(iuis_gen, file="CSV/iuis_arthritis-genes-sub.csv") 

# Find IUIS overlap of categories
iuis <- distinct(iuis)
iuis$Major.category <- str_to_title(iuis$Major.category) %>% trimws(, which = "both") # trim white space
iuis$Subcategory <- str_to_title(iuis$Subcategory) %>% trimws(, which = "both")

# Count by categories for all genes in IUIS table
major_sub <- iuis %>% count(Major.category, Subcategory, name = "n")

# Format for combining
major_sub$m <- gsub("Table", "", major_sub$Major.category) 
major_sub$o <- substr(major_sub$m,1,2) # Category number
major_sub$Subcategory <- gsub("Subtable ", "", major_sub$Subcategory)
major_sub$cat <- paste0(major_sub$o,"_",major_sub$Subcategory) # Add main category number to subcategory label
major_reduced <- major_sub[,c("n", "cat")]

# Correct formatting issue re duplicate category
major_reduced[46,1] <- 41
major_reduced <- major_reduced[-45,]

# Repeat for arthritis-associated genes
# Note that some genes do not have the correct directionality (GOF / LOF), or have duplicate entries, need to manually select these out
# So manually check the included genes for inclusion errors and the correct versions of syndromes

# Available from monogenic sheet 4 in Data_Objects
arthritis.iuis <- read.csv(file="monogenic_sheet_4.csv")

# Remove excluded genes
arthritis.iuis$Genetic.defect <- arthritis.iuis$Genetic.defect %>% trimws(, which = "both")
arthritis.iuis <- arthritis.iuis %>% filter( !(Genetic.defect %in% excluded_genes))

# Reformatting names
arthritis.iuis$Major.category <- str_to_title(arthritis.iuis$Major.category) %>% trimws(, which = "both")
arthritis.iuis$Subcategory <- str_to_title(arthritis.iuis$Subcategory) %>% trimws(, which = "both")

# Count categories represented amongst arthritis genes, having manually selected GOF/LOF disorders
arthritis_sub <- arthritis.iuis %>% count(Major.category, Subcategory, name = "genes")

# Format subcategory names to include main category numbering
arthritis_sub$m <- gsub("Table", "", arthritis_sub$Major.category) 
arthritis_sub$o <- substr(arthritis_sub$m,1,2)
arthritis_sub$Subcategory <- gsub("Subtable ", "", arthritis_sub$Subcategory)
arthritis_sub$cat <- paste0(arthritis_sub$o,"_",arthritis_sub$Subcategory)
arthritis_reduced <- arthritis_sub[,c("genes", "cat")]

# Combine genes counts and total counts into one table
enriched.table <- left_join(major_reduced, arthritis_reduced, by="cat")
enriched.table <- enriched.table[,c(2,1,3)]
enriched.table$genes[is.na(enriched.table$genes)] <- 0

# Determine enrichment of arthritic genes compared to total
enriched_stats <- enriched.table %>% rowwise() %>% mutate(total_genes = sum(enriched.table$n),
    total_subset = sum(enriched.table$genes), expected_prop = n / total_genes, observed_prop = genes / total_subset,
    enrichment = observed_prop / expected_prop,
    pval = fisher.test(matrix(c(genes, total_subset - genes, n - genes, total_genes - n - (total_subset - genes)), nrow = 2))$p.value) %>%
  ungroup() %>% mutate(padj = p.adjust(pval, method = "BH"))

# Scale enrichment (z-score style)
enriched_stats <- enriched_stats %>%
  mutate(enrichment_scaled = scale(enrichment)[, 1],  sig_label = ifelse(padj < 0.05, "*", ""))

# Plot as one row of tiles
ggplot(enriched_stats, aes(x = cat, y = 1, fill = enrichment_scaled)) +
  geom_tile(color = "grey80") + geom_text(aes(label = sig_label), vjust = 0.5, size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Scaled enrichment") +
  labs(x = NULL, y = NULL, title = "Category enrichment") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

enriched_stats$cat <- factor(enriched_stats$cat, levels=rev(enriched_stats$cat))

#######################################################
#######################################################
#######################################################

# Repeat for major IUIS categories
major_sub <- iuis %>% count(Major.category, name = "n")
arthritis_sub <- arthritis.iuis %>% count(Major.category, name = "genes")

# Combine genes and total table
enriched2 <- left_join(major_sub, arthritis_sub, by="Major.category")
enriched2$Major.category <- gsub("Table", "", enriched2$Major.category) 
enriched2$genes[is.na(enriched2$genes)] <- 0

# Add stats
enriched_stats2 <- enriched2 %>% rowwise() %>% mutate(
    total_genes = sum(enriched2$n),
    total_subset = sum(enriched2$genes),
    expected_prop = n / total_genes,
    observed_prop = genes / total_subset,
    enrichment = observed_prop / expected_prop,
    pval = fisher.test(matrix(c(genes, total_subset - genes, n - genes, total_genes - n - (total_subset - genes)), nrow = 2))$p.value) %>%
  ungroup() %>% mutate(padj = p.adjust(pval, method = "BH"))

# Scale enrichment (z-score style)
enriched_stats2 <- enriched_stats2 %>%
  mutate(enrichment_scaled = scale(enrichment)[, 1],  sig_label = ifelse(padj < 0.05, "*", ""))

# Plot as one row of tiles
ggplot(enriched_stats2, aes(x = Major.category, y = 1, fill = enrichment_scaled)) +
  geom_tile(color = "grey80") + geom_text(aes(label = sig_label), vjust = 0.5, size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Scaled enrichment") +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

#######################################################
#######################################################
#######################################################

# Plot immune categories and phenotypic features
# Cluster by IUIS category
syndromes_tb <- syndromes %>% filter(exclude=="N")
syndromes_tb <- syndromes_tb %>% dplyr::select(c(Category, Gene_D))

# Create matrix of column values
syndromes_tb <- syndromes_tb %>% mutate(Category = replace_na(str_trim(na_if(Category "")), "None"))

syndromes_tb <- syndromes_tb %>% pivot_wider(names_from = Category, names_prefix = "", 
                                             values_from = Category, values_fn = ~1, values_fill = 0)

syndromes_tb <- syndromes_tb %>% column_to_rownames("Gene_D")

pheatmap(t(syndromes_tb), color = c("white", "orange"), na_col = "white", 
         fontsize=7, border_color = "lightgrey", legend = F)

# Add annotation
syn_anno <- subset(syndromes, select=rev(c("Gene_D", "Direction","Immunodeficiency", "Autoinflammatory",
                                    "SLE", "IBD","HLH...MAS", "Vasculitis_Panniculitis", "Syndrome")))

# Replace blank with none
syn_anno <- syn_anno %>% mutate(across(Syndrome, ~ replace_na(str_trim(na_if(.x, "")), "None")))
rownames(syn_anno) <-syn_anno$Gene_D
syn_anno$Gene_D <- NULL
syn_anno[, ][is.na(syn_anno[, ])] <- 0 # Convert NA values to 0

# Prepare annotation data frame for pheatmap
syn_anno <- syn_anno[rownames(categories), , drop=FALSE] %>% as.data.frame()
syn_anno$Direction <- as.factor(syn_anno$Direction)
syn_anno$Syndrome <- as.factor(syn_anno$Syndrome)

syn_anno$Vasculitis_Panniculitis <- as.numeric(syn_anno$Vasculitis_Panniculitis)
syn_anno$HLH...MAS <- as.numeric(syn_anno$HLH...MAS)
syn_anno$Immunodeficiency <- as.numeric(syn_anno$Immunodeficiency)
syn_anno$Autoinflammatory <- as.numeric(syn_anno$Autoinflammatory)
syn_anno$IBD <- as.numeric(syn_anno$IBD)
syn_anno$SLE <- as.numeric(syn_anno$SLE)

# Generate distinct colours for each unique level of the annotation factors
direction_colours <- setNames(c("red", "white")[1:length(levels(syn_anno$Direction))], levels(syn_anno$Direction))
syndrome_colours <- c("Hyper IgM" = "purple3", "Wiskott-Aldrich /-like"  = "darkolivegreen2","LSD" = "grey45", "None"="white",
                      "Hyper IgD"  = "deeppink3",  "IPEX /-like"  = "cyan", "Hypomorphic SCID"  = "deepskyblue3",  
                      "AGS /-like"  = "cornsilk3", "FCAS" = "lightcoral",  "EDA-ID /-like" = "forestgreen", "CGD"= "goldenrod3")

# Combine all annotation colours into a list for pheatmap
annotation_colours_list <- list(
  Direction = direction_colours,
  Syndrome = syndrome_colours,
  Immunodeficiency = c("white", "wheat"), 
  Autoinflammatory = c("white", "plum2"), 
  SLE = c("white", "lightsalmon"), 
  HLH...MAS = c("white", "deeppink4"), 
  IBD = c("white", "tan"),
  Vasculitis_Panniculitis = c("white", "gold"))

categories <- syndromes_tb[,c(9:5,3,1,2,4)]
ph <- pheatmap(t(categories), color = c("white", "orange"), na_col = "white", 
         fontsize=7, border_color = "lightgrey", legend = F, cluster_rows = F,
         annotation_col = syn_anno, annotation_colors = annotation_colours_list)

gene_order = rownames(categories)[ ph$tree_col$order ]
gene_order

gene_dend <- gene_order[c(1:46,56:106,47:55)]
categories <- categories[gene_dend,]

pheatmap(t(categories), color = c("white", "orange"), na_col = "white", 
               fontsize=7, border_color = "lightgrey", legend = F, cluster_rows = F,
         cluster_cols=F, annotation_col = syn_anno, annotation_colors = annotation_colours_list)

#######################################################
#######################################################
#######################################################

# Distribution of onset age by immune category, generate medians from 1_gene_summary
age_long -> median_test
colnames(median_test) <- c("Gene", "Age")

# Cut into category groups
median_test <- median_test %>% mutate(age_group = cut(Age, breaks = c(-Inf, 2.9, 6.9, 16.9, Inf), 
                                                labels = c("0-2", "3-6", "7-16", ">17")))

syndromes_tb <- subset(syndromes_inc, select=c(Gene_D, Category))
median_test <- merge(median_test, syndromes_tb, by.x="Gene", by.y ="Gene_D", all=T)
median_test <- na.omit(median_test)
median_test <- median_test %>% mutate(across(Category, ~ replace_na(str_trim(na_if(.x, "")), "None")))

tbl <- table(median_test$age_group, median_test$Category)
chisq.test(tbl)

#######################################################
#######################################################
#######################################################

# Investigate IUIS immune phenotype by onset age
# Create age categories for median onset age of genes
syndromes_1 <- syndromes_inc %>% mutate(age_group = cut(as.numeric(Onset_median), 
                                                    breaks = c(-Inf, 2.9, 6.9, 16.9, Inf), 
                                                    labels = c("0-2", "3-6", "7-16", ">17")))
syndromes_tb <- subset(syndromes_1, select=c(age_group, Category))
syndromes_tb <- syndromes_tb %>% mutate(across(Category, ~ replace_na(str_trim(na_if(.x, "")), "None")))
syndromes_tb$age_group <- factor(syndromes_tb$age_group, levels = rev(c("0-2", "3-6", "7-16", ">17")))

# Count immune phenotype by age category
cat_age <- syndromes_tb %>% count(age_group, Category)
cat_age <- cat_age[!is.na(cat_age$age_group), ]
cat_age$Category <- factor(cat_age$Category, levels=rev(c("Autoinflammatory disorders",                                
                                                          "Combined immunodeficiency",                                 
                                                          "Complement",                                               
                                                          "Congenital defects of phagocytes",                          
                                                          "Defects in intrinsic and innate immunity",                  
                                                          "Immune dysregulation",                                      
                                                          "Immunodeficiencies affecting cellular and humoral immunity",
                                                          "Predominantly antibody deficiencies", "None")))

# Set colours for plotting
category_colours <- c("Autoinflammatory disorders"="deeppink4", "Predominantly antibody deficiencies"= "thistle", 
                      "Immune dysregulation" ="lightcyan2", "Combined immunodeficiency"="darkolivegreen3", "Complement"="coral",
                      "Immunodeficiencies affecting cellular and humoral immunity"="deepskyblue1", "None"="lightgrey",
                      "Congenital defects of phagocytes"="gold3", "Defects in intrinsic and innate immunity"="yellow2")

# Stacked bar plot (counts)
ggplot(cat_age, aes(y = age_group, x = n, fill = Category)) +
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values=category_colours, guide = guide_legend(reverse = TRUE)) +
  theme_bw(base_size = 14) + scale_x_continuous(labels = scales::percent)+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
  labs(fill = "Age group", x="", y="")

#######################################################
#######################################################
#######################################################

# Autoimmunity disorders by age
medians2 <- syndromes_1 %>% filter(exclude == "N") %>% dplyr::select(c(Gene_D, Autoinflammatory, Category, Onset_median, age_group))
medians2$Autoinflammatory <- ifelse(medians2$Autoinflammatory == 1, 1, 0)
medians2$Autoinflammatory[is.na(medians2$Autoinflammatory)] <- 0

# Create groups by age and auotinflammatory status
medians2 <- medians2 %>% mutate(age_group = ifelse(as.numeric(Onset_median) <17, "<17", "≥17"))
med_counts <- medians2 %>% count(age_group, Autoinflammatory) %>% as.data.frame()
med_counts <- med_counts[!is.na(med_counts$age_group),]
med_counts

# Fisher test
med_tab <- table(medians2$age_group, medians2$Autoinflammatory)
ft <- fisher.test(med_tab)
ft

med_counts$Autoinflammatory <- factor(med_counts$Autoinflammatory, levels=rev(c("0", "1")))

# Plot with p-value
ggplot(med_counts, aes(x = age_group, y = n, fill = Autoinflammatory)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw(base_size = 10) + scale_fill_manual(values = c("1"= "deeppink4", "0"= "grey70")) +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))


#######################################################
#######################################################
#######################################################

# Determine representation in metabolic classification (ICIMD)

# Metabolic gene summary  (Ferreira 2021), 
# Read in full classification table, data available in sheet 5 of the Data Object
metab <- read.csv(file="monogenic_sheet_5.csv")

# Find genes in metabolic database
metab_gen <- metab[metab$Gene.name %in% gen,]
metab_gen <- metab_gen[!duplicated(metab_gen$Gene.name), ]
# write.csv(metab_gen, file="CSV/metabolic_arthritis-genes.csv")

metab <- distinct(metab)
metab$X <- metab$X %>% trimws(, which = "both")
metab_gen$X <- metab_gen$X %>% trimws(, which = "both")

# Count by both together for all genes in metab table
major_sub <- metab %>% count(X, name = "n")

# Count represented genes
arthritis_sub <- metab_gen %>% count(X)

# Combine genes and total table
enriched.table <- left_join(major_sub, arthritis_sub, by="X")
enriched.table$n.y[is.na(enriched.table$n.y)] <- 0
colnames(enriched.table) <- c("Pathway", "Total", "Cases")

# compute universe / subset totals once
total_genes  <- sum(enriched.table$Total)
total_subset <- sum(enriched.table$Cases)

enriched_stats <- enriched.table %>% rowwise() %>%
  mutate(expected_prop = Total / total_genes,
         observed_prop = Cases / total_subset,
         enrichment = observed_prop / expected_prop,
         pval = {
           a <- Cases
           b <- Total - Cases
           c <- total_subset - Cases
           d <- total_genes - Total - c
           if(any(c(a,b,c,d) < 0)) NA_real_ else
             fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")$p.value}) %>% ungroup() %>% 
  mutate(padj = p.adjust(pval, method = "BH"))

# Scale enrichment (z-score style)
enriched_stats <- enriched_stats %>%
  mutate(enrichment_scaled = scale(enrichment)[, 1],  sig_label = ifelse(padj < 0.05, "*", ""))

absent <- enriched_stats %>% filter(Cases == 0) %>% pull(Pathway)
enriched <- enriched_stats %>% filter(Cases > 0)

# Plot as one row of tiles
ggplot(enriched, aes(x = Pathway, y = 1, fill = enrichment_scaled)) +
  geom_tile(color = "grey80") + geom_text(aes(label = sig_label), vjust = 0.5, size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Scaled enrichment") +
  labs(x = NULL, y = NULL, title = "Category enrichment") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

# Repeat for subgroup categories
# Read in full classification table, data available in sheet 5 of the Data Object
# Find genes in metabolic database
metab_sub <- read.csv(file="monogenic_sheet_6.csv")
subg <- metab_sub[metab_sub$Gene.name %in% gen,]
subg <- subg[!duplicated(subg$Gene.name), ]

metab_sub <- distinct(metab_sub)
metab_sub$X <- metab_sub$X %>% trimws(, which = "both")
subg$X <- subg$X %>% trimws(, which = "both")

# Count by both together for all genes in metab table
sub_sub <- metab_sub %>% count(X.1, name = "n")

# Count represented genes
case_sub <- subg %>% count(X.1)

# Combine genes and total table
enriched.table <- left_join(sub_sub,case_sub, by="X.1")
enriched.table$n.y[is.na(enriched.table$n.y)] <- 0
colnames(enriched.table) <- c("Pathway", "Total", "Cases")

# compute universe / subset totals once
total_genes  <- 1459
total_subset <- 30

enriched_stats <- enriched.table %>% rowwise() %>%
  mutate(expected_prop = Total / total_genes, 
         observed_prop = Cases / total_subset,
         enrichment = observed_prop / expected_prop,
         pval = {
           a <- Cases
           b <- Total - Cases
           c <- total_subset - Cases
           d <- total_genes - Total - c
           if(any(c(a,b,c,d) < 0)) NA_real_ else
             fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")$p.value }) %>% ungroup() %>%
  mutate(padj = p.adjust(pval, method = "BH"))

# Scale enrichment (z-score style)
enriched_stats <- enriched_stats %>%
  mutate(enrichment_scaled = scale(enrichment)[, 1],  sig_label = ifelse(padj < 0.05, "*", ""))

enriched <- enriched_stats %>% filter(Cases > 0)

ggplot(enriched, aes(x = Pathway, y = 1, fill = enrichment_scaled)) +
  geom_tile(color = "grey80") + geom_text(aes(label = sig_label), vjust = 0.5, size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Scaled enrichment") +
  labs(x = NULL, y = NULL, title = "Category enrichment") +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

subg_subset <- subg %>% dplyr::select(X, X.1, Gene.name)

