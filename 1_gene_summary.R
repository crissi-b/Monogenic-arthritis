
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(scales)
library(binom)
library(scales)
library(patchwork)
library(devtools)
library(disgenet2r)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)
library(ggthemes)
library(stats) 
library(pheatmap)

options(bitmapType='cairo')

#######################################################
#######################################################
#######################################################

syndromes <- read.csv(file="CSV/monogenic_sheet_1.csv", header=T)

# Extract genes to include / exclude
excluded_genes <- syndromes %>% filter(exclude == "Y") %>% pull (Gene)
gen <- syndromes %>% filter(exclude == "N") %>% pull (Gene)
syndromes_inc <- syndromes %>% filter(exclude == "N")

# Data in Data Object sheet 2
case <- read.csv(file="CSV/monogenic_sheet_a1.csv", header=T)

# Pull genes with only 1 case of arthritis identified
case_check <- case %>% filter(X2_cases == "N") %>% dplyr::select(-Cohort_references)

# Calculate estimated penetrance
case_check <- case_check %>% mutate(Arthritis_cases = as.numeric(Arthritis_cases), 
                                    Cohort_totals   = as.numeric(Cohort_totals),
                                    Age.of.Arthritis.Onset = as.numeric(Age.of.Arthritis.Onset))

# Compute proportions and 90% Wilson confidence intervals
ci_tbl <- binom::binom.confint(x = case_check$Arthritis_cases, n = case_check$Cohort_totals,
                               conf.level = 0.90, methods = "agresti-coull") %>% as_tibble() %>% 
                              dplyr::select(mean, lower, upper)

# Ensure bounds are between 0 and 1
ci_tbl <- ci_tbl %>% mutate(lower  = pmax(lower, 0.001), upper = pmax(upper, 1e-6), mean = pmax(mean, 1e-6))

case_ci <- bind_cols(case_check, ci_tbl) 

# Order genes by descending proportion
case_ci <- case_ci %>% mutate(Gene = fct_reorder(Gene, mean, .desc = TRUE))
age_col <- "Age.of.Arthritis.Onset" 

# Main plot 
p_main <- ggplot(case_ci, aes(x = Gene, y = mean)) +
  geom_point(size = 2) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "red") +  
  scale_y_log10(labels = percent_format(accuracy = 0.1), breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), limits = c(0.001, 1.1)) +
  labs(x = NULL, y = "Arthritis penetrance") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"), panel.grid = element_blank())

# Age tile plot with gene names, using the same Gene factor ordering as p_main (case_ci$Gene)
p_age <- ggplot(case_ci, aes(x = Gene, y = 1, fill = .data[[age_col]])) +
  geom_tile() + scale_fill_viridis_c(name = "Age of onset", option = "mako", direction = -1, na.value = "grey80") +
  labs(x = "Gene") + theme_minimal(base_size = 11) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 5)), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 9),
    legend.position = "bottom", panel.grid = element_blank())

p_main / p_age + plot_layout(heights = c(4, 0.8))

# Calculate significance of age
case_test <- case_ci %>% mutate(high_penetrance = lower > 0.01)

# Summarise groups for inspection
case_test %>% group_by(high_penetrance) %>% summarise(n = n(),
    mean_age = mean(Age.of.Arthritis.Onset, na.rm = TRUE), 
    median_age = median(Age.of.Arthritis.Onset, na.rm = TRUE),
    sd_age = sd(Age.of.Arthritis.Onset, na.rm = TRUE))

wilcox.test(Age.of.Arthritis.Onset ~ high_penetrance, data = case_test, exact = FALSE)

# This data is incorporated into the 'exclude' column in Table S5

#######################################################
#######################################################
#######################################################

# Plot numbers of genes by variant type and direction (GOF/LOF)

# Remove genes with low penetrance
syn_brief <- syndromes %>% filter(exclude == "N")
syn_brief <- subset(syn_brief, select=c(Variant.type,Direction), drop=FALSE)
syn_brief$Variant.type <- trimws(syn_brief$Variant.type, which = "both") 
counts <- syn_brief %>% count(Variant.type, Direction)

fill_colours <- c(
  "Biallelic.GOF" = "darkgreen",
  "Biallelic.LOF" = "lightgreen",
  "Biallelic / Het.GOF" = "deeppink4",
  "Biallelic / Het.LOF" = "pink",
  "Heterozygous.GOF" = "purple4",
  "Heterozygous.LOF" = "orchid",
  "X-linked.GOF" = "royalblue3",
  "X-linked.LOF" = "lightblue")

ggplot(counts, aes(x = reorder(Variant.type, n), y = n, fill = interaction(Variant.type, Direction))) +
  geom_bar(stat = "identity") + coord_flip() + 
  scale_fill_manual(name = "Variant Type & Direction", values = fill_colours) +
  labs(x = "Variant type", y = "Number of genes") + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Caclulate GOF enrichment- note that some genes in 'Bialellic / Heterozygous' group not included in this count
counts2 <- counts %>% filter(Variant.type %in% c("Biallelic", "Heterozygous")) %>%
  dplyr::select(Variant.type, Direction, n) %>%
  tidyr::pivot_wider(names_from = Direction, values_from = n, values_fill = 0)

counts2
m <- as.matrix(counts2[, c("GOF", "LOF")])
rownames(m) <- counts2$Variant.type
fisher.test(m)

#######################################################
#######################################################
#######################################################

# DisGeNet representation of genes

# Requires api_key from registration at https://disgenet.com - free academic accounts available
api_key <- # Insert your key here
Sys.setenv(DISGENET_API_KEY= api_key)

# Ensure aliases are captured of first 90 genes, disgenet only allows searches of 100 genes at a time
t90 <- syndromes_inc[1:90,] %>% pull(Gene)

# Convert genes to symbols
alias_to_symbol <- AnnotationDbi::select(org.Hs.eg.db,keys = t90, keytype = "ALIAS",columns = c("SYMBOL", "ENTREZID"))

# Check duplicates
dup <- alias_to_symbol[duplicated(alias_to_symbol$ALIAS), ] 

# Inspect duplicated and remove the wrong gene aliases
remove <- as.character(c(6871, 3183, 1671, 2194, 9027, 6677, 257641, 22883, 22883, 55893))
alias_to_symbol <- alias_to_symbol %>% filter( ! ENTREZID %in% remove )

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(alias_to_symbol$SYMBOL,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

entrez_ids <- gene_entrez$ENTREZID

res <- gene2disease(gene = entrez_ids, 
                    vocabulary = "ENTREZ",
                    database = "CURATED", 
                    score = c( 0.1,1))

qresult <- res@qresult
head(qresult)

# Convert list columns to character
qresult_fixed <- qresult %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(unlist(x), collapse = "; "))))

# Repeat for remaining GOF genes
# Select remaining genes
rest <- syndromes_inc[91:length(rownames(syndromes_inc)),] %>% pull(Gene)

alias_to_symbol <- AnnotationDbi::select(org.Hs.eg.db,keys = rest, keytype = "ALIAS",columns = c("SYMBOL", "ENTREZID"))
dup <- alias_to_symbol[duplicated(alias_to_symbol$ALIAS), ] # Inspect duplications

# Inspect duplicated and remove the wrong gene aliases
remove <- as.character(c(9997))
alias_to_symbol <- alias_to_symbol %>% filter( ! ENTREZID %in% remove )

# Convert gene symbols to Entrezgene_entrez# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(alias_to_symbol$SYMBOL,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

entrez_ids <- gene_entrez$ENTREZID

res2 <- gene2disease(gene = entrez_ids, 
                     vocabulary = "ENTREZ",
                     database = "CURATED", 
                     score = c( 0.1,1))

qresult2 <- res2@qresult
head(qresult2)[1:5,1:5]

# Convert list columns to character
qresult_fixed2 <- qresult2 %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(unlist(x), collapse = "; "))))

# Combine datasets of LOF and GOF genes
combined_dgn <- rbind(qresult_fixed, qresult_fixed2)
counts <- combined_dgn %>% count(disease_name)

# Select immune-mediated inflammatory diseases
disease_type <- c('C0003873', 'C3890205', 'C0157917', 'C0157916', 'C0003872', 'C0409667', 'C3495559', 
                  'C2931171', 'C0003089', 'C0156147', 'C0009324', 'C0375359', 'C0011854', 'C0024141', 
                  'C0026769', 'C0033860', 'C1527336', "C0157918", "C0004943")

# Filter for genes and selected disease types
dgn_select <- combined_dgn %>% filter(gene_symbol %in% gen, diseaseUMLSCUI %in% disease_type) %>%
  distinct(gene_symbol, disease_name)

# Correct names and combine categories
dgn_select$disease_name[dgn_select$disease_name == "Ankylosing spondylitis and other inflammatory spondylopathies"] <- "Ankylosing spondylitis"
dgn_select$disease_name[dgn_select$disease_name == "Diabetes Mellitus, Insulin-Dependent"] <- "Insulin-dependent diabetes mellitus"
dgn_select$disease_name[dgn_select$disease_name == "Lupus Erythematosus, Systemic"] <- "Systemic lupus erythematosus"
dgn_select$disease_name[dgn_select$disease_name == "Chronic left-sided ulcerative colitis"] <- "Ulcerative colitis"
dgn_select$disease_name[dgn_select$disease_name == "Crohn's disease of large bowel"] <- "Crohn's disease"
dgn_select$disease_name[dgn_select$disease_name == "Polyarticular Juvenile Idiopathic Arthritis, Rheumatoid Factor Negative"] <- "Juvenile polyarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Acute polyarticular juvenile rheumatoid arthritis"] <- "Juvenile polyarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Juvenile Chronic Polyarthritis"] <- "Juvenile polyarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Pauciarticular juvenile rheumatoid arthritis"] <- "Juvenile oligoarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Juvenile arthritis"] <- "Juvenile oligoarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Juvenile pauciarticular chronic arthritis"] <- "Juvenile oligoarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Monoarticular juvenile rheumatoid arthritis"] <- "Juvenile oligoarticular arthritis"
dgn_select$disease_name[dgn_select$disease_name == "Arthritis, Psoriatic"] <- "Psoriatic arthritis"
dgn_select$disease_name <- str_to_sentence(dgn_select$disease_name)

unique(dgn_select$disease_name)

# Calculate number of genes per disease
gene_counts <- dgn_select %>% distinct(disease_name, gene_symbol) %>% count(disease_name)

# Horizontal bar plot for gene counts
ggplot(gene_counts, aes(y = reorder(disease_name, n), x = n)) +
  geom_col(fill = "grey50") +
  theme_minimal() + scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of genes", y = NULL) + 
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text = element_text(color="black"))
                                                
# Collapse gene symbols by disease_name
disease_gene_table <- dgn_select %>%
  distinct(disease_name, gene_symbol) %>%
  group_by(disease_name) %>%
  summarise(genes = paste(unique(gene_symbol), collapse = ", ")) %>%
  arrange(disease_name)

#######################################################
#######################################################
#######################################################

# Plot median age of onset for the gene defects
age <- case

# Split out ages by comma
colnames(age) <- c("Gene", "Age")
age_long <- age[,c("Gene", "Age")] %>% separate_rows(Age, sep = ",\\s*") %>%        
  mutate(Age = as.numeric(str_trim(Age)))   

# Remove genes with low penetrance
age_long <- age_long[!(age_long$Gene %in% excluded_genes),]

# Add medians
medians <- age_long %>% group_by(Gene) %>%
  summarise(median_age = median(Age, na.rm = TRUE)) %>% arrange(median_age)

medians <- medians[!(medians$Gene==""), ]

age_long %>% mutate(Gene = factor(Gene, levels = medians$Gene)) %>% # order by median
  ggplot(aes(x = Gene, y = Age)) +
  geom_point(alpha = 0.7, size=2, color = "grey50") +
  geom_line(data = medians, aes(x = Gene, y = median_age, group = 1),
            color = "red", linetype = "dashed", size=0.8) + 
  geom_hline(aes(yintercept = 17), color="blue", linetype="dashed") +
  labs(x = "Gene", y = "Age of arthritis onset (years)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                     text = element_text(size = 10), panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                     axis.text = element_text(color="black"))

# Calculate proportion below 3 year sand 17 years
medians %>%
  summarise(
    total_genes = n(),
    n_below_3 = sum(median_age < 3, na.rm = TRUE),
    prop_below_3 = n_below_3 / total_genes,
    n_below_17 = sum(median_age < 17, na.rm = TRUE),
    prop_below_17 = n_below_17 / total_genes)

#######################################################
#######################################################
#######################################################

# Cumulative LOF from the gnomad constraint database
direct <- syndromes_inc[,c("Gene", "Direction"),drop=FALSE]
lof <- direct %>% filter(Direction == "LOF") %>% pull(Gene)

# Load gnomad- URL for the gnomAD v2.1 constraint metrics file
dest_file <- "Datasets/gnomad.v2.1.1.constraint.txt.bgz" # Filename for the downloaded file
url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
download.file(url, destfile = dest_file, mode = "wb")
constraint <- fread(dest_file)

# Filter for LOF genes 
lof_constraint <- constraint[gene %in% lof]
setdiff(lof, lof_constraint$gene) # check missing

# Add synonyms not captured in the gene list: OTULIN = FAM105B, ADA2 = CECR1, RLTPR = CARMIL2
lof <- c(lof, "CECR1", "FAM105B", "RLTPR")
lof_constraint <- constraint[gene %in% lof]

old <- c("CECR1", "FAM105B", "RLTPR")
new <- c("ADA2",  "OTULIN", "CARMIL2")
lof_constraint$gene[match(old, lof_constraint$gene)] <- new

# replace where matches are found
ggplot(lof_constraint, aes(x = reorder(gene, classic_caf), y = classic_caf)) +
  geom_point(alpha = 1, color = "grey30") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal(base_size = 8) +  geom_hline(yintercept = 1e-2, linetype = "dashed", color = "grey") +
  labs(x = "", y = "Cumulative allele frequency of predicted LOF variants (log10 scale)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Add X-linked annotation to plot
lof_type <- syndromes %>% filter(Direction == "LOF") %>% dplyr::select(c(Gene, Variant.type))

# Order genes by classic_caf 
levels_order <- lof_constraint %>% arrange(classic_caf) %>% pull(gene) %>% as.character() 
all_genes <- tibble(Gene = levels_order)

# add ordered factor to lof_constraint for plotting
lof_constraint <- lof_constraint %>% mutate(gene_f = factor(gene, levels = levels_order))

# prepare lof_type annotation, ensuring same genes and same factor order
lof_type <- lof_type %>% filter(Gene %in% levels_order) %>% mutate(Gene = factor(Gene, levels = levels_order)) %>%
  mutate(tile_fill = if_else(Variant.type == "X-linked", "steelblue", "white"))

# Combine to order
lof_type <- all_genes %>% left_join(lof_type, by = "Gene") %>%
  mutate(Gene = factor(Gene, levels = levels_order),
    tile_fill = if_else(Variant.type == "X-linked", "steelblue", "white", missing = "white"))

p_main <- ggplot(lof_constraint, aes(x = gene_f, y = classic_caf)) +
  geom_point(alpha = 1, color = "grey30") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal(base_size = 8) +  geom_hline(yintercept = 1e-2, linetype = "dashed", color = "grey") +
  labs(x = "", y = "Cumulative allele frequency of predicted LOF variants (log10 scale)") +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text = element_text(color="black"))

# Tile annotation row with gene names as x-axis labels
p_tile <- ggplot(lof_type, aes(x = Gene, y = 1, fill = tile_fill)) +
  geom_tile(color="grey50") + scale_fill_identity() + labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
    panel.grid = element_blank(),axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    legend.position = "none",plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
    axis.text = element_text(color="black"))

combined <- p_main / p_tile + plot_layout(heights = c(4, 0.4))
combined

caf_values <- subset(lof_constraint, select= c("gene", "classic_caf", "oe_lof_upper"), drop=FALSE)

save.image(file="gene-summary-data.RData")

#######################################################
#######################################################
#######################################################

# To determine representation in immune categories of IUIS Inborn Errors of Immunity, first find genes represented in IUIS tables
# Available from monogenic sheet 3 in Data_Objects
iuis <- read.csv(file="CSV/monogenic_sheet_3.csv", header=T)

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

# Available from monogenic sheet 3 in Data_Objects
arthritis.iuis <- read.csv(file="CSV/monogenic_sheet_4.csv")

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

# Determine representation in metabolic classification (ICIMD)

# Metabolic gene summary  (Ferreira 2021), 
# Read in full classification table, data available in sheet 5 of the Data Object
metab <- read.csv(file="CSV/monogenic_sheet_5.csv")

# Find genes in metabolic database
metab_gen <- metab[metab$Gene.name %in% gen,]
metab_gen <- metab_gen[!duplicated(metab_gen$Gene.name), ]
# write.csv(metab_gen, file="CSV/metabolic_arthritis-genes.csv")

metab <- distinct(metab)
metab$X <- metab$X %>% trimws(, which = "both")
metab_gen$X <- metab_gen$X %>% trimws(, which = "both")

syndromes <- syndromes %>%
  left_join(
    metab %>% dplyr::select(Gene.name, X),
    by = c("Gene_D" = "Gene.name"))

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
metab_sub <- read.csv(file="CSV/monogenic_sheet_6.csv")
subg <- metab_sub[metab_sub$Gene.name %in% gen,]
subg <- subg[!duplicated(subg$Gene.name), ]

metab_sub <- distinct(metab_sub)
metab_sub$X <- metab_sub$X %>% trimws(, which = "both")
subg$X <- subg$X %>% trimws(, which = "both")

syndromes <- syndromes %>%
  left_join(
    subg %>% dplyr::select(Gene.name, X.1),
    by = c("Gene_D" = "Gene.name"))

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
