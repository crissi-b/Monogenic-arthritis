
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

# Import Gene Summary data from Data Object sheet 1 and reformat 
syndromes <- read.csv(file="monogenic_sheet_1.csv", header=F)
syndromes <- syndromes[3:nrow(syndromes),]
colnames(syndromes) <- syndromes[1, ]
syndromes <- syndromes[-1, ]

# Extract genes to include / exclude
excluded_genes <- syndromes %>% filter(exclude == "Y") %>% pull (Gene)
gen <- syndromes %>% filter(exclude == "N") %>% pull (Gene)
syndromes_inc <- syndromes %>% filter(exclude == "N")

#######################################################
#######################################################
#######################################################

# Data in Data Object sheet 2
case <- read.csv(file="monogenic_sheet_2.csv", header=F)
case <- case[3:nrow(case),]
colnames(case) <- case[1, ]
case <- case[-1, ]

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
