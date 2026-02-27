
library(Seurat)
library(pheatmap)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(scales)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(dplyr)
library(tibble)
library(KEGGREST)
library(tidyr)
library(stringr)
library(purrr)

formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
colurs <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), brewer.pal(8, "Paired"))

syndromes <- read.csv(file="monogenic_sheet_1.csv", header=F, blank.lines.skip = TRUE)
syndromes <- syndromes[3:nrow(syndromes),]
colnames(syndromes) <- syndromes[1, ]
syndromes <- syndromes[-1, ]

#######################################################
#######################################################
#######################################################

# KEGG, GO, reactome, hallmark pathways
# Ensure aliases are captured
alias_to_symbol <- AnnotationDbi::select(org.Hs.eg.db,keys = gen,keytype = "ALIAS",columns = c("SYMBOL", "ENTREZID"))
check <- alias_to_symbol[duplicated(alias_to_symbol$ALIAS), ] # Deduplicate

# Remove erroneous gene names
alias_to_symbol <- alias_to_symbol[!(alias_to_symbol$SYMBOL %in% c("HNRNPC", "DEFA6", "TADA2A", "FASN", "GBAP1","NAT8", "SPAM1", "NPCA1", "CLSTN1", "PPY", "ZNF395", "SCO2")),]
setdiff(alias_to_symbol$SYMBOL, gen)

# Convert gene symbols to Entrezgene_entrez# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(alias_to_symbol$SYMBOL,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

entrez_ids <- unique(gene_entrez$ENTREZID)

#######################################################

# GO Enrichment
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

go_df <- as.data.frame(ego) %>% mutate(Database = "GO")

#######################################################

# KEGG Enrichment
kegg_links <- keggLink("pathway", "hsa")  # Named vector: names = genes, values = pathways

kegg_df <- data.frame(
  ENTREZID = sub("hsa:", "", names(kegg_links)),
  pathway_id = sub("path:", "", kegg_links),
  stringsAsFactors = FALSE
)

term2gene_kegg <- kegg_df[, c("pathway_id", "ENTREZID")]

pathway_info <- keggList("pathway", "hsa")
term2name_kegg <- data.frame(
  pathway_id = sub("path:", "", names(pathway_info)),
  pathway_name = as.character(pathway_info),
  stringsAsFactors = FALSE
)

ekegg <- enricher(gene = entrez_ids,
                  TERM2GENE = term2gene_kegg,
                  TERM2NAME = term2name_kegg,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

kegg_df <- as.data.frame(ekegg) %>%
  mutate(Database = "KEGG")

kegg_df$Description <- gsub("- Homo sapiens \\(human\\)", "", as.character(kegg_df$Description))

# Function to convert Entrez string to comma-separated symbol string
convert_entrez_to_symbols <- function(entrez_str) {
  entrez_ids <- unlist(strsplit(entrez_str, "/"))
  symbols_df <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)   # Convert to symbols
  symbols_ordered <- symbols_df$SYMBOL[match(entrez_ids, symbols_df$ENTREZID)]
  paste(symbols_ordered, collapse = "/")
}

# Apply function to each row
kegg_df$geneID <- sapply(kegg_df$geneID, convert_entrez_to_symbols)

#######################################################

# Reactome Enrichment
ereact <- enrichPathway(gene = entrez_ids,
                        organism = "human",
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        readable = TRUE)

reactome_df <- as.data.frame(ereact) %>%
  mutate(Database = "Reactome")

#######################################################

# Hallmark Enrichment (via MSigDB)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

ehallmark <- enricher(
  gene = entrez_ids,
  TERM2GENE = hallmark_sets[, c("gs_name", "entrez_gene")])

hallmark_df <- as.data.frame(ehallmark) %>%
  mutate(Database = "Hallmark")

# Apply function to each row
hallmark_df$geneID <- sapply(hallmark_df$geneID, convert_entrez_to_symbols)

# Reformat Hallmark pathways to be the same
hallmark_df$Description <- gsub("HALLMARK_", "", hallmark_df$Description)
hallmark_df$Description <- str_to_sentence(hallmark_df$Description)
hallmark_df$Description <- gsub("_", " ", hallmark_df$Description)

#######################################################

# Combine all into one dataframe
all_enrichment <- bind_rows(go_df, kegg_df, reactome_df, hallmark_df)
all_enrichment$Description <- str_to_sentence(all_enrichment$Description)

#######################################################
#######################################################
#######################################################

# Visualise top expressing pathways

# Identify more specific pathways, using BgRatio- extract the number before the "/"
filtered <- all_enrichment %>% mutate(BgNum = as.numeric(str_extract(BgRatio, "^[0-9]+"))) %>% filter(BgNum < 300)
uni_gene <- filtered %>% distinct(geneID, .keep_all = TRUE) # Unique combinations of genes to reduce redundancy

top20 <- uni_gene %>% arrange(qvalue) %>% slice_head(n = 30)
top20 <- top20[!grepl("Shigell|covid|Epstein|Primary", top20$Description),] # Remove disease terms

top20 <- top20 %>%
   mutate(GeneRatio_num = as.numeric(sapply(str_split(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
    GeneRatio_numer = as.numeric(sapply(str_split(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]))),
    GeneRatio_denom = as.numeric(sapply(str_split(as.character(GeneRatio), "/"), function(x) as.numeric(x[2]))),
    Count = as.numeric(Count),
    negLogQ = -log10(as.numeric(qvalue))) 

# Remove redundant terms
top20 <- top20[!(top20$Description %in% c("Primary immunodeficiency ",
                                          "Negative regulation of defense response",
                            "Positive regulation of nf-kappab transcription factor activity",
                            "Nucleotide-binding oligomerization domain containing 2 signaling pathway",
                            "Non-canonical nf-kappab signal transduction", 
                            "Positive regulation of nf-kappab transcription factor activity",
                            "Regulation of canonical nf-kappab signal transduction")),]

top20 <- top20[!(top20$Description %in% c("Cytoplasmic PRR signaling pathway",
                                          "Nucleotide-binding domain, leucine rich repeat containing receptor (nlr) signaling pathways",
                                          "Regulation of B cell activation",
                                          "Regulation of response to cytokine stimulus")),]

top20$Description <- gsub("Nf-kappa b", "NFkB", top20$Description)
top20$Description <- gsub("Nod", "NOD", top20$Description)
top20$Description <- gsub("Interleukin", "IL", top20$Description)
top20$Description <- gsub("b cell", "B cell", top20$Description)
top20$Description <- gsub("Tumor necrosis factor", "TNF", top20$Description)
top20$Description <- gsub("dna", "DNA", top20$Description)
top20$Description <- gsub("pattern recognition receptor", "PRR", top20$Description, ignore.case=T)
top20$Description <- gsub("transcription factor", "TF", top20$Description, ignore.case=T)
top20$Description <- gsub("TNFR1", "TNFR1", top20$Description, ignore.case=T)

top20 <- top20 %>% mutate(Description_plot = fct_reorder(Description, GeneRatio_num)) 
top20_order <- rev(levels(top20$Description_plot))

ggplot(top20, aes(x = GeneRatio_num, y = Description_plot)) +
  geom_point(aes(size = Count, color = negLogQ), alpha = 0.9) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), name = "Fraction of genes in pathway") +
  scale_size_continuous(name = "Count (genes)") + scale_color_viridis(name = "-log10(qvalue)", option = "D") +
  labs(title = "Top enriched pathways", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 12), plot.title = element_text(face = "bold"), legend.position = "right",
        panel.grid=element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
        
#######################################################
#######################################################
#######################################################

# Age related pathways
syndromes_tb <- syndromes %>% filter(exclude == "N")
syn_long <- syndromes_tb %>% dplyr::select(Gene_D, age_group)
colnames(syn_long) <- c("Gene_D", "Phenotype")
syn_long$Gene_D <- gsub("_LOF|_GOF", "", syn_long$Gene_D)
syn_long <- distinct(syn_long)

# Remove LOF or GOF
syn_long$Gene_D <- gsub("_LOF|_GOF", "", syn_long$Gene_D)

# For each phenotype, list its genes
phenotype_genes <- syn_long %>%
  group_by(Phenotype) %>%
  summarise(GeneSet = list(unique(Gene_D)))

# Convert pathway gene lists (geneID) to character vectors
all_enrichment_clean <- all_enrichment %>%
  mutate(geneID = gsub("STING1", "TMEM173", geneID), GeneList = str_split(geneID, "/")) %>%
  mutate(geneID = gsub("RIGI", "DDX58", geneID), GeneList = str_split(geneID, "/"))

# Increase specificity by excluding pathways with high number of background genes
all_enrichment_clean <- all_enrichment_clean %>% mutate(BgNum = as.numeric(str_extract(BgRatio, "^[0-9]+"))) %>% filter(BgNum < 200)

# For each pathway, get its genes and ID
pathway_genes <- all_enrichment_clean %>%
  dplyr::select(ID, Description, GeneList, Database)

# Define a function to compute overlap between phenotype and pathway gene sets
compute_overlap <- function(pheno, pheno_genes, path_id, path_name, path_genes, db) {
  overlap <- intersect(pheno_genes, path_genes)
  data.frame(Phenotype = pheno, PathwayID = path_id, Pathway = path_name,
             Database = db, OverlapCount = length(overlap), OverlapGenes = paste(overlap, collapse = ","))
}

# Cross all phenotype/pathway pairs and compute overlaps
overlap_df <- map_dfr(1:nrow(phenotype_genes), function(i) {
  pheno <- phenotype_genes$Phenotype[i]
  pheno_genes <- phenotype_genes$GeneSet[[i]]
  
  map_dfr(1:nrow(pathway_genes), function(j) {
    path_id <- pathway_genes$ID[j]
    path_name <- pathway_genes$Description[j]
    path_genes <- pathway_genes$GeneList[[j]]
    db <- pathway_genes$Database[j]
    
    compute_overlap(pheno, pheno_genes, path_id, path_name, path_genes, db)
  })
})

# Keep only overlaps with at least 2 genes
overlap_summary <- overlap_df %>%
  filter(OverlapCount > 1) %>%
  arrange(desc(OverlapCount))

# Top 2 per phenotype
top_pw <- overlap_summary %>%
  group_by(Phenotype) %>%
  arrange(desc(OverlapCount)) %>%
  slice_max(n = 2, order_by=OverlapCount, with_ties=TRUE) %>%   
  ungroup()

top_pw <- top_pw[!is.na(top_pw$Phenotype),]

# Ensure only pathways with highest overlap count are selected if featuring twice
top_pw2 <- top_pw %>% group_by(Pathway) %>% slice_max(order_by = OverlapCount, n = 1, with_ties = FALSE) %>% ungroup()

# Select rows of relevant pathways
included.pathways <- top_pw$PathwayID
combined_pw <- all_enrichment[all_enrichment$ID %in% included.pathways,]
age_pw <- merge(top_pw2, combined_pw, by.x="PathwayID", by.y="ID")
age_pw <- age_pw %>% dplyr::select(-c(ONTOLOGY, pvalue, Pathway, GeneRatio, p.adjust, Database.y))

age_pw$geneID <- gsub("STING1", "TMEM173", age_pw$geneID)
age_pw$geneID <- gsub("RIGI", "DDX58", age_pw$geneID)
gene_sets <- str_split(age_pw$geneID, "/") # Extract gene sets from the top pathways

pf <- age_pw %>% mutate(gene = str_split(geneID, "/")) %>% unnest(gene) %>%
  mutate(gene = str_trim(gene)) %>%
  left_join(syn_long %>% mutate(Gene_D = str_trim(Gene_D)), by = c("gene" = "Gene_D"))

# compute counts and proportions per pathway x phenotype
pf_sum <- pf %>% group_by(Description, Phenotype.y) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Description) %>% mutate(prop = n / sum(n)) %>%
  ungroup() 

pf_sum <- pf_sum[!is.na(pf_sum$Phenotype.y),]

order_m <- unique(pf_sum$Phenotype.y)
order_m <- order_m[c(4,3,2,1)]
pf_sum$Phenotype.y <- factor(pf_sum$Phenotype.y, levels=c(order_m))
age_colours <- c("0-2" = "blue", "3-6" = "lightblue1", "7-16" = "lemonchiffon2", ">17" = "cornsilk3")

ordered <- c("Regulation of response to cytokine stimulus", "Regulation of cytokine-mediated signaling pathway", "B cell differentiation","Regulation of b cell activation", "Myeloid cell homeostasis", "Leukocyte apoptotic process",                               
          "Regulation of leukocyte apoptotic process", "Osteoclast differentiation ", "Nod-like receptor signaling pathway ","Cytoplasmic pattern recognition receptor signaling pathway", "Immune system development", "Nf-kappa b signaling pathway " )

pf_sum$Description <- factor(pf_sum$Description, levels=ordered)

# Stacked proportion bar (one bar per pathway)
ggplot(pf_sum, aes(y = Description, x = n, fill = Phenotype.y)) +
  geom_col(position="dodge") + scale_fill_manual(values=age_colours) + 
  theme_bw() + guides(fill = guide_legend(reverse = TRUE)) +
  theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
  labs(x = NULL, y = "", fill = "Onset age") + 
  scale_x_continuous(breaks = scales::pretty_breaks()) + scale_y_discrete(limits = rev)


#######################################################
#######################################################
#######################################################

# stats
library(broom)
library(tidyr)

table(syndromes_tb$age_group)
age_totals <- c(">17" = 13, "0-2" = 27, "3-6" = 30, "7-16" = 34)

pf_sum$Phenotype.y <- as.character(pf_sum$Phenotype.y)
colnames(pf_sum) <- c("pathway", "age_group", "n", "prob")

N <- sum(age_totals)

pathway_sizes <- aggregate(n ~ pathway, data = df, sum)
names(pathway_sizes)[2] <- "M"

results <- data.frame()

# To determine odds ratio > 1.5
for(i in 1:nrow(pf_sum)){
  pathway   <- pf_sum$pathway[i]
  age       <- pf_sum$age_group[i]
  k         <- pf_sum$n[i]                    # genes in pathway & age
  K         <- age_totals[age]             # total genes in age group
  M         <- pathway_sizes$M[pathway_sizes$pathway == pathway]
  N_total   <- N
  
  # 2x2 contingency table
  mat <- matrix(c(k, K - k, M - k, N_total - K - (M - k)), nrow = 2)
  
  ft <- fisher.test(mat)
  
  oddsr <- rbind(results, data.frame(pathway = pathway, age_group = age,
      k = k, K = K, M = M, odds_ratio = ft$estimate))
}

