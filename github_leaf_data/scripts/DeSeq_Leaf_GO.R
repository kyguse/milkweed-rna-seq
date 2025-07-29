#####Gene Enrichment Analysis for LEAF data (07/29/25) - This code was supported by CHAT GPT when some errors occurred

#### Load libraries
library(topGO)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DESeq2)

###Step 1: Extract DeSeq2 Results
res_leaf <- results(dds_leaf, contrast = c("species", "speciosa", "syriaca"))

####Step 2: Filter and create geneList vector (1 = significant, 0 = not)
geneList <- as.integer(rownames(res_leaf) %in% rownames(res_leaf[res_leaf$padj < 0.05 & !is.na(res_leaf$padj), ]))
names(geneList) <- rownames(res_leaf)

###Step 3: Get the correct annotations
library(dplyr)
library(tidyr)

# Unite the GO sources into one column
annot_go <- annot %>%
  dplyr::select(X.gene_id, gene_ontology_Pfam, gene_ontology_BLASTX, gene_ontology_BLASTP) %>%
  unite("go_terms_all", gene_ontology_Pfam, gene_ontology_BLASTX, gene_ontology_BLASTP, sep = ";", na.rm = TRUE)

# Separate into long format by semicolon
go_long <- annot_go %>%
  separate_rows(go_terms_all, sep = ";") %>%
  rename(go_terms_short = go_terms_all) %>%
  filter(go_terms_short != "")

go_long <- go_long %>%
  mutate(go_terms_short = sub("\\^.*", "", go_terms_short)) %>%
  filter(go_terms_short != "")

go_map <- go_long %>%
  group_by(X.gene_id) %>%
  summarise(go_list = list(unique(go_terms_short)), .groups = "drop") %>%
  deframe()

###Step 4: Create named list of GO terms for each gene
go_ids_only <- go_long %>%
  group_by(X.gene_id) %>%
  summarise(go_list = list(unique(go_terms_short)), .groups = "drop") %>%
  deframe()

###Step 5: Define gene universe and overlap with GO terms
shared_genes <- intersect(names(geneList), names(go_ids_only))
geneList <- geneList[shared_genes]
go_ids_only <- go_ids_only[shared_genes]

###Step 6: Define DE genes by species
speciosa_up <- rownames(res_leaf[
  res_leaf$log2FoldChange > 1 & 
    res_leaf$padj < 0.05 & 
    !is.na(res_leaf$padj) & 
    !is.na(res_leaf$log2FoldChange),
])

### Do the same for syriaca
syriaca_up <- rownames(res_leaf[
  res_leaf$log2FoldChange < -1 & 
    res_leaf$padj < 0.05 & 
    !is.na(res_leaf$padj) & 
    !is.na(res_leaf$log2FoldChange),
])

all_genes <- rownames(res_leaf)

###Step 7: Creating Binary gene list for topGO
geneList_speciosa <- factor(as.integer(all_genes %in% speciosa_up))
names(geneList_speciosa) <- all_genes

geneList_syriaca <- factor(as.integer(all_genes %in% syriaca_up))
names(geneList_syriaca) <- all_genes

### Step 8: Filter GO map to match gene universe
go_map <- go_ids_only[names(go_ids_only) %in% all_genes]

# Step 9: Function to run TopGO
display <- function(geneList, go_map, ontology = "BP") {
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                geneSelectionFun = function(x) x == 1,
                annot = annFUN.gene2GO,
                gene2GO = go_map)
  
  result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  GenTable(GOdata, classicFisher = result, topNodes = 20)
}

###Step 10: Run Enrichment
speciosa_go <- display(geneList_speciosa, go_map)
syriaca_go  <- display(geneList_syriaca,  go_map)

###Step 11: Prepare Data for Barplot
speciosa_go$Species <- "Speciosa"
syriaca_go$Species  <- "Syriaca"

combined_go <- bind_rows(speciosa_go, syriaca_go)
combined_go$log10pval <- -log10(as.numeric(combined_go$classicFisher))

###Step 12: Keep Top 20 by species
top_combined <- combined_go %>%
  group_by(Species) %>%
  arrange(desc(log10pval)) %>%
  slice_head(n = 20)

##Step 13: Create barplot
png("NEW_LEAF_GO_side-by-side_plot.png", height = 6, width = 8, units = 'in', res = 600)

ggplot(top_combined, aes(x = reorder(Term, log10pval), y = log10pval, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  coord_flip() +
  labs(
    title = "Top 20 Enriched GO Terms per Species (Leaf)",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  scale_fill_manual(values = c("Speciosa" = "#FEC44F", "Syriaca" = "deeppink")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 6),
    plot.margin = margin(10, 10, 10, 20),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

dev.off()

#####Now let's do the same looking at Molecular Function:

###Change the ontology to "MF" for molecular function
speciosa_go_MF <- display(geneList_speciosa, go_map, ontology = "MF")
syriaca_go_MF  <- display(geneList_syriaca,  go_map, ontology = "MF")

####Add species label and combine
speciosa_go_MF$Species <- "Speciosa"
syriaca_go_MF$Species  <- "Syriaca"

combined_go_MF <- bind_rows(speciosa_go_MF, syriaca_go_MF)
combined_go_MF$log10pval <- -log10(as.numeric(combined_go_MF$classicFisher))

###Keep top 20 by species
top_combined_MF <- combined_go_MF %>%
  group_by(Species) %>%
  arrange(desc(log10pval)) %>%
  slice_head(n = 20)

###Barplot Molecular Function GO Terms
png("NEW_Leaf_GO_MF_side-by-side_plot.png", height = 6, width = 8, units = 'in', res = 600)

ggplot(top_combined_MF, aes(x = reorder(Term, log10pval), y = log10pval, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  coord_flip() +
  labs(
    title = "Top 20 Enriched GO Terms (Molecular Function, Root)",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  scale_fill_manual(values = c("Speciosa" = "plum2", "Syriaca" = "turquoise")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 6),
    plot.margin = margin(10, 10, 10, 20),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

dev.off()

