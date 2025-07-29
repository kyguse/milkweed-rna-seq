#####Gene Enrichment Analysis (07/29/25) - This code was supported by CHAT GPT when some errors occurred

#### Load libraries
library(topGO)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DESeq2)

###Step 1: Extract DeSeq2 Results
res_root <- results(dds_root, contrast = c("species", "speciosa", "syriaca"))

####Step 2: Filter and create geneList vector (1 = significant, 0 = not)
geneList <- as.integer(rownames(res_root) %in% rownames(res_root[res_root$padj < 0.05 & !is.na(res_root$padj), ]))
names(geneList) <- rownames(res_root)

###Step 3: Clean and extract gene-to-GO mappings from annotated dataframe (go_long)
go_long <- go_long %>%
  mutate(go_terms_short = sub("\^.*", "", go_terms_short)) %>%
  filter(go_terms_short != "")

###Step 4: Create named list of GO terms for each gene
go_ids_only <- go_long %>%
  group_by(gene_id) %>%
  summarise(go_list = list(unique(go_terms_short)), .groups = "drop") %>%
  deframe()

###Step 5: Define gene universe and overlap with GO terms
shared_genes <- intersect(names(geneList), names(go_ids_only))
geneList <- geneList[shared_genes]
go_ids_only <- go_ids_only[shared_genes]

###Step 6: Define DE genes by species
speciosa_up <- rownames(res_root[
  res_root$log2FoldChange > 1 & 
    res_root$padj < 0.05 & 
    !is.na(res_root$padj) & 
    !is.na(res_root$log2FoldChange),
])

### Do the same for syriaca
syriaca_up <- rownames(res_root[
  res_root$log2FoldChange < -1 & 
    res_root$padj < 0.05 & 
    !is.na(res_root$padj) & 
    !is.na(res_root$log2FoldChange),
])

all_genes <- rownames(res_root)

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
png("Root_GO_side-by-side_plot.png", height = 6, width = 8, units = 'in', res = 600)

ggplot(top_combined, aes(x = reorder(Term, log10pval), y = log10pval, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  coord_flip() +
  labs(
    title = "Top 20 Enriched GO Terms per Species (Root)",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  scale_fill_manual(values = c("Speciosa" = "salmon2", "Syriaca" = "salmon4")) +
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
png("Root_GO_MF_side-by-side_plot.png", height = 6, width = 8, units = 'in', res = 600)

ggplot(top_combined_MF, aes(x = reorder(Term, log10pval), y = log10pval, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  coord_flip() +
  labs(
    title = "Top 20 Enriched GO Terms (Molecular Function, Root)",
    x = "GO Term",
    y = "-log10(p-value)"
  ) +
  scale_fill_manual(values = c("Speciosa" = "mediumpurple4", "Syriaca" = "lightgoldenrod4")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 6),
    plot.margin = margin(10, 10, 10, 20),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

dev.off()

