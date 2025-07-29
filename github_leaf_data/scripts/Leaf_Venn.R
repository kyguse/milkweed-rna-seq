####Venn Diagram of Leaf Data - Recreated on 07/29/25

library(dplyr)
library(purrr)
library(VennDiagram)

# Recreate binary gene lists
geneList_speciosa <- factor(as.integer(rownames(dds_leaf) %in% speciosa_up))
names(geneList_speciosa) <- rownames(dds_leaf)

geneList_syriaca <- factor(as.integer(rownames(dds_leaf) %in% syriaca_up))
names(geneList_syriaca) <- rownames(dds_leaf)

# Create topGOdata objects
library(topGO)

GOdata_speciosa <- new("topGOdata",
                       ontology = "BP",
                       allGenes = geneList_speciosa,
                       geneSelectionFun = function(x) x == 1,
                       annot = annFUN.gene2GO,
                       gene2GO = go_map)

GOdata_syriaca <- new("topGOdata",
                      ontology = "BP",
                      allGenes = geneList_syriaca,
                      geneSelectionFun = function(x) x == 1,
                      annot = annFUN.gene2GO,
                      gene2GO = go_map)

# Run enrichment tests
result_speciosa <- runTest(GOdata_speciosa, algorithm = "classic", statistic = "fisher")
result_syriaca  <- runTest(GOdata_syriaca,  algorithm = "classic", statistic = "fisher")

# Generate tables
speciosa_go <- GenTable(GOdata_speciosa, classicFisher = result_speciosa, topNodes = 100)
syriaca_go  <- GenTable(GOdata_syriaca,  classicFisher = result_syriaca,  topNodes = 100)

# Step 1: Extract GO IDs of enriched terms
speciosa_terms <- speciosa_go$GO.ID
syriaca_terms  <- syriaca_go$GO.ID

# Step 2: Reverse GO map (GO term → gene ID)
go_to_gene <- map(names(go_map), function(gene) {
  setNames(rep(gene, length(go_map[[gene]])), go_map[[gene]])
}) %>% unlist()

# Step 3: Create named vector GO.ID → gene_id
names(go_to_gene) <- names(go_to_gene)  # GO:XXXXXXX → gene_id

# Step 4: Get gene IDs linked to enriched GO terms
speciosa_genes <- unique(go_to_gene[names(go_to_gene) %in% speciosa_terms])
syriaca_genes  <- unique(go_to_gene[names(go_to_gene) %in% syriaca_terms])

library(VennDiagram)

# Save to file
png("LEAF_Venn_fixed.png", width = 6, height = 6, units = "in", res = 300)

draw.pairwise.venn(
  area1 = length(speciosa_genes),
  area2 = length(syriaca_genes),
  cross.area = length(intersect(speciosa_genes, syriaca_genes)),
  category = c("Speciosa functional genes", "Syriaca functional genes"),
  fill = c("yellow", "orchid"),
  cex = 2,              # Number size
  cat.cex = 1.2,        # Category label size
  cat.pos = c(0, 0),    # Puts labels on top of each circle
  cat.dist = c(0.1, 0.1),  # Controls distance of category labels from circles
  lty = "blank",
  main = "Overlap of Functionally Annotated Genes (Leaf)",
  main.cex = 1
)

dev.off()
