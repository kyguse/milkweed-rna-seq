##########Top 10 upregulated and downregulated genes for volcano plot annotation

library(dplyr)
library(ggplot2)
library(ggrepel)

# Pick most extreme DEGs
top_n <- 10
top_up <- de_filtered_clean %>%
  filter(sig == "Significant") %>%
  arrange(desc(log2FoldChange)) %>%
  head(top_n)

top_down <- de_filtered_clean %>%
  filter(sig == "Significant") %>%
  arrange(log2FoldChange) %>%
  head(top_n)

top_extreme <- bind_rows(top_up, top_down)

# Shorten label to just symbol (add accession in hover if using plotly)
de_filtered_clean <- de_filtered_clean %>%
  mutate(label = ifelse(gene_id %in% top_extreme$gene_id, symbol, NA))

# Volcano plot with cleaned annotations
png("ROOT_VOLCANO_Significant_DE_Genes_NEW.png", height = 5, width = 7, units = 'in', res = 600)
ggplot(de_filtered_clean, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("gray70", "dodgerblue")) +
  geom_text_repel(aes(label = label), 
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  force = 3,
                  segment.color = "gray40") +
  theme_minimal() +
  labs(title = "Volcano Plot: Annotated Top DEGs",
       x = "Log2 Fold Change", y = "-log10 Adjusted p-value",
       color = "Significance")
dev.off()

###Clean up a bit - let's collapse by Symbol
###First, write.csv the de_filtered_clean file

write.csv(de_filtered_clean,"de_filtered_clean.csv")
