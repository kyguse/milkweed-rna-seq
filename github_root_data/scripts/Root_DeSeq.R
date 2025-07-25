#####ALL ROOT (from pdata) - SPECIOSA/Syriaca Analysis 06/30/2025 --- Including gene name in annotation file

# Load libraries
library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(ggplot2)

# Load transcript-to-gene mapping
tx2gene <- read.delim("trinity_out_dir.Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- tx2gene[, c(2,1)]  # Flip columns: transcript, gene
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Define the directory containing the 7 folders with quant.sf
base_dir <- "salmon_output"  # Change this if it's different

# Automatically list sample directories
sample_dirs <- list.dirs("salmon_output", recursive = FALSE)

# Construct paths to quant.sf files inside each sample folder
files <- file.path(sample_dirs, "quant.sf")

# Name the files vector using the folder names (sample names)
names(files) <- basename(sample_dirs)

print(files)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi_counts<-as.data.frame(txi)
colnames(txi_counts)

###Creating dataframe with just the count data
count_df<-txi_counts[, grepl("counts", names(txi_counts))]
colnames(count_df)
count_df$countsFromAbundance<-NULL

###Filtering Count data
filtered_counts <- count_df[rowSums(count_df) > 0, ]

###Let's filter it again to reduce noise  - keep transcipts with at least 10 counts and in at least 2 samples
filtered_filtered_counts <- filtered_counts[rowSums(filtered_counts >= 10) >= 2, ]

##Now let's collapse isoforms for gene level analyses

##Check out transcript ids
head(rownames(filtered_filtered_counts))

###Extract gene ID by removing isoform part (e.g., "_i1")
gene_ids <- sub("_i[0-9]+$", "", rownames(filtered_filtered_counts))

###Collapse count by genes
gene_counts <- rowsum(filtered_filtered_counts, group = gene_ids)

###Check before and after - in this case genes are already collapsed
nrow(filtered_filtered_counts)   # Number of isoforms
nrow(gene_counts)       # Number of genes

###Load metadata
samples <- read.csv("samples_file_root.txt", sep="\t")

###Renaming column data to match samples file
names(filtered_filtered_counts) <- c("sp_rep1", "sp_rep2", "sp_rep3", "sp_rep4", "sp_rep5", "sp_rep6", "sy_rep1", "sy_rep2", "sy_rep3", "sy_rep4", "sy_rep5", "sy_rep6")

##Since filtered_filtered_counts came from the tximport$counts we need to round 9
###(since these can have decimals)
filtered_filtered_counts <- round(filtered_filtered_counts)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_filtered_counts,
                              colData = samples,
                              design = ~ species)
# Run DESeq2
dds <- DESeq(dds)

# Results for syriaca vs speciosa
res <- results(dds, contrast = c("species", "speciosa", "syriaca"))

# View summary
summary(res)

# Order by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Write results
write.csv(as.data.frame(resOrdered), file = "DESeq2_results_ROOT_ALL_FINAL.csv")

library(stringr)
# Load pfam2go
pfam2go_raw <- readLines("pfam2go.txt")
pfam2go_clean <- grep("^Pfam:", pfam2go_raw, value = TRUE)

# Parse Pfam ID and GO ID from each line
pfam_entries <- str_match(pfam2go_clean, "^Pfam:(PF\\d+)\\s+.*>.*GO:.*;\\s*(GO:\\d+)$")

# Remove lines that didn't match
pfam_entries <- pfam_entries[!is.na(pfam_entries[, 2]), ]

# Create mapping table
pfam_map <- data.frame(
  Pfam_ID = pfam_entries[, 2],
  GO_ID   = pfam_entries[, 3],
  stringsAsFactors = FALSE
)

###Check It
head(pfam_map)

resOrdered <- as.data.frame(res[order(res$padj), ])
resOrdered$gene_id <- rownames(resOrdered)  # Store rownames as a column

head(resOrdered)

library(dplyr)
library(tidyverse)

# Load annotation file
anno <- read_tsv("trinotate_annotation_report.tsv", guess_max = 10000)

# Select relevant columns
anno_clean <- anno %>%
  select(
    gene_id = `#gene_id`,
    transcript_id,
    blastp = sprot_Top_BLASTP_hit,
    go_pfam = gene_ontology_Pfam
  )

# Extract Accession and Symbol from blastp field
anno_clean <- anno_clean %>%
  mutate(blastp_short = str_extract(blastp, "sp\\|[^|]+\\|[^|]+")) %>%
  separate(blastp_short, into = c("db", "accession", "symbol"), sep = "\\|", remove = FALSE)

# Parse and shorten GO terms (take first 2)
anno_clean <- anno_clean %>%
  mutate(
    go_terms = str_split(go_pfam, " "),
    go_terms_short = map_chr(go_terms, ~ paste(head(.x, 2), collapse = "; "))
  )

# Anno Final
anno_final <- anno_clean %>%
  select(gene_id, transcript_id, accession, symbol, go_terms_short)

# View preview
print(head(anno_final, 10))

# Save to file
write_tsv(anno_final, "cleaned_trinotate_annotations.tsv")

###Now you can merge
merged_genes <- merge(resOrdered, anno_final, by = "gene_id", all.x = TRUE)

# Filter to retain rows where at least one of the annotation columns is not missing
resFiltered <- merged_genes %>%
  filter(!(is.na(accession) & is.na(symbol) & is.na(go_terms_short)))

resFiltered %>% drop_na()

############################################################
######Let's try to some visualization with this merged data
############################################################
library(ggplot2)
library(tibble)
###Raw count distribution
ggplot(filtered_filtered_counts, aes(x = log1p(sp_rep3))) +
  geom_histogram() +
  xlab("Log-transformed raw expression counts - Root Data") +
  ylab("Number of genes")


###Make sure the row names in the metadata file are in the same order as the
####column names in the count_data file
samples<-column_to_rownames(samples, var=names(samples)[1])
rownames(samples)
colnames(filtered_filtered_counts)
all(rownames(samples) == colnames(filtered_filtered_counts))

###Count normalization
dds_root<-estimateSizeFactors(dds)

norm_root_counts<-counts(dds_root, normalized=TRUE)
view(norm_root_counts)

###Unsupervised clustering: log transform the data
vst_root<-vst(dds_root, blind=TRUE)

###Extract the vst matrix from the object
vs_mat_root<-assay(vst_root)
###Compute the pairwise correlation values
vst_cor_root<-cor(vs_mat_root)
View(vst_cor_root)

###Make the heatmap
library(pheatmap)
png("sample-to-sample distances.png", height = 5, width = 7, units = 'in', res = 600)
pheatmap(vst_cor_root, annotation = select(samples, species))
dev.off()

###Make the PCA
library(ggplot2)
png("root_PCA.png", height = 4, width = 6, units = 'in', res = 600)
p<-plotPCA(vst_root, intgroup="species")
p + scale_color_manual(values = c("speciosa"="tan", "syriaca" = "brown"))
dev.off()

mean_counts<-apply(filtered_filtered_counts[,1:6],1,mean)
variance_counts<-apply(filtered_filtered_counts[,1:6], 1,var)

mean_var <- data.frame(mean_counts, variance_counts)

png("new_mean_variance_plot_ROOT.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(mean_var) + geom_point(aes(x=mean_counts, y=variance_counts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("mean counts per gene") +
  ylab("Variance per gene")
dev.off()

#Plot using variance-stabilizing transformation
vsd <- vst(dds_root, blind = TRUE)
vst_counts <- assay(vsd)

mean_counts <- rowMeans(vst_counts)
variance_counts <- apply(vst_counts, 1, var)
mean_var <- data.frame(mean_counts, variance_counts)

png("new_mean_variance_plot_ROOT_normalized.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(mean_var, aes(x = mean_counts, y = variance_counts)) +
  geom_point(alpha = 0.3) +
  xlab("VST mean per gene") +
  ylab("VST variance per gene")
dev.off()


###Plot dispersion estimates
png("new_Dispersion_ROOT.png", height = 4, width = 6, units = 'in', res = 600)
plotDispEsts(dds_root)
dev.off()

new_res<-results(dds_root, contrast=c("species", "speciosa", "syriaca"), alpha=0.05)
new_res

###Explore with MA plot
png("MA_PLOT_FINAL_ROOT.png", height = 4, width = 6, units = 'in', res = 600)
plotMA(new_res, ylim=c(-8,8))
dev.off()

###Subset for Significant DE genes
root_res_sig<-subset(resFiltered, padj <0.05)
root_res_sig<-root_res_sig %>%
  arrange(padj)
View(root_res_sig)
View(resFiltered)
write.table(root_res_sig, "Merged_DE_GENES_SIGNIFICANT_ROOT.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Create heatmap with significant genes
library(pheatmap)
library(dplyr)
norm_root_counts<-norm_root_counts[root_res_sig$gene_id,]
library(RColorBrewer)
display.brewer.all()
heat_colors<-brewer.pal(6, "PuOr")


##Run heatmap
png("heatmap_root_significant_DE_genes.png", height = 5, width = 5, units = 'in', res = 600)
pheatmap(norm_root_counts, color=heat_colors, cluster_rows = T, show_rownames = F, annotation = select(samples, species),scale = "row")
dev.off()

###Run Volcano Plot

###Have to remove NAs
de_filtered <- root_res_sig %>%
  filter_at(vars(accession, symbol, go_terms_short), any_vars(!is.na(.) & . != ""))

# Remove rows where all three annotation fields are NA or blank
de_filtered_clean <- de_filtered %>%
  filter(
    !(is.na(accession) & is.na(symbol) & (is.na(go_terms_short) | go_terms_short == "."))
  )

###Make the volcano plot
library(ggplot2)

# Define thresholds for coloring
de_filtered_clean <- de_filtered_clean %>%
  mutate(sig = case_when(
    padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant",
    TRUE ~ "Not Significant"
  ))

# Basic volcano plot
png("volcano_root_DE_padj_0.05.png", height = 5, width = 7, units = 'in', res = 600)
ggplot(de_filtered_clean, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "#0072B2", "Not Significant" = "lightgray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Root Transcriptome",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "DEG Status"
  )
dev.off()

####Now trying to label the top 15 most DE genes
library(dplyr)
# Clean and filter for annotated significant genes - and going back to re-generate for only unique annotated gene symbols
# Get top 15 DEGs with unique annotated gene symbols
library(dplyr)

# Ensure it's a data.frame
de_filtered_clean <- as.data.frame(de_filtered_clean)

# Clean and filter for annotated significant genes
top15_genes <- de_filtered_clean %>%
  filter(
    sig == "Significant",
    !is.na(symbol),
    symbol != "NA",
    symbol != "."
  ) %>%
  mutate(symbol_clean = gsub("\\^.*", "", symbol)) %>%  # remove trailing ^sp etc.
  arrange(padj) %>%
  distinct(symbol_clean, .keep_all = TRUE) %>%
  slice(1:15)

# View selected gene symbols
top15_genes$symbol_clean

###Generate the volcano plot - code was generated with some help from Chat GPT
library(ggplot2)
library(ggrepel)

png("volcano_root_DE_padj_0.05_top15.png", height = 5, width = 7, units = 'in', res = 600)
ggplot(de_filtered_clean, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "#0072B2", "Not Significant" = "lightgray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top15_genes,
                  aes(label = symbol_clean),
                  size = 2.5,
                  max.overlaps = 15, color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Root Transcriptome",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "DEG Status"
  )
dev.off()
