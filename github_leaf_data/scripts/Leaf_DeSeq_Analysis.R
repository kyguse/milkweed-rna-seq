#####ALL Leaf (combined) - SPECIOSA/Syriaca Analysis 04/29/2025 --- Including gene name in annotation file
###Re-working script to include proper filtering

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

# Define the parent directory containing the 7 sample folders
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
samples <- read.csv("samples_file.txt", sep="\t")

###Renaming column data to match samples file
names(filtered_filtered_counts) <- c("sp_rep1", "sp_rep2", "sp_rep3", "sp_rep4", "sp_rep5", "sp_rep6", "sp_rep7", "sp_rep8", "sy_rep1", "sy_rep2", "sy_rep3", "sy_rep4", "sy_rep5", "sy_rep6", "sy_rep7")

##Since filtered_filtered_counts came from the tximport$counts we need to round 
###(since these can have decimals)
filtered_filtered_counts <- round(filtered_filtered_counts)

# Create DESeq2 dat
dds <- DESeqDataSetFromMatrix(countData = filtered_filtered_counts, colData = samples, design = ~ species)

# Run DESeq2
dds <- DESeq(dds)

# Results for syriaca vs speciosa
res <- results(dds, contrast = c("species", "speciosa", "syriaca"))

# View summary
summary(res)

# Order by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Write results
write.csv(as.data.frame(resOrdered), file = "DESeq2_results_LEAF_ALL_07_29_25.csv")

library(stringr)
# Load pfam2go, skipping comments
pfam2go_raw <- readLines("pfam2go copy.txt")
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

###Load annotation file
annot <- read.delim("trinotate_annotation_report.tsv",sep = "\t", header = TRUE)
colnames(annot)

head(annot$X.gene_id)
head(annot$Pfam)
head(annot$sprot_Top_BLASTP_hit)

annot$Pfam_IDs <- str_extract_all(annot$Pfam, "PF\\d+\\.\\d+")
annot$Pfam_IDs <- lapply(annot$Pfam_IDs, function(x) gsub("\\.\\d+$", "", x))

annot$sprot_Top_BLASTP_hit<-str_extract_all(annot$sprot_Top_BLASTP_hit, "sp\\|([^|]+)\\|([^\\^|]+)")

# Check the result
head(annot$Pfam_IDs, 5)
head(annot$sprot_Top_BLASTP_hit, 5)

library(stringr)
###Continue to clean up string in sprot_Top_BLASTP_hit
annot$sprot_Top_BLASTP_hit<-str_replace(annot$sprot_Top_BLASTP_hit, ",.*", "" )

# Remove 'sp|' from the beginning and '"' from the end
annot$sprot_Top_BLASTP_hit <- gsub('^sp\\|', '', annot$sprot_Top_BLASTP_hit)
annot$sprot_Top_BLASTP_hit <- gsub('"$', '', annot$sprot_Top_BLASTP_hit)

# Check again
head(annot$Pfam_IDs, 5)
head(annot$sprot_Top_BLASTP_hit)

library(tidyr)
library(dplyr)

####Cleaning dataframe to use only sprot_Top_BLASTP_hit
annot_clean <- annot %>%
  select(X.gene_id, transcript_id, sprot_Top_BLASTP_hit)

# Convert list-columns to character (flatten character(0) to "")
annot_clean <- annot_clean %>%
  mutate(sprot_Top_BLASTP_hit = sapply(sprot_Top_BLASTP_hit, function(x) if(length(x) == 0) NA else x))

# Now apply filtering again
annot_clean <- annot_clean %>%
  select(X.gene_id, transcript_id, sprot_Top_BLASTP_hit) %>%
  filter(
    !is.na(X.gene_id),
    !is.na(transcript_id),
    !is.na(sprot_Top_BLASTP_hit),
    sprot_Top_BLASTP_hit != "character(0)"
  )

####Still need to clean up sprot_Top_BlASTP_hit column
annot_clean$sprot_Top_BLASTP_hit <- gsub('^c\\("sp\\||\\"|\\)', '', annot_clean$sprot_Top_BLASTP_hit)

###Now you can start to merge

###have to rename column in resOrdered
names(resOrdered)[names(resOrdered) == "gene_id"] <- "X.gene_id"

merged_genes <- merge(resOrdered, annot_clean, by = "X.gene_id", all.x = TRUE)

##Remove NA values 
library(tidyverse)
merged_genes_new <- merged_genes %>% 
  drop_na(sprot_Top_BLASTP_hit)

nrow(merged_genes_new)
write.table(merged_genes_new, "Merged_DE_GENES_LEAF_07_29_25.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

############################################################
######Let's try to some visualization with this merged data
############################################################
library(ggplot2)
library(tibble)
###Raw count distribution
ggplot(filtered_filtered_counts, aes(x = log1p(sp_rep1))) +
  geom_histogram() +
  xlab("Log-transformed raw expression counts - Leaf Data") +
  ylab("Number of genes")


###Make sure the row names in the metadata file are in the same order as the
####column names in the count_data file
samples<-column_to_rownames(samples, var=names(samples)[1])
rownames(samples)
colnames(filtered_filtered_counts)
all(rownames(samples) == colnames(filtered_filtered_counts))

###Count normalization
dds_leaf<-estimateSizeFactors(dds)

norm_leaf_counts<-counts(dds_leaf, normalized=TRUE)
view(norm_leaf_counts)

###Unsupervised clustering: log transform the data
vst_leaf<-vst(dds_leaf, blind=TRUE)

###Extract the vst matrix from the object
vs_mat_leaf<-assay(vst_leaf)
###Compute the pairwise correlation values
vst_cor_leaf<-cor(vs_mat_leaf)
View(vst_cor_leaf)

###Make the heatmap
library(pheatmap)
png("NEW_LEAF_sample-to-sample_distances.png", height = 5, width = 7, units = 'in', res = 600)
pheatmap(vst_cor_leaf, annotation = select(samples, species))
dev.off()

###Make the PCA
library(ggplot2)
png("NEW_leaf_PCA_GENES.png", height = 4, width = 6, units = 'in', res = 600)
p<-plotPCA(vst_leaf, intgroup="species")
p + scale_color_manual(values = c("speciosa"="olivedrab", "syriaca" = "limegreen"))
dev.off()

mean_counts<-apply(filtered_filtered_counts[,1:7],1,mean)
variance_counts<-apply(filtered_filtered_counts[,1:7], 1,var)

mean_var<-data.frame(mean_counts, variance_counts)

png("NEW_mean_variance_plot_LEAF.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(mean_var) + geom_point(aes(x=mean_counts, y=variance_counts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("mean counts per gene") +
  ylab("Variance per gene")
dev.off()

###Plot dispersion estimates
png("New_Dispersion_LEAF.png", height = 4, width = 6, units = 'in', res = 600)
plotDispEsts(dds_leaf)
dev.off()

new_res<-results(dds_leaf, contrast=c("species", "speciosa", "syriaca"), alpha=0.05)
new_res

###Explore with MA plot
png("NEW_MA_PLOT_0.05_LEAF.png", height = 4, width = 6, units = 'in', res = 600)
plotMA(new_res, ylim=c(-8,8))
dev.off()

###Subset for Significant DE genes
leaf_res_sig<-subset(merged_genes_new, padj <0.05)
leaf_res_sig<-leaf_res_sig %>%
  arrange(padj)
View(leaf_res_sig)
View(merged_genes_new)
write.table(leaf_res_sig, "Merged_DE_GENES_SIGNIFICANT_LEAF_07_29_25.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Create heatmap with significant genes
library(pheatmap)
library(dplyr)
norm_leaf_counts<-norm_leaf_counts[leaf_res_sig$X.gene_id,]
library(RColorBrewer)
display.brewer.all()
heat_colors<-brewer.pal(6, "YlGnBu")


##Run heatmap
png("NEW_LEAF_heatmap_Significant_DE_Genes.png", height = 5, width = 5, units = 'in', res = 600)
pheatmap(norm_leaf_counts, color=heat_colors, cluster_rows = T, show_rownames = F, annotation = select(samples, species),scale = "row")
dev.off()

###Run Volcano Plot

###Have to collapse duplicates by keeping the most significant entry

library(dplyr)
library(dplyr)
df_volcano <- data.frame(X.gene_id = paste0("Gene", 1:1000),log2FoldChange = rnorm(1000),padj = runif(1000))

##Take out any NA padj values
df_volcano <- df_volcano %>% filter(!is.na(padj))

###Try to make gene_id_unique
df_volcano$X.gene_id<-make.unique(as.character(df_volcano$X.gene_id))

df_volcano <- merged_genes_new %>%
  mutate(sig = case_when(
    padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant",
    TRUE ~ "Not Significant"
  ))

volcano_plot<-ggplot(df_volcano, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "firebrick2", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

print(volcano_plot)

png("NEW_DE_LEAF_Volcano_genes.png", height = 4, width = 5, units = 'in', res = 600)
print(volcano_plot)
dev.off()

###Plotting by significance (Top15)
df_volcano$gene_name<-rownames(df_volcano)

top15_genes<-df_volcano[order(df_volcano$padj), ] [1:15, ]

top15_genes

write.csv(top15_genes, "top_15_genes_LEAF_07_29_25.csv", row.names=FALSE)

###Take out Accession number for plotting
# Split by '|' and keep the second part
top15_genes$gene_clean <- sapply(strsplit(top15_genes$sprot_Top_BLASTP_hit, "\\|"), `[`, 2)
print(top15_genes$gene_clean)

###Filter out any NAs
df_volcano <- df_volcano %>% filter(!is.na(padj), !is.na(log2FoldChange))

####Okay now need to filter for most significant transcript per gene
top_unique_genes <- top15_genes %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  group_by(gene_clean) %>%
  slice_min(padj, n = 1) %>%  # keep most significant transcript per gene
  ungroup()

###Plot Volcano Plot
library(ggrepel)

png("NEW_DE_LEAF_volcano_top15_labeled_geneNAMES.png", height = 5, width = 7,units = 'in', res = 600)
ggplot(df_volcano, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = case_when(
    is.na(padj) | is.na(log2FoldChange) ~ "Missing",
    padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
    TRUE ~ "Not significant"
  )), alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_unique_genes,
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_clean),  # <-- need x and y here too!
                  size = 2, max.overlaps = 30,
                  box.padding = 0.5, point.padding = 0.3, 
                  segment.color = "lightgrey", min.segment.length = 0) +
  scale_color_manual(values = c(
    "Significant" = "firebrick2",
    "Not significant" = "grey",
    "Missing" = "black"
  )) +
  labs(color = "Significance") +
  theme_minimal()

dev.off()

