library(DESeq2)

#  count 
setwd("C:\\Users\\msc2\\Desktop\\Do_not_open\\TP\\DeNovo_Transcriptomics")
count<- as.matrix(read.csv("count_data.csv", sep = "," , row.names = "Transcripts"))

head(count)
dim(count)
meta<- as.matrix(read.csv("metadata.csv" , sep="," , row.names = "samples"))
head(meta)
dim(meta)

all(rownames(meta) == colnames(count))
dds <- DESeqDataSetFromMatrix(countData = count, colData =meta, design = ~condition)
dds
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]
keep

#Run DESeq function and Perform normalization using vst function:
dds <- DESeq(dds)
class(dds)
sum( rowMeans( counts(dds, normalized=TRUE)) > 5 ) #removing the features less than 5
vsd <- vst(dds, blind= TRUE,nsub = 262)
class(vsd)
head(vsd, 3)


# Result of Differential expression analysis
res <- results(dds, contrast=c("condition","treated","untreated"), alpha=0.05)
summary(res)
write.csv(res,"DEG.csv")
f<-read.csv("DEG.csv")



#Heatmap generation:
library("pheatmap")
library("RColorBrewer")

# Select top 20 highly expressed genes (normalized counts)
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:20]

# Use only 'condition' for column annotation (optional)
df <- as.data.frame(colData(dds))[, "condition", drop = FALSE]

# Create the heatmap
pheatmap(assay(vsd)[select, ], 
         scale = "row",                   # z-score normalization across samples
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,            # Set to FALSE if still too crowded
         show_colnames = TRUE,
         annotation_col = df,
         fontsize_row = 6,                # Smaller font for gene names
         angle_col = 45,                  # Rotate sample labels for better readability
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(255))  # nice diverging color palette




# Load libraries
library(ggplot2)
library(ggrepel)

# Convert results to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove NA
res_df <- na.omit(res_df)

# Define regulation categories
res_df$regulation <- "Not Significant"

res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")





# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Convert DESeq2 results to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove NA values
res_df <- na.omit(res_df)

# Define significance
res_df$significant <- "No"
res_df$significant[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Yes"

# Create volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

# Optional: label top significant genes
top_genes <- res_df[order(res_df$padj), ][1:10, ]

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot with Labels",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

