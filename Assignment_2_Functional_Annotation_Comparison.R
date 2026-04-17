
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("ggplot2")
install.packages("tidyverse")


BiocManager::install("AnnotationDBI")
BiocManager::install("apeglm")
BiocManager::install("clusterProfiler")

BiocManager::install("DOSE")
BiocManager::install("enrichplot")
BiocManager::install("GenomicFeatures")

BiocManager::install("txdbmaker")
BiocManager::install("org.Hs.eg.db")



library("AnnotationDbi")
library(clusterProfiler)
library(DOSE)

library(enrichplot)
library("dplyr")
library("GenomicFeatures")

library(ggplot2)
library(org.Hs.eg.db)
library("pheatmap")

library(tidyverse)
library("txdbmaker")
library("tximport")







gtfFile <- file.path("GCF_000146045.2_R64_genomic.gtf")



#Make txDb and gene table for use with tximport. Note: need to explicitly call AnnotationDbi::select to avoid conflicts with dplyr.
txdb <- makeTxDbFromGFF(gtfFile, format = "gtf")

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")


#Load file names from sample table template, then run tximport.


sampleTable <- read.csv(file.path("s_cerevisae_biofilm_accession_lookup.csv"),row.names=1)


files <- file.path("quant_output",sampleTable$SRA,"quant.sf")

txi <- tximport(files, type="salmon", tx2gene=tx2gene)


# DESeq2
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~Stage)
dds <- DESeq(dds)

#Check comparison names
resultsNames(dds)

#Generate a results table for biofilm_stage
res <- results(dds, name="Stage_Mature_biofilm_vs_Early_biofilm")

#Take a look at the table:
res

# Shrinkage and plots
resLFC <- lfcShrink(dds, coef="Stage_Mature_biofilm_vs_Early_biofilm", type="apeglm")

#MA plot without shrinkage
plotMA(res, ylim=c(-2,2))

#With shrinkage
plotMA(resLFC, ylim=c(-2,2))

#Volcano plot

#We need to mark genes as upregulated, downregulated, or not significant to colour them in ggplot
#We'll also exclude genes with under 2-fold change (log2FoldChange < 1)
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                             ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "Not Sig")
res_df <- na.omit(res_df)


# We plot log2foldchange against -log10(adjusted p value)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
       title = "Volcano Plot: Mature vs Early Biofilm") +
  theme(legend.position = "right")

#Heatmap
# Remove NA values
resLFC <- na.omit(resLFC)

# Select top 20 genes by p values from non-NA genes (could also sort by logfoldchange)
top_genes <- head(order(abs(resLFC$padj), decreasing = FALSE), 20)
gene_names <- rownames(resLFC)[top_genes]

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)
# Store counts in a matrix for the heatmap
mat <- assay(vsd)[gene_names, ]

# Add annotation for cell lines and treatment
annotation_df <- sampleTable[, c("Sample_ID", "Stage")]
colnames(annotation_df) <- c("Sample_ID", "Biofilm_Stage")

# Create heatmap
pheatmap(mat, 
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         show_colnames = FALSE,
)

#PCA Plot

# Same VST used above for the heatmap
vsd <- vst(dds)

# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = c("Stage", "Sample_ID"), returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# GGplot code to display samples by colour, and growth stage by shape
ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample_ID, shape = Stage)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples")
coord_fixed()



# Make a dataframe containing only our results table post-shrinkage, and with NA values trimmed

res_df <- as.data.frame(resLFC)


# Convert Ensembl IDs to Entrez IDs
# Remove version numbers (e.g., .9 from ENSG00000189221.9)
ensembl_ids <- rownames(res_df)
ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)

# Map to Entrez IDs
gene_map <- bitr(ensembl_ids_clean, 
                 fromType = "ENSEMBL", 
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)

# Add the mapping to results
res_df$ENSEMBL <- sub("\\..*", "", rownames(res_df))
res_df <- merge(res_df, gene_map, by = "ENSEMBL", all.x = TRUE)