

#Package installation

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

#Note: AnnotationDBI checks for updates each time it is run, consider skipping it's install check if recently updated
if (!require("AnnotationDBI", quietly = TRUE))
  BiocManager::install("AnnotationDBI")
if (!require("apeglm", quietly = TRUE))
  BiocManager::install("apeglm")
if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!require("DOSE", quietly = TRUE))
  BiocManager::install("DOSE")
if (!require("enrichplot", quietly = TRUE))
  BiocManager::install("enrichplot")
if (!require("GenomicFeatures", quietly = TRUE))
  BiocManager::install("GenomicFeatures")


if (!require("pheatmap", quietly = TRUE))
  BiocManager::install("pheatmap")
if (!require("txdbmaker", quietly = TRUE))
  BiocManager::install("txdbmaker")
if (!require("tximport", quietly = TRUE))
  BiocManager::install("tximport")

#yeast annotation database
if (!require("org.Sc.sgd.db", quietly = TRUE))
  BiocManager::install("org.Sc.sgd.db")

#Package loading
library("AnnotationDbi")
library(clusterProfiler)
library(DOSE)

library(enrichplot)
library("dplyr")
library("GenomicFeatures")

library(ggplot2)
library(org.Sc.sgd.db)
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

res2 <- results(dds, name="Stage_Thin_biofilm_vs_Early_biofilm")

#Take a look at the table:
res

# Shrinkage and plots
resLFC <- lfcShrink(dds, coef="Stage_Mature_biofilm_vs_Early_biofilm", type="apeglm")

resLFC2 <- lfcShrink(dds, coef="Stage_Thin_biofilm_vs_Early_biofilm", type="apeglm")

#MA plot without shrinkage
plotMA(res, ylim=c(-2,2))

#With shrinkage
plotMA(resLFC, ylim=c(-2,2))

plotMA(resLFC2, ylim=c(-2,2))

#Volcano plot

#We need to mark genes as upregulated, downregulated, or not significant to colour them in ggplot
#We'll also exclude genes with under 2-fold change (log2FoldChange < 1)
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                             ifelse(res_df$log2FoldChange > 0, "Up", "Down"), "Not Sig")
res_df <- na.omit(res_df)

#As above, but for thin biofilm vs early comparison
res_df2 <- as.data.frame(resLFC2)
res_df2$gene <- rownames(res_df2)
res_df2$significant <- ifelse(res_df2$padj < 0.05 & abs(res_df2$log2FoldChange) > 1, 
                             ifelse(res_df2$log2FoldChange > 0, "Up", "Down"), "Not Sig")
res_df2 <- na.omit(res_df2)



# We plot log2foldchange against -log10(adjusted p value)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
       title = "Volcano Plot: Mature vs Early Biofilm") +
  theme(legend.position = "right")

ggplot(res_df2, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
       title = "Volcano Plot: Thin vs Early Biofilm") +
  theme(legend.position = "right")

#Heatmap
# Remove NA values
resLFC <- na.omit(resLFC)
resLFC2 <- na.omit(resLFC2)

# Select top 20 genes by p values from non-NA genes (could also sort by logfoldchange)
top_genes <- head(order(abs(resLFC$padj), decreasing = FALSE), 20)
gene_names <- rownames(resLFC)[top_genes]

top_genes2 <- head(order(abs(resLFC2$padj), decreasing = FALSE), 20)
gene_names2 <- rownames(resLFC2)[top_genes2]

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds)
# Store counts in a matrix for the heatmap
mat <- assay(vsd)[gene_names, ]
mat2 <- assay(vsd)[gene_names2, ]

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

pheatmap(mat2, 
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



#Functional Annotation Comparison section

# Make a dataframe containing only our results table post-shrinkage, and with NA values trimmed

res_df <- as.data.frame(resLFC)

res_df2 <- as.data.frame(resLFC2)

# Convert ORFs to Entrez IDs
# Remove version numbers (e.g., .9 from ENSG00000189221.9)
ORF_ids <- rownames(res_df)
ORF_ids_clean <- sub("\\..*", "", ORF_ids)

ORF_ids2 <- rownames(res_df2)
ORF_ids_clean2 <- sub("\\..*", "", ORF_ids2)

# Map to Entrez IDs
gene_map <- bitr(ORF_ids_clean, 
                 fromType = "ORF", 
                 toType = c("ENTREZID", "GENENAME"),
                 OrgDb = org.Sc.sgd.db)

gene_map2 <- bitr(ORF_ids_clean2, 
                 fromType = "ORF", 
                 toType = c("ENTREZID", "GENENAME"),
                 OrgDb = org.Sc.sgd.db)

# Add the mapping to results
res_df$ORF <- sub("\\..*", "", rownames(res_df))
res_df <- merge(res_df, gene_map, by = "ORF", all.x = TRUE)

res_df2$ORF <- sub("\\..*", "", rownames(res_df2))
res_df2 <- merge(res_df2, gene_map2, by = "ORF", all.x = TRUE)



# Define significant genes (Here we chose p < 0.05 and foldchange > 2, but other values could be justified)
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

sig_genes2 <- res_df2 %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

# Define our background list of genes to compare to
# Rembember that ORA needs an "interesting gene" set and a background set.
all_genes <- res_df %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

all_genes2 <- res_df2 %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

# Here we'll do a GO analysis for only Biological Process 
# (Try Molecular Function or Cellular Component instead and see what you get!)
ego_bp <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)
head(as.data.frame(ego_bp))

#Performing seperate GO analysis for Thin Biofilm comparison
ego_bp2 <- enrichGO(gene = sig_genes2,
                   universe = all_genes2,
                   OrgDb = org.Sc.sgd.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)

#Need to define significant and background gene list by ORF to function with KEGG

sig_genes_ORF <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

all_genes_ORF <- res_df %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

sig_genes_ORF2 <- res_df2 %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

all_genes_ORF2 <- res_df2 %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

# Now we'll do a KEGG analysis
kegg_enrich <- enrichKEGG(gene = sig_genes_ORF,
                          organism = 'sce',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

head(as.data.frame(kegg_enrich))


kegg_enrich2 <- enrichKEGG(gene = sig_genes_ORF2,
                          organism = 'sce',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

head(as.data.frame(kegg_enrich2))


# Dot plots
dotplot(ego_bp, showCategory = 20, title = "GO Biological Process")

dotplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment")


# Bar plot
barplot(ego_bp, showCategory = 15, title = "GO Biological Process")

barplot(kegg_enrich, showCategory = 15, title = "KEGG Biological Process")

# Enrichment map (Displays linked GO terms)
emapplot(pairwise_termsim(ego_bp), showCategory = 30)

# Notice that none of these plots split our genes by upregulated/downregulated?
# We would need to split those out ourselves as our gene set of interest
# Here's an example for upregulation:

upregulated_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_up <- enrichGO(gene = upregulated_genes,
                      universe = all_genes,
                      OrgDb = org.Sc.sgd.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = FALSE)


upregulated_genes2 <- res_df2 %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_up2 <- enrichGO(gene = upregulated_genes2,
                      universe = all_genes2,
                      OrgDb = org.Sc.sgd.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = FALSE)


#Plots of upregulated genes:

dotplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")

barplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")


#Kegg plot of upregulated genes

upregulated_genes_ORF <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

kegg_up <- enrichKEGG(gene = upregulated_genes_ORF,
                          organism = 'sce',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

upregulated_genes_ORF2 <- res_df2 %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

kegg_up2 <- enrichKEGG(gene = upregulated_genes_ORF2,
                      organism = 'sce',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

barplot(kegg_up, showCategory = 15, title = "KEGG - Upregulated Genes")

dotplot(kegg_up, showCategory = 15, title = "KEGG - Upregulated Genes")


#Downregulate gene plots
downregulated_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_down <- enrichGO(gene = downregulated_genes,
                        universe = all_genes,
                        OrgDb = org.Sc.sgd.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = FALSE)

downregulated_genes2 <- res_df2 %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_down2 <- enrichGO(gene = downregulated_genes2,
                        universe = all_genes2,
                        OrgDb = org.Sc.sgd.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = FALSE)

dotplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes")

barplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes")


# as above, but for KEGG downregulation
downregulated_genes_ORF <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

kegg_down <- enrichKEGG(gene = downregulated_genes_ORF,
                      organism = 'sce',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

downregulated_genes_ORF2 <- res_df2 %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

kegg_down2 <- enrichKEGG(gene = downregulated_genes_ORF2,
                        organism = 'sce',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)

barplot(kegg_down, showCategory = 15, title = "KEGG - Downregulated Genes")

dotplot(kegg_down, showCategory = 15, title = "KEGG - Downregulated Genes")


# All plots for Thin-biofilm functional annotation anlysis:

emapplot(pairwise_termsim(ego_bp2), showCategory = 30)

dotplot(ego_bp_up2, showCategory = 15, title = "GO BP - Thin Biofilm Upregulated Genes")

dotplot(ego_bp_down2, showCategory = 15, title = "GO BP - Thin Biofilm Downregulated Genes")

dotplot(kegg_up, showCategory = 15, title = "KEGG - Thin Biofilm Upregulated Genes")

dotplot(kegg_down2, showCategory = 15, title = "KEGG - Thin Biofilm Downregulated Genes")



