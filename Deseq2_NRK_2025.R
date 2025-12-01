library(DESeq2)
library(ggplot2)
library(ggrepel)
library(readr)
library(pheatmap)
library(RColorBrewer)
library (vsn)
library(dplyr)
library(gplots)

#this dataset was analyzed with mRatBN7.2.114 so the map is from  GTF file and is now an excel file named GRCr8.114.csv

GRCr8_114 <- read_csv("GRCr8.114.csv")


#DESEQ2
countData <- read.csv("~counts.csv", header = TRUE, row.names = 1)
countData <-as.matrix(countData)
head(countData)
metaData <- read.csv("~metadata.csv", header = TRUE)
metaData

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~treatment)
dds

#prefilter low count genes
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]

#assign treatment levels.  This will make cisplatin the numerator and saline the denominator
dds$treatment <- factor(dds$treatment, levels = c("HDAC1","control", "HDAC1TGF", "VTGF"))

#DEG
dds <-DESeq(dds)

result <-results(dds, contrast=c("treatment", "HDAC1", "control"), alpha = 0.05)
head(result)
summary(result, )


plotMA(result, main=paste0('Condition: Control vs. Treatment'), ylim=c(-5,5))
rld <- rlogTransformation(dds, blind = TRUE)
PCA <-plotPCA(rld, intgroup = c("treatment"))
PCA + geom_label(aes(label = treatment))


#Plot counts for a single gene. Below is the plot for the gene with the lowest p-value:
plotCounts(dds, gene=which.min(result$padj), intgroup='treatment', pch = 19)

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
vsd <- vst(dds, blind=FALSE)

meanSdPlot(assay(vsd))


#heatmap of samples
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("treatment", "id")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap and sample to sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$id, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#look for sample outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)


## Merge with normalized count data
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"

## Annotate Ensembl ID using Ensembl database  rat
resdata$GeneName <- with(resdata, GRCr8_114$GENENAME[match(EnsemblID, GRCr8_114$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)

resdata$GeneType <- with(resdata, GRCr8_114$GENEBIOTYPE[match(EnsemblID, GRCr8_114$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="~/NRK_HDAC1.csv")



#DEG for TGFbeta

result <-results(dds, contrast=c("treatment", "VTGF", "control"), alpha = 0.05)
head(result)
summary(result, )
## Merge with normalized count data
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"

## Annotate Ensembl ID using Ensembl database v111 rat
resdata$GeneName <- with(resdata, GRCr8_114$GENENAME[match(EnsemblID, GRCr8_114$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)

resdata$GeneType <- with(resdata, GRCr8_114$GENEBIOTYPE[match(EnsemblID, GRCr8_114$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="~/NRK_TGFb.csv")


#DEG for HDAC1 and HDAC1+TGFb
result <-results(dds, contrast=c("treatment", "HDAC1TGF", "VTGF"), alpha = 0.05)
head(result)
summary(result, )
## Merge with normalized count data
resdata <- merge(as.data.frame(result), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "EnsemblID"

## Annotate Ensembl ID using Ensembl database v111 rat
resdata$GeneName <- with(resdata, GRCr8_114$GENENAME[match(EnsemblID, GRCr8_114$GENEID)])
resdata <- resdata %>% relocate(GeneName, .after = EnsemblID)

resdata$GeneType <- with(resdata, GRCr8_114$GENEBIOTYPE[match(EnsemblID, GRCr8_114$GENEID)])
resdata <- resdata %>% relocate(GeneType, .after = GeneName)

write.csv(as.data.frame(resdata), 
          file="~/HDAC1_TGFb.csv")


#########################visualizations###########################
# genes ranked by fold change, and volcano plots were made using Graphpad prism.

#linking DEGs to pathways/GO terms and visualizing results
# Required packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)

setwd("~/_2025/DEG")


# Input: genes stored in a matrix

de_matrix <- read_csv("~_2025/DEG/TGF_GO.csv")   

#filter genes up in HDAC1 treatment
de_filtered <- subset(de_matrix, log2FoldChange > 0.58)
de_genes <- de_filtered$symbol

# Convert gene symbols to Entrez IDs

entrez_ids <- bitr(de_genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

# GO enrichment (Molecular Function)
go_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# KEGG pathway enrichment (optional)
kegg_res <- enrichKEGG(gene = entrez_ids$ENTREZID,
                       organism = "mmu",
                       pvalueCutoff = 0.05)

# Dotplot for GO terms
dotplot(go_mf, showCategory = 5) + ggtitle("GO Molecular Function Enrichment")

# Barplot for KEGG pathways
barplot(kegg_res, showCategory = 5) + ggtitle("KEGG Pathway Enrichment")

# Enrichment map (GO)
emapplot(pairwise_termsim(go_mf))

# Compute pairwise term similarity
term_sim <- pairwise_termsim(ego_selected)

# Enrichment map (GO) w636 h482
emapplot(pairwise_termsim(ego_selected))

# Save results
write.csv(as.data.frame(go_mf), "go_mf_results_UP.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_res), "kegg_results_UP.csv", row.names = FALSE)
#plotted top GO/pathways with Graphpad Prism

#filter genes up in Control treatment
de_filtered <- subset(de_matrix, log2FoldChange < -0.58)
de_genes <- de_filtered$symbol

# Convert gene symbols to Entrez IDs

entrez_ids <- bitr(de_genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

# GO enrichment (Molecular Function)
go_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# KEGG pathway enrichment (optional)
kegg_res <- enrichKEGG(gene = entrez_ids$ENTREZID,
                       organism = "mmu",
                       pvalueCutoff = 0.05)

# Visualization examples
# Dotplot for GO terms
dotplot(go_mf, showCategory = 10) + ggtitle("GO Molecular Function Enrichment")
dotplot(go_mf, showCategory = 5) + ggtitle("GO Molecular Function Enrichment")

# Barplot for KEGG pathways
barplot(kegg_res, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")

# Enrichment map (GO)
emapplot(pairwise_termsim(go_mf))

# Save results
write.csv(as.data.frame(go_mf), "go_mf_results_down.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_res), "kegg_results_down.csv", row.names = FALSE)
#plotted top GO/pathways with Graphpad Prism


#proliferation and cell cycle gene summaries
# Input: genes stored in a matrix

de_matrix <- read_csv("~/_2025/DEG/GO/HDAC1.csv")   
de_matrix <- de_matrix[!is.na(de_matrix$GeneName) & de_matrix$GeneName != "", ]

cellcycle <-read_csv("~/_2025/DEG/GO/proliferation_genes.csv")

#extract gene list
deg_genes <- cellcycle[[1]] 

#Step 2: Subset expression matrix for your genes of interest
subset_expr <- de_matrix[de_matrix$GeneName %in% deg_genes, ]

subset_expr <- subset_expr[order(subset_expr$GeneName), ]
expr_only <- subset_expr[, sapply(subset_expr, is.numeric)]

# Convert to matrix
expr_mat <- as.matrix(expr_only)

# Assign gene names as rownames
rownames(expr_mat) <- subset_expr$GeneName

# 1. Remove rows with any NAs
expr_mat_clean <- expr_mat[apply(expr_mat, 1, function(x) all(!is.na(x))), ]

# 2. Z-score
z_expr <- t(scale(t(expr_mat_clean)))

# 3. Remove rows with NaN/Inf caused by sd=0
z_expr <- z_expr[apply(z_expr, 1, function(x) all(is.finite(x))), ]

# ---- HEATMAP ----
library(pheatmap)
pheatmap(z_expr,
         cluster_rows=FALSE,
         cluster_cols=TRUE,
         show_rownames=TRUE,
         show_colnames=TRUE,
         color=colorRampPalette(c("blue","white","red"))(50),
         main="Heatmap of Cell Cycle Genes (Z-score)")

# Convert to data frame
z_df <- as.data.frame(z_expr)

# Add GeneName column
z_df$GeneName <- rownames(z_expr)

# Put gene names first
z_df <- z_df[, c("GeneName", setdiff(colnames(z_df), "GeneName"))]

# Save to CSV
write.csv(z_df, "zscore_cellcycle_genes.csv", row.names = FALSE)

#proliferation

cellcycle <-read_csv("~/_2025/DEG/GO/proliferation_genes.csv")

#extract gene list
deg_genes <- cellcycle[[2]] 

#Step 2: Subset expression matrix for your genes of interest
subset_expr <- de_matrix[de_matrix$GeneName %in% deg_genes, ]

subset_expr <- subset_expr[order(subset_expr$GeneName), ]
expr_only <- subset_expr[, sapply(subset_expr, is.numeric)]

# Convert to matrix
expr_mat <- as.matrix(expr_only)

# Assign gene names as rownames
rownames(expr_mat) <- subset_expr$GeneName

# 1. Remove rows with any NAs
expr_mat_clean <- expr_mat[apply(expr_mat, 1, function(x) all(!is.na(x))), ]

# 2. Z-score
z_expr <- t(scale(t(expr_mat_clean)))

# 3. Remove rows with NaN/Inf caused by sd=0
z_expr <- z_expr[apply(z_expr, 1, function(x) all(is.finite(x))), ]

# ---- HEATMAP ----
pheatmap(z_expr,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_rownames=TRUE,
         show_colnames=TRUE,
         color=colorRampPalette(c("blue","white","red"))(50),
         main="Heatmap of proliferation genes (Z-score)"). #800 by 1200

# Convert to data frame
z_df <- as.data.frame(z_expr)

# Add GeneName column
z_df$GeneName <- rownames(z_expr)

# Put gene names first
z_df <- z_df[, c("GeneName", setdiff(colnames(z_df), "GeneName"))]

# Save to CSV
write.csv(z_df, "zscore_proliferation
          _genes.csv", row.names = FALSE)
