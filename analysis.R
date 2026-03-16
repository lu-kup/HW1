library(DESeq2)
library(dplyr)
library(tibble)
library(pheatmap)

counts_raw <- read.delim(
    "outputs/feature_counts.out",
    sep = "\t",
    check.names = FALSE,
    comment.char = "#"
)

gene_table <- read.delim(
    "outputs/gene_table.out",
    sep = ",",
    header = FALSE
)

colnames(gene_table) <- c("geneName", "geneType", "ENSG")
countData <-  cbind(ENSG = counts_raw$Geneid,counts_raw[,7:ncol(counts_raw)])

geneData <- merge(countData, gene_table, by = "ENSG")

geneData$ENSG <- gsub("\\..*$", "", geneData$ENSG)
dim(geneData)

codingGeneData <- geneData[geneData$geneType == "protein_coding", ]
dim(codingGeneData)

meta <- read.delim(
    "metadata.csv",
    sep = ",",
    header = FALSE
)
colnames(meta) <- c("Sample", "treatment", "replicate")

bam_cols <- grep("SRR", colnames(codingGeneData))
colnames(codingGeneData)[bam_cols] <- sub(".*/(SRR[0-9]+).*", "\\1", colnames(codingGeneData)[bam_cols])

head(codingGeneData)

# Preparing dataset for DESeq
counts_subset <- codingGeneData[,2:7]
rownames(counts_subset) <- codingGeneData[,1]

dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = meta,
    design = ~ treatment
)

dds$treatment <- relevel(dds$treatment, ref = "shCtrl")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Genes remaining after low-count filter:", nrow(dds), "\n")

# Run DESeq itself
dds <- DESeq(dds)
res <- results(dds,
    contrast = c("treatment", "shESCCAL-1", "shCtrl"),
    alpha = 0.05)

summary(res)

res_shrunk <- lfcShrink(dds,
    coef = "treatment_shESCCAL.1_vs_shCtrl",
    type = "apeglm")

summary(res_shrunk)

res_increased <- res_shrunk[res_shrunk$log2FoldChange > 1, ]
res_decreased <- res_shrunk[res_shrunk$log2FoldChange < -1, ]
res_threshold <- rbind(res_decreased, res_increased)
dim(res_threshold)

write.csv(res_threshold, "outputs/res_threshold.csv")
write.csv(res_shrunk, "outputs/res_full.csv")


dds2 <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = meta,
    design = ~ treatment
)
vst_object <- vst(dds2, blind = TRUE)
vst_res <- assay(vst_object)
head(vst_res)

write.csv(vst_res, "outputs/res_vst.csv")

# Q1
total_reads_per_sample <- colSums(geneData[,2:7])
write.csv(total_reads_per_sample, "outputs/total_reads_per_sample.csv")

total_reads_per_sample_coding <- colSums(codingGeneData[,2:7])
write.csv(total_reads_per_sample, "outputs/total_reads_per_sample_coding.csv")

namedCodingGeneData <- codingGeneData[2:7]
rownames(namedCodingGeneData) <- codingGeneData[,1]
codingGeneDataFiltered <- namedCodingGeneData[rowSums(namedCodingGeneData) > 0, ]
write.csv(dim(codingGeneDataFiltered), "outputs/number_coding_genes.out")

# Q2
cor_mat  <- cor(vst_res, method = "pearson")

anno_df <- meta %>%
  column_to_rownames("Sample") %>%
  dplyr::select(treatment)

anno_colors <- list(
  treatment = c(
    shCtrl = "#66C2A5",
    shESCCAL1 = "#E78AC3"
  )
)

pheatmap(
  cor_mat,
  annotation_col = anno_df,
  annotation_row = anno_df,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("white", "#2166AC"))(50),
  border_color = NA,
  main = "Sample-to-Sample Pearson Correlation",
  fontsize = 10,
  filename = "outputs/sample_sample_heatmap.png"
)

# Q3
# PCA
pca <- prcomp(t(vst_res), scale. = FALSE)
percentVar <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df <- merge(pca_df, meta, by = "Sample")
treatment_colors <- anno_colors$treatment

png("outputs/pca_plot.png", width = 800, height = 600)

plot(
  pca_df$PC1,
  pca_df$PC2,
  col = treatment_colors[pca_df$treatment],
  pch = 19,
  cex = 3,
  xlab = paste0("PC1 (", round(percentVar[1]*100,1), "%)"),
  ylab = paste0("PC2 (", round(percentVar[2]*100,1), "%)"),
  main = "PCA of Samples"
)

grid()

# add labels to each point
text(
  pca_df$PC1,
  pca_df$PC2,
  labels = pca_df$Sample,
  pos = 3,    
  cex = 0.8
)

legend(
  "topright",
  legend = names(treatment_colors),
  col = treatment_colors,
  pch = 19
)

dev.off()

# Q4
png("outputs/gene_count_boxplot.png")

boxplot(
  counts_subset,
  las = 2,              
  col = "lightblue",     
  ylab = "Raw Gene Counts",
  xlab = "Samples",
  main = "Distribution of Gene Counts Across Samples"
)

grid()
dev.off()

png("outputs/gene_count_boxplot_norm.png")

boxplot(
  vst_res,
  las = 2,              
  col = "lightblue",     
  ylab = "Normalized Gene Counts",
  xlab = "Samples",
  main = "Distribution of Gene Counts Across Samples"
)

grid()
dev.off()

# Q5
png("outputs/volcano_plot.png", width = 900, height = 700)

volcano_df <- as.data.frame(res_shrunk)
volcano_df <- volcano_df[!is.na(volcano_df$padj), ]

volcano_df$significance <- "NS"
volcano_df$significance[volcano_df$log2FoldChange > 1 & volcano_df$padj < 0.05] <- "Up"
volcano_df$significance[volcano_df$log2FoldChange < -1 & volcano_df$padj < 0.05] <- "Down"

volcano_colors <- c(
  "Up" = "#E41A1C",
  "Down" = "#377EB8",
  "NS" = "grey70"
)

plot(
  volcano_df$log2FoldChange,
  -log10(volcano_df$padj),
  col = volcano_colors[volcano_df$significance],
  pch = 19,
  cex = 1.2,
  xlab = "log2 Fold Change",
  ylab = "-log10 adjusted p-value",
  main = "Volcano Plot"
)

abline(v = c(-1, 1), lty = 2, col = "black")
abline(h = -log10(0.05), lty = 2, col = "black")

legend(
  "topright",
  legend = c("Upregulated", "Downregulated", "Not significant"),
  col = volcano_colors,
  pch = 19
)

top_genes <- head(volcano_df[order(volcano_df$padj), ], 10)

text(
  top_genes$log2FoldChange,
  -log10(top_genes$padj),
  labels = rownames(top_genes),
  pos = 3,
  cex = 0.7
)

grid()

dev.off()

# Q6
png("outputs/MA_plot.png", width = 900, height = 700)

plotMA(
  res_shrunk,
  ylim = c(-6, 6),
  main = "MA Plot"
)

abline(h = 0, col = "red", lwd = 2)

dev.off()

# Q7
res_df <- as.data.frame(res_shrunk)
res_df <- res_df[order(res_df$padj), ]
top_genes <- rownames(res_df)[1:30]
heatmap_mat <- namedCodingGeneData[top_genes, ]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

anno_df <- meta %>%
  column_to_rownames("Sample") %>%
  dplyr::select(treatment)

png("outputs/top_DEG_heatmap.png", width = 900, height = 800)

pheatmap(
  heatmap_mat_scaled,
  annotation_col = anno_df,
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(50),
  main = paste("Top 30 Differentially Expressed Genes")
)

dev.off()

# Q8
write.csv(summary(res_df[,4:5]), "outputs/pvalue_summary.csv")

png("outputs/pvalue_boxplot.png")

boxplot(
  res_df[,4:5],
  las = 2,              
  col = "lightblue",     
  main = "Distribution of p-value and adj p-value"
)

grid()
dev.off()

# Biological analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)

# Q1
res_df <- res_shrunk %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"  # Not Significant
    ),
    significance = factor(significance, levels = c("Up", "Down", "NS"))
  )

# Count how many in each category
table(res_df$significance)

top_upregulated <- res_df %>% filter(significance == "Up") %>% head(3) %>%
  dplyr::select(gene, log2FoldChange, padj)

top_downregulated <- res_df %>% filter(significance == "Down") %>% head(3) %>%
  dplyr::select(gene, log2FoldChange, padj)

write.csv(top_upregulated, "outputs/top_upregulated.csv")
write.csv(top_downregulated, "outputs/top_downregulated.csv")

# Q2
# Prepare data
sig_up   <- res_df %>% filter(significance == "Up")   %>% pull(gene)
sig_down <- res_df %>% filter(significance == "Down")  %>% pull(gene)
sig_all  <- res_df %>% filter(significance != "NS")    %>% pull(gene)

cat("Upregulated:", length(sig_up), "genes\n")
cat("Downregulated:", length(sig_down), "genes\n")
cat("Total significant:", length(sig_all), "genes\n")

ranked_list <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange)

gene_ranks <- setNames(ranked_list$log2FoldChange, ranked_list$gene)
cat("\nRanked list length:", length(gene_ranks), "\n")

cat("Range of log2FC:", round(range(gene_ranks), 2), "\n")

# Map gene symbols to Entrez IDs (required by KEGG)
sig_entrez <- bitr(
  sig_all,
  fromType = "ENSEMBL",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
)

cat("Mapped", nrow(sig_entrez), "of", length(sig_all), "genes to Entrez IDs\n")

# ORA GO
ego <- enrichGO(
  gene = sig_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",          # BP, MF, CC all included
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE       # convert back to gene symbols for visualization
)

png("outputs/go_barplot.png", width = 1000, height = 800)
barplot(ego, showCategory = 20, title = "Top GO terms")
dev.off()

png("outputs/go_dotplot.png", width = 1000, height = 800)
dotplot(ego, showCategory = 20, title = "GO Enrichment Dotplot")
dev.off()

# Q3
ranked_df <- data.frame(gene = names(gene_ranks), log2FC = gene_ranks) %>%
  arrange(desc(log2FC))

gene_entrez_map <- bitr(
  ranked_df$gene,
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

ranked_df <- merge(ranked_df, gene_entrez_map, by.x = "gene", by.y = "ENSEMBL")  %>%
  arrange(desc(log2FC))
gene_ranks_entrez <- setNames(ranked_df$log2FC, ranked_df$ENTREZID)

# GSEA GO
gsea_go <- gseGO(
  geneList = gene_ranks_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",             # BP, MF, CC
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

png("outputs/gsea_go_top.png", width = 800, height = 600)
gseaplot2(
  gsea_go,
  geneSetID = gsea_go@result$ID[1],
  title = gsea_go@result$Description[1]
)
dev.off()

png("outputs/gsea_go_ridge.png", width = 900, height = 700)
ridgeplot(gsea_go, showCategory = 20)
dev.off()

# GSEA MSigDB
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
msig_list <- split(msig_h$entrez_gene, msig_h$gs_name)

msig_df <- msig_h %>%
  dplyr::select(gs_name, entrez_gene) %>%  # columns: TERM, GENE
  dplyr::rename(TERM = gs_name, GENE = entrez_gene)

head(msig_df)

gsea_msig <- GSEA(
  gene_ranks_entrez,
  TERM2GENE = msig_df,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

png("outputs/gsea_msig_dotplot.png", width = 900, height = 700)
dotplot(gsea_msig, showCategory = 20, title = "Top GSEA MSigDB Hallmark Terms")
dev.off()
