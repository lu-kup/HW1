library(DESeq2)

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