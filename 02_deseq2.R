#!/usr/bin/env Rscript
# Differential expression with DESeq2.
#
# Usage:
#   Rscript 02_deseq2.R count_matrix.tsv samplesheet.csv results/deseq2
#
# count_matrix.tsv : output of 01_rnaseq_counts.sh (gene_id + one column per sample)
# samplesheet.csv  : columns sample_id,condition (condition = control | treatment)
#
# Optional environment variables: PADJ (default 0.05), LFC (default 1)

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 02_deseq2.R count_matrix.tsv samplesheet.csv out_dir")
out_dir <- args[3]
padj_cut <- as.numeric(Sys.getenv("PADJ", "0.05"))
lfc_cut  <- as.numeric(Sys.getenv("LFC", "1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(out_dir, x)

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
ss     <- read.csv(args[2], stringsAsFactors = FALSE)
counts <- counts[, ss$sample_id]                       # same sample order as the sheet
ss$condition <- relevel(factor(ss$condition), ref = "control")
rownames(ss) <- ss$sample_id

dds <- DESeqDataSetFromMatrix(counts, colData = ss, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= min(table(ss$condition)), ]   # drop low-count genes
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treatment", "control"), alpha = padj_cut)
res <- lfcShrink(dds, coef = "condition_treatment_vs_control", res = res, type = "apeglm")
summary(res)

tab <- data.frame(gene_id = rownames(res), as.data.frame(res))
tab <- tab[order(tab$padj), ]
sig <- subset(tab, !is.na(padj) & padj < padj_cut & abs(log2FoldChange) > lfc_cut)

write.table(tab, out("all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig, out("DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(sig, log2FoldChange > 0), out("DEGs_up.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(sig, log2FoldChange < 0), out("DEGs_down.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Plots -------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)

pdf(out("deseq2_plots.pdf"), width = 8, height = 7)
plotDispEsts(dds)
plotMA(res, ylim = c(-5, 5), main = "MA plot (shrunken log2 fold change)")

print(plotPCA(vsd, intgroup = "condition") + theme_bw() + ggtitle("PCA (VST counts)"))

d <- as.matrix(dist(t(assay(vsd))))
pheatmap(d, clustering_distance_rows = dist(t(assay(vsd))),
         clustering_distance_cols = dist(t(assay(vsd))), main = "Sample distances")

tab$status <- ifelse(!is.na(tab$padj) & tab$padj < padj_cut & tab$log2FoldChange >  lfc_cut, "Up",
              ifelse(!is.na(tab$padj) & tab$padj < padj_cut & tab$log2FoldChange < -lfc_cut, "Down", "NS"))
print(ggplot(subset(tab, !is.na(padj)), aes(log2FoldChange, -log10(padj), colour = status)) +
          geom_point(size = 0.8, alpha = 0.7) +
          scale_colour_manual(values = c(Down = "#2b6cb0", NS = "grey70", Up = "#c53030")) +
          geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
          geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
          theme_bw() + ggtitle("Volcano plot"))
dev.off()

message(nrow(sig), " DEGs (padj < ", padj_cut, ", |log2FC| > ", lfc_cut, ") written to ", out_dir)
