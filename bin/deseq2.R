#!/usr/bin/env Rscript
# Differential expression with DESeq2 (called by the Nextflow pipeline, also usable on its own).
#
# Usage:
#   deseq2.R count_matrix.tsv samplesheet.csv out_dir [reference] [target]
#
# count_matrix.tsv : gene_id + one column per sample (from merge_counts.py)
# samplesheet.csv  : needs columns `sample` and `condition`
# reference/target : condition levels to compare (default: control, treatment)
#
# Environment variables: PADJ (default 0.05), LFC (default 1)

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: deseq2.R count_matrix.tsv samplesheet.csv out_dir [reference] [target]")
out_dir   <- args[3]
ref_level <- if (length(args) >= 4) args[4] else "control"
tgt_level <- if (length(args) >= 5) args[5] else "treatment"
padj_cut  <- as.numeric(Sys.getenv("PADJ", "0.05"))
lfc_cut   <- as.numeric(Sys.getenv("LFC", "1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(out_dir, x)
write_tsv <- function(x, f) write.table(x, out(f), sep = "\t", quote = FALSE, row.names = FALSE)

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
ss     <- read.csv(args[2], stringsAsFactors = FALSE, strip.white = TRUE)
ss     <- ss[ss$condition %in% c(ref_level, tgt_level), c("sample", "condition")]
if (!all(c(ref_level, tgt_level) %in% ss$condition)) {
    stop("Samplesheet must contain both conditions: ", ref_level, " and ", tgt_level)
}
missing <- setdiff(ss$sample, colnames(counts))
if (length(missing)) stop("Samples missing from count matrix: ", paste(missing, collapse = ", "))

counts <- counts[, ss$sample, drop = FALSE]
ss$condition <- relevel(factor(ss$condition), ref = ref_level)
rownames(ss) <- ss$sample

dds <- DESeqDataSetFromMatrix(counts, colData = ss, design = ~ condition)
# genes with almost no reads can't be tested reliably and only add to multiple-testing burden
keep <- rowSums(counts(dds) >= 10) >= min(table(ss$condition))
dds <- dds[keep, ]
message(sum(keep), " genes pass the low-count filter")
dds <- DESeq(dds)

coef <- resultsNames(dds)[2]                                       # condition_<target>_vs_<reference>
res  <- results(dds, name = coef, alpha = padj_cut)
# shrink fold changes of noisy low-count genes towards 0, so rankings aren't driven by noise
res  <- lfcShrink(dds, coef = coef, res = res, type = "apeglm")
summary(res)

tab <- data.frame(gene_id = rownames(res), as.data.frame(res))
tab <- tab[order(tab$padj), ]
sig <- subset(tab, !is.na(padj) & padj < padj_cut & abs(log2FoldChange) > lfc_cut)

write_tsv(tab, "all_genes.tsv")
write_tsv(sig, "DEGs.tsv")
write_tsv(subset(sig, log2FoldChange > 0), "DEGs_up.tsv")
write_tsv(subset(sig, log2FoldChange < 0), "DEGs_down.tsv")
write_tsv(data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE),
          "normalised_counts.tsv")

# ---- Plots -------------------------------------------------------------------
# variance-stabilised counts are used only for plots (PCA, distances), not for testing.
# vst() needs >= 1000 genes; small panels and test data fall back to the full transform.
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)

pdf(out("deseq2_plots.pdf"), width = 8, height = 7)
plotDispEsts(dds)
plotMA(res, ylim = c(-5, 5), main = paste("MA plot:", tgt_level, "vs", ref_level))
print(plotPCA(vsd, intgroup = "condition") + theme_bw() + ggtitle("PCA (variance-stabilised counts)"))

d <- dist(t(assay(vsd)))
pheatmap(as.matrix(d), clustering_distance_rows = d, clustering_distance_cols = d,
         main = "Sample-to-sample distances")

plot_tab <- subset(tab, !is.na(padj))
plot_tab$status <- with(plot_tab, ifelse(padj < padj_cut & log2FoldChange >  lfc_cut, "Up",
                                  ifelse(padj < padj_cut & log2FoldChange < -lfc_cut, "Down", "NS")))
print(ggplot(plot_tab, aes(log2FoldChange, -log10(padj), colour = status)) +
          geom_point(size = 0.8, alpha = 0.7) +
          scale_colour_manual(values = c(Down = "#2b6cb0", NS = "grey70", Up = "#c53030")) +
          geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
          geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
          theme_bw() + ggtitle(paste("Volcano plot:", tgt_level, "vs", ref_level)))
invisible(dev.off())

message(nrow(sig), " DEGs (padj < ", padj_cut, ", |log2FC| > ", lfc_cut, ") written to ", out_dir)
