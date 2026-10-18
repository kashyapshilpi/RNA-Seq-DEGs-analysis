#!/usr/bin/env Rscript
# Differential expression with edgeR (quasi-likelihood F-test), as a cross-check on DESeq2.
#
# Usage:
#   Rscript 03_edger.R count_matrix.tsv samplesheet.csv results/edger
#
# Same inputs as 02_deseq2.R. Optional environment variables: FDR (default 0.05), LFC (default 1)

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 03_edger.R count_matrix.tsv samplesheet.csv out_dir")
out_dir <- args[3]
fdr_cut <- as.numeric(Sys.getenv("FDR", "0.05"))
lfc_cut <- as.numeric(Sys.getenv("LFC", "1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(out_dir, x)

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
ss     <- read.csv(args[2], stringsAsFactors = FALSE)
counts <- counts[, ss$sample_id]
group  <- relevel(factor(ss$condition), ref = "control")

y <- DGEList(counts, group = group)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                       # TMM normalisation

design <- model.matrix(~ group)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)              # treatment vs control

tab <- topTags(qlf, n = Inf)$table
tab <- data.frame(gene_id = rownames(tab), tab)
sig <- subset(tab, FDR < fdr_cut & abs(logFC) > lfc_cut)

write.table(tab, out("all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig, out("DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(sig, logFC > 0), out("DEGs_up.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(sig, logFC < 0), out("DEGs_down.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(out("edger_plots.pdf"), width = 8, height = 7)
plotMDS(y, col = as.integer(group), main = "MDS plot")
plotBCV(y, main = "Biological coefficient of variation")
plotQLDisp(fit)
col <- ifelse(tab$FDR < fdr_cut & tab$logFC >  lfc_cut, "#c53030",
       ifelse(tab$FDR < fdr_cut & tab$logFC < -lfc_cut, "#2b6cb0", "grey70"))
plot(tab$logFC, -log10(tab$FDR), pch = 20, cex = 0.5, col = col,
     xlab = "log2 fold change", ylab = "-log10 FDR", main = "Volcano plot")
abline(v = c(-lfc_cut, lfc_cut), h = -log10(fdr_cut), lty = 2)
dev.off()

message(nrow(sig), " DEGs (FDR < ", fdr_cut, ", |log2FC| > ", lfc_cut, ") written to ", out_dir)
