#!/usr/bin/env Rscript
# Differential expression with edgeR (quasi-likelihood F-test), as a cross-check on DESeq2.
#
# Usage:
#   edger.R count_matrix.tsv samplesheet.csv out_dir [reference] [target]
#
# Same inputs as deseq2.R. Environment variables: FDR (default 0.05), LFC (default 1)

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: edger.R count_matrix.tsv samplesheet.csv out_dir [reference] [target]")
out_dir   <- args[3]
ref_level <- if (length(args) >= 4) args[4] else "control"
tgt_level <- if (length(args) >= 5) args[5] else "treatment"
fdr_cut   <- as.numeric(Sys.getenv("FDR", "0.05"))
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
counts <- counts[, ss$sample, drop = FALSE]
group  <- relevel(factor(ss$condition), ref = ref_level)

y <- DGEList(counts, group = group)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                       # TMM normalisation
message(nrow(y), " genes pass filterByExpr")

design <- model.matrix(~ group)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)              # target vs reference

tab <- topTags(qlf, n = Inf)$table
tab <- data.frame(gene_id = rownames(tab), tab)
sig <- subset(tab, FDR < fdr_cut & abs(logFC) > lfc_cut)

write_tsv(tab, "all_genes.tsv")
write_tsv(sig, "DEGs.tsv")
write_tsv(subset(sig, logFC > 0), "DEGs_up.tsv")
write_tsv(subset(sig, logFC < 0), "DEGs_down.tsv")

pdf(out("edger_plots.pdf"), width = 8, height = 7)
plotMDS(y, col = as.integer(group), main = "MDS plot")
plotBCV(y, main = "Biological coefficient of variation")
plotQLDisp(fit)
col <- ifelse(tab$FDR < fdr_cut & tab$logFC >  lfc_cut, "#c53030",
       ifelse(tab$FDR < fdr_cut & tab$logFC < -lfc_cut, "#2b6cb0", "grey70"))
plot(tab$logFC, -log10(tab$FDR), pch = 20, cex = 0.5, col = col,
     xlab = "log2 fold change", ylab = "-log10 FDR",
     main = paste("Volcano plot:", tgt_level, "vs", ref_level))
abline(v = c(-lfc_cut, lfc_cut), h = -log10(fdr_cut), lty = 2)
invisible(dev.off())

message(nrow(sig), " DEGs (FDR < ", fdr_cut, ", |log2FC| > ", lfc_cut, ") written to ", out_dir)
