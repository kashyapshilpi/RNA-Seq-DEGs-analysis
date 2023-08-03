#import library
library(limma)
library(DESeq2)
library(ggplot2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(cancerTiming)
library(geneplotter)
args <- commandArgs(trailingOnly = TRUE)
library(Rsubread)
library(dplyr)
library(genefilter)
library(stringr)
library(LSD)
#To import the file into R and name variable for it
counts = read.table("/home/scbb/abhijit/new_all_gene_count.txt", sep="\t", header=T,row.names=c(1))
#to extract the columns from the main file
x = counts[, c(8,14,15,1,12,13)]
#define condition for it
condition<-c("C","C","C","T","T","T")
#define metadata
coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
filter <- res[(abs(res$log2FoldChange)>2 & res
$padj>0),]
#The resulting table can be written as text file for further analyses
write.table(filter,"regulated.txt", quote=F,sep="\t", col.names = NA)
#for counting up-regulated &down-regulated genes
resSig <- subset(filter, padj >0)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_4vsbtp_5_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_4vsbtp_5_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
#to plot boxplot of Normalized and Un-normalized data
dds = estimateSizeFactors(dds)
pdf("boxplot_normalized.pdf",width=10,height=10)
colors = c(rep("blue",3), rep("green",3))
boxplot(log2(counts(dds, normalized=TRUE)+1), col=colors, outline = FALSE, main="Box-plot of Normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("
), fill=c("green","blue"), cex=0.8)
dev.off()
pdf("boxplot_unnormalized.pdf",width=10,height=10)
boxplot(log2(x+1), col=colors, outline = FALSE, main="Box-plot of Un-normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("btp_1","btp_4"), fill=c("green","blue"), cex=0.8)
dev.off()
#scatter plot of gene Dispersion data
dds <- estimateDispersions(dds)
pdf("Dispersion.pdf",width=10,height=10)
plotDispEsts(dds)
dev.off()
#MA plot of log2Foldchange and Mean normalized count to visualize differences between measurements taken in two samples
pdf("MA.pdf",width=10,height=10)
plotMA(results)
dev.off()
#Volcano plot
pdf("volcano.pdf",width=10,height=10)
with(filter, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot", col="grey", xlim=c(-10,10)))
with(subset(filter,padj<0.01 & log2FoldChange>0),points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(filter,padj<0.01 & log2FoldChange<0),points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(h=2, lty=2)
abline(v=-1, lty=2)
abline(v=1, lty=2)
dev.off()
#Heatmap plotting
rld = rlogTransformation(dds, blind=T)
pdf("sample_heatmap1.pdf",width=10,height=10)
distRL = dist(t(assay(rld)))
mat=as.matrix(distRL)
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
pdf("sample_heatmap2.pdf",width=10,height=10)
pheatmap(mat, clustering_distance_rows=distRL, clustering_distance_cols=distRL, col=hmcol)
dev.off()
#PCA plot to check Batch effect
pdf("pca.pdf")
pca = DESeq2::plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
p<-ggplot(pca,aes(x=PC1,y=PC2,color=group, label=row.names(pca) ))
p<-p+geom_point()
p
dev.off()
#Multi-ecdf
pdf("ecdf.pdf",width=10,height=10)
multiecdf(counts(dds, normalized=TRUE)[,], xlab="Meancounts", xlim=c(0,1000))
dev.off()
