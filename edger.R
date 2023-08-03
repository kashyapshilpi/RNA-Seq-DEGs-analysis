#import library
library("limma")
library("edgeR")
#import the data
count_matrix <- as.matrix(read.table("/media/scbb/data3/shilpi_trainee/shilpi_RNA/Final_count.txt",header = TRUE, sep = "\t", row.names = "Gene"))
#check if it's import properly
head(count_matrix, 2)
#get the desire column from the dataset           #column of contol and test data only 
x = count_matrix[, c(7,8,9,10,11,12)]
#set the condition
sample_info <- c("ctr", "ctr", "ctr", "trt", "trt", "trt")
#create DGEList data class for count and sample information           
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes                      #CPM and TMM 
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
gene <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = gene, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]
#summary
summary(decideTests(object = gene, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="condition_vascular_vs_non_vascular_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="vascular_vs_non_vascular_up_regulated.csv")
write.csv(as.data.frame(down_regulated), file="vascular_vs_non_vascular_down_regulated.csv")


#scatter plot of gene Dispersion data
dge <- estimateDisp(y = dge)
pdf("vascular_vs_non_vascular_Dispersion.pdf",width=10,height=10)
plotBCV(dge)
dev.off()

#Volcano plot
pdf("vascular_vs_non_vascular_volcano.pdf",width=10,height=10)
with(DEG, plot(logFC, -log10(FDR), main="Volcano Plot", col="grey"))
with(subset(DEG,FDR <0.01 & logFC>0),points(logFC, -log10(FDR),pch = 20, col="green"))
with(subset(DEG,FDR<0.01 &logFC<0),points(logFC, -log10(FDR),pch = 20, col="red"))
abline(h=2, lty=2)
abline(v=-1, lty=2)
abline(v=1, lty=2)
dev.off()


