# installing bioconductor and DESeq2
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("DESeq2")
boicLite("genefilter")
install.packages("pheatmap")
install.packages("ggplot2")

library("DESeq2")
library("pheatmap")
library("ggplot2")

# Change this to the directory on your computer where you put the counts directory
setwd("~/EnGen2018/rnaseq_analysis")

#### LOAD IN SAMPLE TABLE FOR DESeq2
# read in sample table from a text file
# sampleTable.txt must be located in R current working directory
# excel can Save As > format > Tab-delimited .txt to get excel data to text file
sampleTable = read.table("SampleTable.txt", sep='\t', header=TRUE, colClasses=c(rep("character",2),rep("factor",8),rep("numeric",4)))
                         
# this line ensures that our Control exposure is the base level in DESeq
sampleTable$Exposure = relevel(sampleTable$Exposure, ref="Control")
                         
#### USE DESeqDataSetFromHTSeqCount FUNCTION OF DESeq TO IMPORT DATA ####
directory = "counts"

# create DESeqDataSet
# This is where we design our linear model formula
raw_data = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~Exposure + Clade + Clade:Exposure)
#raw_data = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~Exposure)
raw_data

# perform differential expression analysis linear model
data = DESeq(raw_data)

# list different results we can look at based on our model
resultsNames(data)

### RESULTS
# sig. differential genes due to Cadmium exposure
# because Clade is in the model, these are only significant if consistent DE across both clades
# running the model with just ~Exposure
results_exposure = results(data, contrast=c("Exposure", "Cadmium", "Control"), alpha=0.05, cooksCutoff = TRUE)
summary(results_exposure)

# Expression differences between the clades (genetic differences)
results_clade = results(data, contrast=c("Clade", "adapted", "nonadapted"), alpha=0.05)
summary(results_clade)

# Genes that responded differently in response to exposure between the two clades
results_interaction = results(data, contrast=list( c("Exposure_Cadmium_vs_Control","ExposureCadmium.Cladenonadapted")), alpha=0.05)
summary(results_interaction)

# order results by adjusted p-value
resultsFDRordered = as.data.frame(results_exposure[order(results_exposure$padj),])
resultsFDRordered[1:20,]

# order by Fold Change - Rank List
rankList = as.data.frame(results_exposure[order(abs(results_exposure$log2FoldChange), decreasing=T),])
rankList[1:10,]

# write results_exposure to file
write.table(resultsFDRordered, file="rnaseq_DE_results_exposure.txt", sep='\t', quote=F)

# add normalized counts to results_exposure
results_exposure_counts <- merge(resultsFDRordered, as.data.frame(counts(data, normalized=TRUE)), by="row.names", sort=FALSE)
results_exposure_counts[20,]
boxplot(as.numeric(results_exposure_counts[20,8:67]) ~ sampleTable$Exposure, main="Metallothionein Expression")

# FDR threshold significant genes
sigGenesFDR = which(results_exposure$padj < 0.05)
length(sigGenesFDR)
write.table(as.data.frame(results_exposure[sigGenesFDR,]), file="exposure_sig_DE_genes.txt", sep='\t', quote=F)

# Fold Change + pvalue threshold significant genes
sigGenesFCpval = which(abs(results_exposure$log2FoldChange) > 1 & results_exposure$pvalue < 0.01)
length(sigGenesFCpval)

# MA plot
options(scipen=1)
plotMA(results_exposure, main="Cd vs Cntl MA Plot", ylim=c(-5.1, 5.1), colSig='green2', cex=0.6)
plot(results_exposure$log2FoldChange ~ results_exposure$baseMean,log='x', ylab="log2 Fold Change", xlab="Mean Expression", pch=16, cex=0.6, col="gray32", main="Cd vs Cntl MA Plot")
abline(h=0, col='red')
points(results_exposure$log2FoldChange[sigGenesFCpval] ~ results_exposure$baseMean[sigGenesFCpval], col = "green2", pch = 16, cex=0.7)

high_exp = subset(as.data.frame(results_exposure),baseMean > 10000)
high_exp = high_exp[order(high_exp$baseMean, decreasing=T),]
high_exp

# Volcano plot
volcanoplot <- function (results_exposure, LFCthresh=1, pvalthresh=0.01, FDRthresh=NULL, main="Volcano Plot", legendpos="bottomright", textcx=1, ...) {
  with(results_exposure, plot(log2FoldChange, -log10(pvalue), pch=16, cex=0.6, main=main, ...))
  if(is.null(FDRthresh)){
    with(subset(results_exposure, abs(log2FoldChange) > LFCthresh & pvalue < pvalthresh), points(log2FoldChange, -log10(pvalue), pch=16, cex=0.7, col="green2", ...))
    abline(h=-log10(pvalthresh), lty=2)
    abline(v=LFCthresh, lty=2)
    abline(v=-LFCthresh, lty=2)
  } else {
    with(subset(results_exposure, padj < FDRthresh ), points(log2FoldChange, -log10(pvalue), pch=16, cex=0.7, col="green2", ...))
  }
}
volcanoplot(results_exposure, FDRthresh=0.1, main="Cd vs Cntl Volcano Plot")
volcanoplot(results_exposure, LFCthresh=1, pvalthres=0.01, main="Cd vs Cntl Volcano Plot")

# PCA of 500 most variable genes among all samples
data_transform = varianceStabilizingTransformation(data)
plotPCA(data_transform, intgroup=c("Clade"))
plotPCA(data_transform, intgroup=c("Clade","Exposure"))

# Heatmap of Significantly DE Genes counts
annot <- as.data.frame(colData(data_transform)[,c("Exposure","Clade")])
sig_genes_transform = assay(data_transform)[sigGenesFDR,]
pheatmap(sig_genes_transform, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=annot, fontsize_col=8)


sampleDists <- dist(t(assay(data_transform)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(data_transform$Clade, data_transform$Isolate, data_transform$Exposure, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

