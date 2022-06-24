# LIBRARIES, LOG FILES
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library(purrr)
library(tidyverse)
samples <- read.table(snakemake@input[['samples']], header=TRUE, row.names="sample", check.names=FALSE)
coldata <- samples

# CREATE COUNT MATRIX
f_files <- snakemake@input[['counts']] #list.files("results/superEnhancer/rna_expression/MSTC", pattern = "featureCount.txt", full.names = T)


read_in_feature_counts<- function(file){
        cnt<- read_tsv(file, col_names =T, comment = "#")
        cnt<- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
        return(cnt)
}

raw_counts <- map(f_files, read_in_feature_counts)
cts <- purrr::reduce(raw_counts, inner_join)
cts_df <- as.data.frame(cts)
cts <- cts_df[,-1]
rownames(cts) <- cts_df[,1]
colnames(cts) <- rownames(coldata) #samples$sample

# DESEQ
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
norm_counts = counts(dds, normalized=T)
contrast <- c("condition", "normal", "tumor") #TODO
res <- results(dds, contrast=contrast)
res <- lfcShrink(dds, contrast=contrast, res=res)#, type="ashr")
res <- res[order(res$padj),] # sort by p-value

#PCA plot
# print(nrow(dds)) # 1799 rows  = transcripts
# vsd <- vst(dds, blind=FALSE, nsub=3) # gives : timateDispersionsFit(object.sub, fitType = fitType, quiet = TRUE) :  all gene-wise dispersion estimates are within 2 orders of magnitude
# svg(snakemake@output[["pca"]])
# plotPCA(vsd, intgroup="condition")

ldat <- normTransform(dds)
svg(snakemake@output[["pca"]])
plotPCA(ldat)

#RESULTS
dev.off()
write.table(data.frame("gene"=rownames(res),res), file=snakemake@output[["table"]], row.names=FALSE, sep='\t')

