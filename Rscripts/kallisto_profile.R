library(DESeq2)
library(tximport)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)

samples <- read.csv("units.tsv",sep="\t",header=TRUE)
samples <- subset(samples, select = -c(fq1, fq2) )
exprnames <- do.call(paste,c(samples[c("sample","unit")],sep="."))
files <- file.path("results",exprnames,"abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- exprnames
colnames(txi.kallisto$abundance) <- exprnames
write.csv(txi.kallisto$abundance,"reports/kallisto.TPM.csv")
write.csv(txi.kallisto$counts,"reports/kallisto.counts.csv")

# DEseq2 analyses

geno = factor (samples$sample)

sampleTable <- data.frame(genotype = geno)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable, ~ genotype )

#nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
#nrow(dds)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_kallisto.pdf")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:60]
labels <- as.data.frame(as_tibble(colData(dds))) %>% add_column(sampid=exprnames) %>% column_to_rownames(var="sampid")
colnames(labels) = c("genotype")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=labels,main="VSD")


pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=labels,main="RLD")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$Rep,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
dev.off()
