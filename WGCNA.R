# script to perform WGCNA

# Display current working directory
getwd();

# workingDir = "/Users/trishna/Dropbox/00_DROPBOX/GitHub/WGCNA";
# setwd(workingDir) ;

# load packages

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

# The following setting is important, do not omit.

# options(stringsAsFactors = FALSE);

# skip multithreding line if you run RStudio or other third-party R environments

# enableWGCNAThreads()   # allow multithreading(optional)

# fetch data

data <- read.delim("/Users/trishna/Dropbox/00_DROPBOX/GitHub/WGCNA/cts_merged.tsv", sep ="\t", header = T )

# QC - outlier detection
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove outlier genes
data1 <- data[gsg$goodGenes == TRUE, ]

# detect outlier samples - hierarchical clustering - method1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# pca - method2

pca <- prcomp(t(data1))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.var.percent

pca.dat <- as.data.frame(pca.dat)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2:', pca.var.percent[2], ' %'))

# exclude outlier sample

samples.to.be.excluded <- c('SRR13147312', 'SRR13147316', 'SRR13147318')
data.subset <- data1[,!(colnames(data1) %in% samples.to.be.excluded)]

df <- data.subset + 1
# Normalization
# create deseq2 dataset
coldata_merged <- read.delim("/Users/trishna/Dropbox/00_DROPBOX/GitHub/WGCNA/coldata_merged.tsv", sep ="\t", header = T )

colData <- coldata_merged %>%
  filter(!row.names(.) %in% samples.to.be.excluded)

names(colData)

# making sure the rownames and column names are identical 
all(rownames(colData) %in% colnames(df))
all(rownames(colData) == colnames(df))

# create Deseq2 dataset

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = colData,
                              design = ~ condition) #not specifying model

# remove all genes with counts < 15 in more than 75% of the samples (num. of obs in pheno or coldata*0.75 = number of samples must have >= 15 counts)
#suggested by WGCNA on RNAseq FAQ
#(240*0.75 = 180)

dds75 <- dds[rowSums(counts(dds) >= 5) >= 180,]
nrow(dds75)

# perform variance stabilization

dds_norm <- varianceStabilizingTransformation(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>%
  t()

# soft threshloding
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <- sft$fitIndices
# visulization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

#convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 5
temp_cor <- cor
cor <- WGCNA::cor

#memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)
cor <- temp_cor

#module eigengenes
module_eigengenes <-bwnet$MEs

#print a preview
head(module_eigengenes)

#get number of genes for each module
table(bwnet$colors)

#plot the dendogram and module colours before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


