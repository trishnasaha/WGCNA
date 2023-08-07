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



