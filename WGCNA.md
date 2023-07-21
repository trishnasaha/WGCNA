### WGCNA is an unsupervised method to cluster genes based on their expression profile. 

### WGCNA cluster genes that are functionally associated 

### Identified modules can be used for several downstrem analysis. For example

#### 1. Enrichment analysis
#### 2. Correlation of modules with phenotype or traits
#### 3. Identify driver gene or hub gene
#### 4. Identify regulatory network
#### 5. correlation between similar type of matrix

### GOOD practice criteria
#### 1. More number of samples (at least 20 samples)
#### 2. Genes should not be filtered by differential expression.
#### 3. For RNA-seq data:

- Remove all features that have consistently low counts
- Data should be normalized (variance-stabilizing transformation or log(RPKM/FPKM+1))
- Any normalization method can be used. But all sample need to processed the same way.
- If multiple Batch are present then batch effects need to be adjusted before running the analysis.