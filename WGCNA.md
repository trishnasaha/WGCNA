### Weighted Gene Co-expression Network Analysis (WGCNA) is an unsupervised method to cluster genes based on their expression profile. 

### WGCNA cluster genes that are functionally associated in same pathway or regulatory mechanism, part of the same complex or influence each other or may be influenced by the same underlying mechanism.

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

### Hard threshold

- unweighted network
- It only idicates if a pair of genes are connected or not connected.
- It does not indicate how strong or weak the connection is between a pair of genes.

### Soft Threshold

- Weighted network
-correlation between pair of genes are raised to a power term. 

- R package has pickSoft threshold() function. It pick the appropriate power term. It selects the ideal threshold/power term that results a network with scale free topology.

- scale free network identify few highly connected nodes also known as hub, that intereact with lot of other nodes.

- PPI network and gene regulatory netowrk are scale free network.


### network proximity measure for clustering gene.

- measure of similarity = topological overlap measure (TOM)

- dissimilarity value is the distance of a gene from everyother gene in the system by substracting
the value from 1
1-TOM 

- Linkage hierarchical clustering methods use the similarity and dissimilarity measure to construct gene dendogram 

- cuttreeDynamic()
to calculate pairwise eigengene of similar module and calculate dissimilarity measurement by 1- cor of similarity eigengene.

- add line at the hieght of 0.25, which corresponds to correlation value of 75% 

- any module greater than 75% correlation are correlated so we can merge them.

### module trait association

- Module eigengene is the standardized gene expression profile for a given module.

- we can find correlation between module eigengene( ME ) and external traits/phenotypes

- hub or driver gene can be identified using module membership or highly connected genes present in the module. this can be calculated based on the correlation between expressions of the genes with the module eigengene. gene with high module membership is likely to be hub gene.

### R packages:
- WGCNA
- DESeq2
- GEOquery
- tidyverse
- CorLevelPlot
- gridExtra


### Resources

https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

https://peterlangfelder.com/articles/


