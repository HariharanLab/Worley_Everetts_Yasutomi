# Coded in R markdown 

# To create a gene signature of developmental progression / cellular maturity - we first looked for genes with differential expression during normal development
# (between mid and late 3rd instar)

# this set of genes is called "temporalGenes" 



```{r}
library(dplyr)
library(Seurat)
library(umap)
library(Matrix)
library(wordspace)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(ggcorrplot)
library(RColorBrewer)
library(viridis)
```

```{r}

# temporalGenes file was created by comparing the Epithelial cells from different time points in development (96h and 120h AEL, or mid and later 3rd instar)

temporalGenes = read.csv( file = "~/Epithelum_temporalGenes_265.csv", stringsAsFactors = FALSE)

temporalGenes = temporalGenes$x

# An example of a gene with differential expression during normal development: 
temporalGenes[temporalGenes == "br"]
```

# Next we ran Principal Component Analysis (PCA) on the single-cell data only considering these temporal genes.
# We investigating how the principal components seperated the data

# We converted the PC1 into a score of cellular maturity or developmental progression

```{r}


# This creats a new Seurat Object that runs PCA based on a set of genes that change during normal developmet: 

CellularMaturity = RunPCA( ScVI_intergratedObject , features = temporalGenes )

# this allows us to see how individual cells from the distint batches express these different temporal genes. 
# first we checkout how the PC1 and PC2 seperates the data:
DimPlot(CellularMaturity, reduction = "pca")

# Then we visualized the data in different ways to determine the relative cellular maturity score for distinct batchs and clusters 
RidgePlot(CellularMaturity, "PC_1", group.by = "Batch_Idents")
VlnPlot( CellularMaturity, "PC_1", pt.size = 0,  group.by = "Batch_Idents")
VlnPlot( CellularMaturity, "PC_1", pt.size = 0,  group.by = "named_cluster_ID")

```






