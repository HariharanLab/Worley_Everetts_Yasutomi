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
temporalGenes = read.csv( file = "~/Epithelum_temporalGenes_265.csv", stringsAsFactors = FALSE)

temporalGenes = temporalGenes$x

# An example of a gene with differential expression during normal development: 

temporalGenes[temporalGenes == "br"]
```

# Next we ran Principal Component Analysis (PCA) on the single-cell data only considering these temporal genes.
# We investigating how the principal components seperated the data

# We converted the PC1 into a score of cellular maturity or developmental progression

