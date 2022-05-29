# 1. Load the gene-disease network and Archs4 gene-gene similarity network
# 2. For each gene in Step 1 get genes with at least 0.5 similarity

library(readr)
library(igraph)
library(sqldf)
library(dplyr)

PTH = "/path/to/data"

# ------------------------  Step 1 -------------------------------------------------------------

# Load G_tag
gene_disease <- read_csv(paste0(PTH,"G_tag.csv"))
gene_disease<-gene_disease[gene_disease$type == 1,]
names(gene_disease) = c("gene","disease")

# load human gene-gene correlation matrix from Archs4
load(paste0(PTH,"human_correlation.rda"))
genes <- row.names(cc)

# ------------------------  Step 2 -------------------------------------------------------------
# for each gene in gene_disease get the most similar genes from Archs4
new_edges <- data.frame()
genes_unique<- unique(gene_disease$gene)
archs4_genes<-unique(rownames(cc))
for(gene in genes_unique ) {
  print(gene)
  if(gene %in% archs4_genes){
    x = cc[gene,]
    x = x[(x>0.5 | x< -0.5)& x<1]
    if(length(x)!=0){
      similar_genes <- names(x)
      new_edges <- rbind(data.frame(from=names(x), to = gene ), new_edges)
    }
  }
}


