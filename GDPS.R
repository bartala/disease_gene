# # work plan
# [R]
# 1. Load the gene-disease network and Archs4 gene-gene similarity network
# 2. For each gene in Step 1 get genes with at least 0.5 similarity
# 3. Original network - create a network (edgelist) with two types of edges:
#   * gene-disease
# * gene-gene similarity


library(readr)
library(igraph)
library(sqldf)
library(dplyr)

PTH = "/path/to/data"

# ------------------------  Step 1 -------------------------------------------------------------

gene_disease <- read_csv(paste0(PTH,"edgelist_OMIM_Expanded.csv.gz"), col_types = cols(X1 = col_skip()))

# load human gene-gene correlation matrix from Archs4
load(paste0(PTH,"human_correlation.rda"))
genes <- row.names(cc) # 26,415 genes

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

# ------------------------  Step 3 -------------------------------------------------------------
# create edgelist
new_edges$type = '2' # meaning gene-gene similarity score
names(gene_disease)<-c("from","to")
gene_disease$type = '1' # meaning disease

edges <- rbind(gene_disease,new_edges)
write.csv(edges,file="~/zalon/gene_disease/data/combined_edges.csv")
