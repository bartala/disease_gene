# 1. Load the gene-disease network and Archs4 gene-gene similarity network
# 2. For each gene in Step 1 get genes with at least 0.5 similarity

library(readr)
library(igraph)
library(sqldf)
library(dplyr)

PTH = "/path/to/data"

# ------------------------ load data -------------------------------------------------------------

# Load G_tag
G <- read_csv(paste0(PTH,"g_d_t.csv"))
G <- G[G$type == 1,]
names(G) = c("gene","disease")


# load human gene-gene correlation matrix from Archs4
load(paste0(PTH,"human_correlation.rda"))
genes <- row.names(cc)

# ------------------------  GDPS -------------------------------------------------------------

GDPS <- function(disease){         
          z_tag <-  G_tag[ , G_tag$disease == disease,]$genes ] # select gene columns in cc
          GDPS_vec <- rowMeans(z_tag)
          return(GDPS_vec)
}



diseaeses <- unique(G$disease)

GDPS_matrix = data.frame()

for(disease in diseaeses){
  
        GDPS_vec <- GDPS(disease)
        GDPS_matrix <- cbind(GDPS_matrix,GDPS)
}

write.csv(GDPS_matrix, file = paste0(PTH,"gene_disease_GDPS_matrix"))




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


