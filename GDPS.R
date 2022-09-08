library(readr)

PTH = "/path/to/data"

# Load G
G <- read_csv(paste0(PTH,"G_tag.csv")) # or load d_g_t.csv
G <- G[G$type == 1,] # gene-disease edges
names(G) = c("gene","disease","label")


# Load human gene-gene correlation matrix from Archs4
load(paste0(PTH,"human_correlation.rda")) # this loades a `dataset` object
genes_archs4 <- names(dataset)


GDPS <- function(disease){    
  genes <- unique(G[G$disease == disease,]$gene)
  genes <- genes[genes %in% genes_archs4]
  z_tag <- dataset[ ,genes] # select gene columns of the gene-gene matrix
  z_tag <- data.frame(z_tag)
  GDPS_vec <- rowMeans(z_tag)
  return(GDPS_vec)
}



diseaeses <- unique(G$disease)

GDPS_matrix <- list()
i = 1
for(disease in diseaeses){
  
  GDPS_matrix[[i]] <- GDPS(disease)
  i <- i+1
  print(i)
  
}

GDPS_matrix <- data.frame(GDPS_matrix)
row.names(GDPS_matrix) <- genes_archs4
names(GDPS_matrix) <- diseaeses

write_tsv(GDPS_matrix, file = paste0(PTH,"disease_prediction.tsv.gz"))
