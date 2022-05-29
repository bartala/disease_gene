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
