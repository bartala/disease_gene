library(igraph)
library(readr)

PTH = '/users/alon/desktop/github/disease_gene/'

#---- load data --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# load ARCHS4 matrix
load(file.path(PTH, "data/human_correlation.rda"))

# load gene_disease edgelist from Leskovec
gene_disease <- read_delim("Desktop/GitHub/disease_gene/data/DG-Miner_miner-disease-gene.tsv.gz", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

# read gene disease 
gmt<-GSA.read.gmt('/users/alon/desktop/github/disease_gene/data/OMIM_Expanded.gmt')
# create edgelist for gene-disease
edges<-data.frame()
for(i in 1:length(gmt$geneset.names)){
  print(i)
  df<-data.frame(disease = gmt$geneset.names[i], gene = gmt$genesets[[i]])
  edges<-rbind(edges,df)
}
edges<-edges[edges$gene!="",]
edges<-edges[edges$disease!="",]
z <- gzfile(file.path(PTH, "data/edgelist_OMIM_Expanded.csv.gz"))
write.csv(edges, z )

G <- graph.data.frame(edges,directed=FALSE)
A <- as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)


# load gene_disease edgelist from enrichr
gene_disease_omn = read(file.path(PTH, "data/OMIM_Expanded.txt"), compression = 'gzip')



# --- crete a graph from ARCHS4 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# compute threshold
stdv = sd(cc)
UCL = 0 + 0.5*stdv
LCL = 0 - 0.5*stdv
cc[cc>LCL & cc<UCL]=0

# delete edges with low correlation
g<-graph_from_adjacency_matrix(cc, mode = 'upper', diag = FALSE, weighted = TRUE)
rm(cc)
# get edgelist
edgelist_ARCHS4 = get.data.frame(g)
z <- gzfile(file.path(PTH, "data/edgelist_ARCHS4.csv.gz"))
write.csv(edgelist_ARCHS4, z )

# for gene in V(g): get corresponding adjacency vector from gene_disease





