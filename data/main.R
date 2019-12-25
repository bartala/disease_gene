library(igraph)
PTH = '/users/alon/desktop/github/disease_gene/'

#---- load data --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# load ARCHS4 matrix
load(file.path(PTH, "data/human_correlation.rda"))

# load gene_disease edgelist from Leskovec
gene_disease = read(file.path(PTH, "data/DG-Miner_miner-disease-gene.tsv.gz"), compression = 'gzip')

# load gene_disease edgelist from enrichr
gene_disease_omn = read(file.path(PTH, "data/OMIM_Expanded.txt"), compression = 'gzip')


# --- crete a graph from ARCHS4 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
g = graph_from_adjacency_matrix(cc, mode = 'upper', weighted = TRUE, diag = FALSE )

# compute threshold
weights = E(g)$weight
stdv = sd(weights)
UCL = 0 + 0.5*stdv
LCL = 0 - 0.5*stdv


# delete edges
g1 = delete.edges(g, which( weights < UCL & weights > LCL ))