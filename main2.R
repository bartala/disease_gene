library(readr)
library(sqldf)
library(igraph)
library(data.table)
library(foreach)
library(dplyr)

# =========================================================================================================================
# STEP1: Create gene-disease network
# =========================================================================================================================

# load gene-disease data from https://www.disgenet.org/downloads

g_d_edges<-read_delim("/users/alon/desktop/gd/all_gene_disease_associations.tsv.gz", "\t", escape_double = FALSE, trim_ws = TRUE)
g_d_edges<-g_d_edges[,c("geneSymbol","diseaseId")]
# remove duplicates
g_d_edges<-sqldf("select geneSymbol,diseaseId, count(1) as weight from g_d_edges group by geneSymbol,diseaseId")
g_d_edges<-g_d_edges[,c("geneSymbol","diseaseId")] # Edges: 628,668 ; Genes 17,545 ; Diseases: 24,166
keep<-sqldf("select diseaseId,count(1) from g_d_edges group by diseaseId")
keep<-keep[keep$`count(1)`>100,] # keep only dieseases that have at least 100 genes
length(unique(keep$diseaseId)) # number of diseases
g_d_edges<-g_d_edges[g_d_edges$diseaseId %in% keep$diseaseId,]
# Diseases: 6342; Genes: 17403
write.csv(g_d_edges,file = "/users/alon/desktop/gd/g_d_edges.csv")

# =========================================================================================================================
# STEP2: gene-tissue --> gene-gene network
# =========================================================================================================================

# load gene-tissue data from https://www.proteinatlas.org/about/download

normal_tissue_tsv <- read_delim("/users/alon/desktop/gd/normal_tissue.tsv.zip", "\t", escape_double = FALSE, trim_ws = TRUE)
tis<-normal_tissue_tsv[normal_tissue_tsv$Reliability %in% c('Approved','Enhanced'),]
tis<-tis[tis$Level %in% c( 'High','Medium'),]
g<-graph_from_data_frame(tis[,c('Gene name', 'Tissue')])
g2M <-g
V(g2M)$type <- bipartite_mapping(g2M)$type 
g_el <- as_edgelist(g)
colnames(g_el) <- c("gene", "tissue")
V(g)$type <- ifelse(V(g)$name %in% g_el[,"gene"], TRUE, FALSE)
projected_g <- bipartite_projection(g, multiplicity = TRUE)
gene_gene_net <- projected_g$proj2 # nodes:10722, edges:52686716
gene_gene_net_edges <- cbind( get.edgelist(gene_gene_net) , E(gene_gene_net)$weight )
gene_gene_net_edges<-data.table(gene_gene_net_edges)
write.csv(gene_gene_net_edges,file="/users/alon/desktop/gd/gene_gene_net_edges_normal_tissue.csv")
#fwrite(gene_gene_net_edges, "/users/alon/desktop/gene_gene_net_edges_normal_tissue.csv")

# =========================================================================================================================
# STEP3: enrich the gene-disease network with gene-gene edges from Step 2
# =========================================================================================================================

g_d_edges<-fread("/users/alon/desktop/gd/g_d_edges.csv",header=TRUE)
g_d_edges$V1<-NULL
g_d_edges$type<-1

gene_gene <- fread("/users/alon/desktop/gd/gene_gene_net_edges_normal_tissue.csv",header=TRUE)
gene_gene$V1<-NULL
names(gene_gene)<-c("from","to","weight")
# keep only gene-gene tissue edges with weight above 100
g_g_t<-gene_gene[gene_gene$weight>100,]
g_g_t$type<-2
g_g_t$weight<-NULL

# from g_g_t (tissue info) edges, keep only edges that at least 1 node of an edge (genes) is in gene-disease (g_d_edges)
genes<-unique(g_d_edges$geneSymbol)
keep_gg_edges<-g_g_t[g_g_t$from %in% genes | g_g_t$to %in% genes,]

# add edges to creare gene-disease network with gene-gene edges from gene-tissue network
g_d_edges<-data.table(g_d_edges)
names(g_d_edges)<-names(keep_gg_edges)
g_d_t<-rbind(g_d_edges,keep_gg_edges) # gene=disease:628668; gene-gene: 615903 edges

write.csv2(g_d_t,"/users/alon/desktop/gd/g_d_t.csv")

genes<- c(g_d_t[g_d_t$type==1,]$from, g_d_t[g_d_t$type==2,]$to, g_d_t[g_d_t$type==2,]$from)
genes<-unique(genes)

# =========================================================================================================================
# STEP4: Create positive (`POS`) and negative (`NEG`) examples.
# =========================================================================================================================

# 4.1 `POS` Examples: sample 20% of *gene-disease* edges (type 1) and delete them from the original network (g_d_t.csv) 
# while ensuring that the original network obtained after edge removals is connected
g_d_t <- fread('/users/alon/desktop/gd/g_d_t.csv')
g_d_t$V1<-NULL

gene_disease <- g_d_t[g_d_t$type==1,c(1,2)]
gene_disease<-data.frame(gene_disease)
names(gene_disease)<-c("gene","disease")
g<-graph_from_data_frame(gene_disease,directed = FALSE)
diseases<-unique(gene_disease$disease) # 24,166 diseases
clu <- components(g)

# this is a try to sample 20% positive gene-disease edges while ignoring number of components
pos_edges<-foreach(i = 1:length(diseases),.combine = rbind)%do%{
  print(i)
  df<-gene_disease[gene_disease$disease==diseases[i],]
  N = round(0.5*nrow(df))
  edges<-df[sample(1:nrow(df),N),] # sample N edges
  return(edges)
}

names(pos_edges)<-c("from","to")

# delete pos edges from the graph to create train edges for the N2V embeddings
new_network<-anti_join(g_d_t, pos_edges, by=c("from","to"))

# remove only edges that their removal their nodes are still in the new network
pos_edges<-pos_edges[pos_edges$from %in% union(new_network$from,new_network$to) & 
                       pos_edges$to %in% union(new_network$from,new_network$to),]

new_network<-anti_join(g_d_t, pos_edges, by=c("from","to"))
table( union(pos_edges$from, pos_edges$to) %in% union(new_network$from,new_network$to))

write.csv(pos_edges,file="/users/alon/desktop/gd/pos_edges.csv")
write.csv(new_network,file="/users/alon/desktop/gd/new_network.csv")


# `NEG`: sample the same number as `POS` examples - randomly select pairs of nodes that have no edges in Step 3.
#* 50% of gene-gene edges
#* 50% of gene-disease edges

#* 50% of gene-disease edges
new_network<-fread('/users/alon/desktop/gd/new_network.csv',header = TRUE)
names(new_network)<-c("index","from","to")
genes<-genes[genes %in% union(new_network$from, new_network$to)]

N = round(nrow(pos_edges)/2)
neg_edges <- data.frame( from=sample(genes, N,replace = TRUE), to=sample(diseases, N,replace = TRUE) )
# keep only non existing edges
neg_edges$from<-as.character(neg_edges$from)
neg_edges$to<-as.character(neg_edges$to)
library(dplyr)
neg_edges <- anti_join(neg_edges,g_d_t, by=c("from","to"))

# check that NEG edges do not exist in real edges
tmp<-data.frame()
while(N-nrow(neg_edges)>nrow(tmp)){
  n = N-nrow(neg_edges)-nrow(tmp)
  print(n)
  tmp<-data.frame( from=sample(genes, n,replace = TRUE), to=sample(diseases, n,replace = TRUE) )
  tmp$from<-as.character(tmp$from)
  tmp$to<-as.character(tmp$to)
  tmp <- anti_join(tmp,g_d_t, by=c("from","to"))
  tmp <- anti_join(tmp,neg_edges, by=c("from","to"))
}
neg_edges<-rbind(neg_edges,tmp)
neg_edges <- anti_join(neg_edges,g_d_t, by=c("from","to"))


#* 50% of gene-gene edges

neg_gene_gene_edges<-data.frame( from=sample(genes, N,replace = TRUE), to=sample(genes, N,replace = TRUE) )
neg_gene_gene_edges <- anti_join(neg_gene_gene_edges,g_d_t, by=c("from","to"))
neg_gene_gene_edges <- anti_join(neg_gene_gene_edges,neg_edges, by=c("from","to"))
n = N-nrow(neg_gene_gene_edges)
tmp<-data.frame( from=sample(genes, n,replace = TRUE), to=sample(genes, n,replace = TRUE) )
tmp <- anti_join(tmp,g_d_t, by=c("from","to"))
tmp <- anti_join(tmp,neg_edges, by=c("from","to"))
neg_gene_gene_edges <- rbind(neg_gene_gene_edges,tmp)

neg_edges<-rbind(neg_edges,neg_gene_gene_edges)
neg_edges <- anti_join(neg_edges,g_d_t, by=c("from","to"))

table( union(neg_edges$from, neg_edges$to) %in% union(new_network$from,new_network$to))
write.csv(neg_edges,file="/users/alon/desktop/gd/neg_edges.csv")


#4.3 Shuffle and split the examples (`training_testing_edges.csv`) into training (80%) and testing (20%)
pos_neg_edges<-rbind(neg_edges,pos_edges)
rows <- sample(nrow(pos_neg_edges))
df <- pos_neg_edges[rows, ]

# test that all genes are in the gene_disease_tissue network
table( union(df$from, df$to) %in% union(new_network$from,new_network$to))

write.csv(a,file="/users/alon/desktop/gd/pos_neg_edges.csv")



