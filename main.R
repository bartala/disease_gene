# work plan
# 1. Load the gene-disease network and Archs4 gene-gene similarity network
# 2. For each gene in Step 1 get the genes with at least 0.5 similarity
# 3. Create a network (edgelist) with two types of edges: 1) gene-disease, and 2) gene-gene-similarity
# 4. Create positive and negative examples.
#    4.1 POS: sample 20% of gene-disease edges and delete them from the original network (Step 3)
#    4.2 NEG: sample the same number of POS examples - randomly select pairs of node-disease that have no edges in Step 3.
# 5. Train Node2Vec on the new network in Step 5 [python]
# 6. Two predistion approches:
#     6.1 Node classification
#         Label nodes with disease affiliation
#         Get node embeddings from the Node2Vec model in 5
#         For each disease cluster- classify node embeddings as belong/not to the disease
#     6.2 Edge prediction
#         Get node embeddings from the Node2Vec model in 5
#         For each disease node - classify sum(embeddings gene, embeddings disease) as link (i.e. type=1) or not (i.e. type=0)


library(readr)
library(igraph)
library(sqldf)
library(dplyr)

# ------------------------  Step 1 -------------------------------------------------------------
# load the human gene-gene correlation matrix from Archs4
load("~/zalon/gene_disease/data/human_correlation.rda")
genes <- row.names(cc) # 26,415 genes

# load gene-disease data from enrichr
# 187 diseases, 2178 genes, 16,628 edges
gene_disease <- read_csv("zalon/gene_disease/data/edgelist_OMIM_Expanded.csv.gz", col_types = cols(X1 = col_skip()))

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

# ------------------------  Step 4 -------------------------------------------------------------
# Create positive and negative examples.

#    4.1 POS: Randomly sample 20% of edges from each gene-disease pair BUT keep the graph connected!

g<-graph_from_data_frame(edges[,c("from","to")])
clu<-components(g)
POS<-data.frame()
df <- edges[edges$type==1,]
disease<-unique(df$from)
for(i in 1:length(disease)){
  print(i)
  e<-get.data.frame(g)
  e<-e[e$from==disease[i],]
  k = as.integer(0.2*nrow(e)) # sample 20% of edges from each disease
  for(j in 1:k){
    eg<-e[sample(nrow(e),1),]
    g_tmp<-delete_edges(g, paste0(eg$from[1],"|",eg$to[1]))
    while(components(g_tmp)$no != clu$no){
      eg<-e[sample(nrow(e),1),]
      g_tmp<-delete_edges(g, paste0(eg$from[1],"|",eg$to[1]))
    }
    g<-delete_edges(g, paste0(eg$from[1],"|",eg$to[1]))
    e<-get.data.frame(g)
    e<-e[e$from==disease[i],]
    POS<-rbind(POS,eg)
  }
}

# save the graph for training (187 diseases and 414 genes)
edgelist_node2vec <- get.data.frame(g)
write.csv(edgelist_node2vec,file="~/zalon/gene_disease/data/edgelist_node2vec.csv")
write.csv(edges,file="~/zalon/gene_disease/data/original_edges.csv")


#    4.2 NEG: the same number of POS examples
#       4.2.1) Randomly select N/2 pairs of gene-gene that have no edges in Step 3
#       4.2.2) Randomly select N/2 pairs of gene-disease that have no edges in Step 3

#       sub-step 4.2.1 ------------------------------------------------------- 
N = nrow(POS)
only_genes <- union(edges[edges$type==2,]$from, edges[edges$type==2,]$to)
only_genes <- union(only_genes, edges[edges$type==1,]$to)

# step 4.1
NEG<-data.frame()
df1<-df2<-data.frame()
flag = TRUE
for(i in 1:round(N/2)){
  print(i)
  while((nrow(df1)>0) | (nrow(df2)>0) | flag){
    flag = FALSE
    fr <- sample(only_genes, 1)
    to <- sample(only_genes[only_genes!=fr], 1)
    df1 <- edges[(edges$from == fr & edges$to == to) | (edges$from == to & edges$to == fr),]
    df2<- NEG[(NEG$from == fr & NEG$to == to) | (NEG$from == to & NEG$to == fr),]
  }
  flag = TRUE
  NEG<-rbind(NEG,data.frame(from=fr,to=to))
}

#     step 4.2  ------------------------------------------------------------------------------------
df = edges[edges$type==1,] # take only gene-disease edges
NEG_df <- data.frame()
N = nrow(POS)-nrow(NEG)
df1<-df2<-data.frame()
diseases <- unique(df$from)
gene_disease_false <- unique(df$to)
flag = TRUE
for(i in 1:N){
  print(i)
  while((nrow(df1)>0) | (nrow(df2)>0) | flag){
    flag = FALSE
    fr <- sample(diseases, 1) # sample disease
    to <- sample(gene_disease_false, 1) # sample a gene that is affiliated with some disease/s
    df1 <- df[(df$from == fr & df$to == to) | (df$from == to & df$to == fr),]
    df2<- NEG_df[(NEG_df$from == fr & NEG_df$to == to) | (NEG_df$from == to & NEG_df$to == fr),]
  }
  flag = TRUE
  NEG_df<-rbind(NEG_df,data.frame(from=fr,to=to))
}

# combin negative examples
NEG<-rbind(NEG,NEG_df)

# save positive and negative examples
POS$type<-1
NEG$type<-0

Training_testing<- rbind(POS,NEG)
# suffle rows
Training_testing<-Training_testing[sample(nrow(Training_testing)),]
Training_testing[Training_testing$from %in% ,]$is_disease


write.csv(Training_testing,file="~/zalon/gene_disease/data/training_testing_edges.csv")
