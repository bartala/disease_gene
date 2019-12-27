library(readr)
library(igraph)
library(sqldf)

# work plan
# 1. Load the gene-disease network and Archs4 gene-gene similarity network
# 2. For each gene in Step 1 get the genes with at least 0.5 similarity
# 3. Create a network (edgelist) with two types of edges: 1) gene-disease, and 2) gene-gene-similarity
# 4. Create positive and negative examples.
#    POS: 50% of gene-disease edges ---> delete them from the network in Step 3.
#    NEG: the same number of POS examples - randomly select pairs of node-disease that have no edges in Step 3.
# 5. Delete POS examples from the original network in Step 3.
# 6. Train Node2Vec on the new network in Step 5
# 7. Two predistion approches:
#     7.1 Node classification
#         Label nodes with disease affiliation
#         Get node embeddings from the Node2Vec model in 6
#         For each disease cluster- classify node embeddings as belong/not to the disease
#     7.2 Edge prediction
#         Get node embeddings from the Node2Vec model in 6
#         For each disease node- classify sum(embeddings gene, embeddings disease) as form/not

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
#    POS: Randomly sample 50% of edges from each gene-disease pair
library(dplyr)
df <- edges[edges$type==1,]
POS <- df %>% group_by(from) %>% sample_frac(0.2)
POS$type<-NULL

#    NEG: the same number of POS examples -
#       4.1) Randomly select N/2 pairs of gene-gene that have no edges in Step 3
#       4.2) Randomly select N/2 pairs of gene-disease that have no edges in Step 3

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

# step 4.2
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
write.csv(POS,file="~/zalon/gene_disease/data/POS_edges.csv")
write.csv(NEG,file="~/zalon/gene_disease/data/NEG_edges.csv")

# ------------------------  Step 5 -------------------------------------------------------------
# delete positive examples (POS edges) from the original network 
edges<-data.frame(edges)
POS<-data.frame(POS)
POS$type<-3 # positive gene-disease edge
train_edges <- rbind(edges,POS)
train_edges$type<-NULL
train_edges <- train_edges[!duplicated(train_edges,fromLast = FALSE)&!duplicated(train_edges,fromLast = TRUE),] 
write.csv(train_edges,file="~/zalon/gene_disease/data/train_edges.csv")

# number of diseases: 187
length(unique(train_edges[train_edges$type!=2 ,]$from))

# number of unique genes 2365
genes<-union(train_edges[train_edges$type!=2 ,]$to , c(train_edges[train_edges$type==1 ,]$to,train_edges[train_edges$type==1 ,]$from))
genes<-union(genes, c(train_edges[train_edges$type==3 ,]$to,train_edges[train_edges$type==3 ,]$from ))
length(unique(genes))

g<-graph_from_data_frame(train_edges,directed=FALSE)
components(g)
