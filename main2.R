library(readr)
library(sqldf)
library(igraph)
library(data.table)
library(foreach)

# =========================================================================================================================
# STEP1: Create gene-disease network
# =========================================================================================================================

# load disease data from https://www.disgenet.org/downloads

g_d_edges<-read_delim("/users/alon/desktop/gd/all_gene_disease_associations.tsv.gz", "\t", escape_double = FALSE, trim_ws = TRUE)

# remove duplicates
g_d_edges[is.na(g_d_edges$YearInitial),'YearInitial']<-0
g_d_edges$YearInitial<-as.numeric(g_d_edges$YearInitial)

# keep edges from 2018 out of the picture
# e_d_2016<-g_d_edges[g_d_edges$YearFinal<=2016,]
# e_d_2017<-g_d_edges[g_d_edges$YearFinal==2017,]
# e_d_2018<-g_d_edges[g_d_edges$YearFinal==2018,]

# g_d_edges <- e_d_2016
# g_d_edges <- e_d_2017
# g_d_edges <- e_d_2018

x<-sqldf("select geneSymbol,diseaseId, count(1), min(YearInitial) as weight from g_d_edges group by geneSymbol,diseaseId")

g_d_edges<-x[,c("geneSymbol","diseaseId","weight")]

g_d_edges<-g_d_edges[,c("geneSymbol","diseaseId")] # Edges: 628,668 ; Genes 17,545 ; Diseases: 24,166
keep<-sqldf("select diseaseId,count(1) from g_d_edges group by diseaseId")
keep<-keep[keep$`count(1)`>10,]
g_d_edges<-g_d_edges[g_d_edges$diseaseId %in% keep$diseaseId,]
# Diseases: 6342; Genes: 17403


write.csv(g_d_edges,file = "/users/alon/desktop/gd/g_d_edges.csv")
#write.csv(g_d_edges,file = "/users/alon/desktop/e_d_2016.csv")
#write.csv(g_d_edges,file = "/users/alon/desktop/e_d_2017.csv")
#write.csv(g_d_edges,file = "/users/alon/desktop/e_d_2018.csv")


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
write.csv(gene_gene_net_edges,file="/users/alon/desktop/gene_gene_net_edges_normal_tissue.csv")
#fwrite(gene_gene_net_edges, "/users/alon/desktop/gene_gene_net_edges_normal_tissue.csv")

# =========================================================================================================================
# STEP3: enrich the gene-disease network with gene-gene edges from Step 2
# =========================================================================================================================

g_d_edges<-fread("/users/alon/desktop/gd/g_d_edges.csv",header=TRUE)
# g_d_edges<-fread("/users/alon/desktop/e_d_2016.csv",header=TRUE)
g_d_edges$V1<-NULL
g_d_edges$type<-1

gene_gene <- fread("/users/alon/desktop/gd/gene_gene_net_edges_normal_tissue.csv",header=TRUE)
gene_gene$V1<-NULL
names(gene_gene)<-c("from","to","weight")

# keep only gene-gene tissue edges with weight above 100
g_g_t<-gene_gene[gene_gene$weight>100,]
g_g_t$type<-2
g_g_t$weight<-NULL

# from g_g_t (tissue info) edges, keep only edges that their nodes (genes) are in gene-disease (g_d_edges)
genes<-unique(g_d_edges$geneSymbol)
keep_gg_edges<-g_g_t[g_g_t$from %in% genes | g_g_t$to %in% genes,]

# add edges to creare gene-disease network with gene-gene edges from gene-tissue network
g_d_edges<-data.table(g_d_edges)
names(g_d_edges)<-names(keep_gg_edges)
g_d_t<-rbind(g_d_edges,keep_gg_edges) # gene=disease:628668; gene-gene: 615903 edges

write.csv(g_d_t,"/users/alon/desktop/gd/g_d_t.csv")
# write.csv(g_d_t,"/users/alon/desktop/g_d_t_m3.csv")

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
diseases<-unique(gene_disease$disease) # 6342 diseases
clu <- components(g)


# for(disease in diseases){
#   print(disease)
#   df<-gene_disease[gene_disease$disease==disease,]
#   N = round(0.2*nrow(df))
#   smpedges<-data.frame()
#   # sample gene-edges
#   i=1
#   while(nrow(smpedges)<N){
#     print(i)
#     edge<-df[sample(1:nrow(df),1),] # sample 1 edge
#     new_edges <- gene_disease[gene_disease$gene!=edge$gene & gene_disease$disease!=edge$disease,] # all edges without the sampled edge
#     g_tmp<-graph_from_data_frame(new_edges)
#     if(components(g)$no <= clu$no){
#       smpedges <- rbind(smpedges,edge)
#       gene_disease <- new_edges
#       i=i+1
#     }
#   }
# }


# this is a try to sample 20% positive gene-disease edges while ignoring number of components
pos_edges<-foreach(i = 1:length(diseases),.combine = rbind)%do%{
  print(i)
  df<-gene_disease[gene_disease$disease==diseases[i],]
  N = round(0.2*nrow(df))
  edges<-df[sample(1:nrow(df),N),] # sample N edges
  return(edges)
}


# delete pos edges from the graph to create train edges for the N2V embeddings
library(dplyr)
names(pos_edges)<-c("from","to")
dis<-anti_join(g_d_t, pos_edges, by=c("from","to"))
write.csv(dis,file="/users/alon/desktop/gd/new_network.csv")


pos_edges<-pos_edges[union(pos_edges$from,pos_edges$to) %in% union(dis$from,dis$to),]

# pos_edges has 12,648 diseases (out of ~24k since some disease have only 1-2 gene)
# round(0.2*N=2) or round(0.2*N=1) is 0! so only diseases with more than 2 genes are sampled for pos examples
write.csv(pos_edges,file="/users/alon/desktop/gd/pos_edges.csv")


# `NEG`: sample the same number as `POS` examples - randomly select pairs of nodes that have no edges in Step 3.
#* 50% of gene-gene edges
#* 50% of gene-disease edges

#* 50% of gene-disease edges
N = round(nrow(pos_edges)/2)
genes<- unique(gene_disease$gene) # ~17k unique genes
neg_edges <- data.frame( from=sample(genes, N,replace = TRUE), to=sample(diseases, N,replace = TRUE) )
# keep only non existing edges
neg_edges <- anti_join(neg_edges,g_d_t, by=c("from","to"))

# check that NEG edges do not exist in real edges
tmp<-data.frame()
while(N-nrow(neg_edges)>nrow(tmp)+1){
  n = N-nrow(neg_edges)-nrow(tmp)
  tmp<-data.frame( from=sample(genes, n,replace = TRUE), to=sample(diseases, n,replace = TRUE) )
  tmp <- anti_join(tmp,g_d_t, by=c("from","to"))
  tmp <- anti_join(tmp,neg_edges, by=c("from","to"))
}
neg_edges<-rbind(neg_edges,tmp)
neg_edges<-neg_edges[union(neg_edges$from,neg_edges$to) %in% union(dis$from,dis$to),]
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
neg_edges<-neg_edges[union(neg_edges$from,neg_edges$to) %in% union(dis$from,dis$to),]

write.csv(neg_edges,file="/users/alon/desktop/gd/neg_edges.csv")


#4.3 Shuffle and split the examples (`training_testing_edges.csv`) into training (80%) and testing (20%)
neg_edges<-neg_edges <- read_csv("Desktop/gd/neg_edges.csv", col_types = cols(X1 = col_skip()))
pos_edges<-pos_edges <- read_csv("Desktop/gd/pos_edges.csv", col_types = cols(X1 = col_skip()))
neg_edges$type<-0
pos_edges$type<-1
names(pos_edges)<-names(neg_edges)
pos_neg_edges<-rbind(neg_edges,pos_edges)
rows <- sample(nrow(pos_neg_edges))
df <- pos_neg_edges[rows, ]
write.csv(df,file="/users/alon/desktop/gd/pos_neg_edges.csv")

#--------------------------------------- Statistics -------------------------------------------------------------------

# network
dis<-new_network[new_network$type==1,]
dis<-unique(dis$to) # 6342 diseases

genes<-new_network[new_network$type==1,]$from
genes<-c(genes,new_network[new_network$type==2,]$from)
genes<-c(genes,new_network[new_network$type==2,]$to)
genes<-unique(genes) # 17213 genes

# original files disgenet
g_d_edges<-read_delim("/users/alon/desktop/gd/all_gene_disease_associations.tsv.gz", "\t", escape_double = FALSE, trim_ws = TRUE)
length(unique(g_d_edges$diseaseId)) # disease 24166

# proteinatlas
normal_tissue_tsv <- read_delim("/users/alon/desktop/gd/normal_tissue.tsv.zip", "\t", escape_double = FALSE, trim_ws = TRUE)
genes<-unique(g_d_edges$geneSymbol) # 17545
genes<-union(genes, normal_tissue_tsv$`Gene name`) # 20751 genes

# Archs4
gene_names <- read_csv("archs4/gene_names.csv")
length(unique(gene_names$x)) # 26415

# genes in Archs4 but not in both datasets (9227)
table(gene_names$x %in% genes)
FALSE  TRUE 
9227 17188 

# genes in datasets but not in Archs4 (3563)
table(genes %in% gene_names$x)
FALSE  TRUE 
3563 17188 







x<-union(x, g_d_t_m3[g_d_t_m3$type==2,]$from)

