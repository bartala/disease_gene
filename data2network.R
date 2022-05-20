library(readr)
library(sqldf)
library(igraph)
library(data.table)
library(foreach)

PTH = "/path/to/data/" # enter your data path

# =========================================================================================================================
# STEP1: Create gene-disease graph (using DS1)
# =========================================================================================================================

# load disease data from https://www.disgenet.org/downloads

g_d_edges<-read_delim(paste0(PTH,"all_gene_disease_associations.tsv.gz"), "\t", escape_double = FALSE, trim_ws = TRUE)

g_d_edges<-g_d_edges[g_d_edges$diseaseType == 'disease',]
g_d_edges<-g_d_edges[g_d_edges$diseaseSemanticType == 'Disease or Syndrome',]

# update missing year values
g_d_edges[is.na(g_d_edges$YearInitial),'YearInitial']<-0
g_d_edges$YearInitial<-as.numeric(g_d_edges$YearInitial)

# create an edge-list
x<-sqldf("select geneSymbol,diseaseId, count(1) as weight, min(YearInitial) from g_d_edges group by geneSymbol,diseaseId")

g_d_edges<-x[,c("geneSymbol","diseaseId","weight")]
g_d_edges<-g_d_edges[,c("geneSymbol","diseaseId")]

# remove diseases with fewer than 2 genes
keep<-sqldf("select diseaseId, count(1) from g_d_edges group by diseaseId")
keep<-keep[keep$`count(1)`>2,]
g_d_edges<-g_d_edges[g_d_edges$diseaseId %in% keep$diseaseId,]

write.csv(g_d_edges,file = paste0(PTH,"g_d_edges.csv"))

# =========================================================================================================================
# STEP2: gene-tissue --> gene-gene graph (using DS2)
# =========================================================================================================================

# load gene-tissue data from https://www.proteinatlas.org/about/download
# load gene-tissue data from https://www.proteinatlas.org/about/download

library(readr)
library(sqldf)
library(igraph)
library(data.table)
library(foreach)

PTH = '/home/bartalab/Desktop/'

normal_tissue_tsv <- read_delim(paste0(PTH,"normal_tissue.tsv.zip"), "\t", escape_double = FALSE, trim_ws = TRUE)

tis<-normal_tissue_tsv[normal_tissue_tsv$Reliability %in% c('Approved'),]

tis<-tis[tis$Level %in% c( 'High'),]
tis<-data.frame(tis)

# gene-disease associations with more than a single cell type
tis <- sqldf("select *, count(1) as weight from tis group by `Gene.name`, Tissue having weight > 1")

length(unique(tis$Gene.name))
length(unique(tis$Tissue))

g<-graph_from_data_frame(tis[,c('Gene.name', 'Tissue')])
g2M<-g
V(g2M)$type <- bipartite_mapping(g2M)$type 

g_el <- as_edgelist(g)

colnames(g_el) <- c("gene", "tissue")

V(g)$type <- ifelse(V(g)$name %in% g_el[,"gene"], TRUE, FALSE)

projected_g <- bipartite_projection(g, multiplicity = TRUE)

gene_gene_net <- projected_g$proj2

gene_gene_net_edges <- cbind( get.edgelist(gene_gene_net) , E(gene_gene_net)$weight )

gene_gene_net_edges<-data.table(gene_gene_net_edges)

gene_gene_net_edges$V3 <- as.numeric(gene_gene_net_edges$V3)

write.csv(gene_gene_net_edges,file=paste0(PTH,"gene_gene_net_edges_normal_tissue.csv"))

# =========================================================================================================================
# STEP3: enrich the gene-disease network with gene-gene edges from Step 2
# Create a Combined Network (g_d_t.csv) with two types of edges by adding gene-gene edges in Step #2 into the gene-disease
# network in Step #1
# =========================================================================================================================

g_d_edges<-fread(paste0(PTH,"g_d_edges.csv"),header=TRUE)
g_d_edges$V1<-NULL
g_d_edges$type<-1  # 13k unique genes
g_d_edges<-data.frame(g_d_edges)


gene_gene <- fread(paste0(PTH,"gene_gene_net_edges_normal_tissue.csv"),header=TRUE)
gene_gene$V1<-NULL
names(gene_gene)<-c("from","to","weight") # 1659 unique genes
g_g_t$type<-2

# from g_g_t (tissue info) edges, keep only edges that their nodes (genes) are in gene-disease (g_d_edges)
genes<-unique(g_d_edges$geneSymbol)
keep_gg_edges<-g_g_t[g_g_t$from %in% genes | g_g_t$to %in% genes,]

# add edges to create gene-disease network with gene-gene edges from gene-tissue network
names(g_d_edges)<-names(keep_gg_edges)
g_d_t<-rbind(g_d_edges,keep_gg_edges)

write.csv(g_d_t,paste0(PTH,"g_d_t.csv"))


length(unique(g_d_t[g_d_t$type==1,]$to))
length(union(g_d_t[g_d_t$type==2,]$to, g_d_t$from))

# =========================================================================================================================
# STEP4: Create positive (`POS`) and negative (`NEG`) examples.
# =========================================================================================================================

# 4.1 `POS` Examples: sample 20% of *gene-disease* edges (type 1) and delete them from the original network (g_d_t.csv) 
# while ensuring that the original network obtained after edge removals is connected
g_d_t <- fread(paste0(PTH,'g_d_t.csv'))
g_d_t$V1<-NULL

gene_disease <- g_d_t[g_d_t$type==1,c(1,2)]
gene_disease<-data.frame(gene_disease)
names(gene_disease)<-c("gene","disease")
g<-graph_from_data_frame(gene_disease,directed = FALSE)
diseases<-unique(gene_disease$disease)
clu <- components(g)

# sample positive edges
for(disease in diseases){
  print(disease)
  df<-gene_disease[gene_disease$disease==disease,]
  N = round(0.2*nrow(df))
  smpedges<-data.frame()
  # sample gene-edges
  i=1
  while(nrow(smpedges)<N){
    print(i)
    edge<-df[sample(1:nrow(df),1),] # sample 1 edge
    new_edges <- gene_disease[gene_disease$gene!=edge$gene & gene_disease$disease!=edge$disease,] # all edges without the sampled edge
    g_tmp<-graph_from_data_frame(new_edges)
    if(components(g)$no <= clu$no){
      smpedges <- rbind(smpedges,edge)
      gene_disease <- new_edges
      i=i+1
    }
  }
}


# delete pos edges from the graph to create train edges for the N2V embeddings
library(dplyr)
names(pos_edges)<-c("from","to")
dis<-anti_join(g_d_t, pos_edges, by=c("from","to"))
write.csv(dis,file=paste0(PTH,"/new_network.csv"))


pos_edges<-pos_edges[union(pos_edges$from,pos_edges$to) %in% union(dis$from,dis$to),]

# round(0.2*N=2) or round(0.2*N=1) is 0! so only diseases with more than 2 genes are sampled for pos examples
write.csv(pos_edges,file=paste0(PTH,"pos_edges.csv"))


# `NEG`: sample the same number as `POS` examples - randomly select pairs of nodes that have no edges in Step 3.
#* 50% of gene-gene edges
#* 50% of gene-disease edges

#* 50% of gene-disease edges
N = round(nrow(pos_edges)/2)
genes<- unique(gene_disease$gene)
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


# 50% of gene-gene edges
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

write.csv(neg_edges,file=paste0(PTH,"neg_edges.csv"))


#4.3 Shuffle and split the examples (`training_testing_edges.csv`) into training (70%) and testing (30%)
neg_edges<-neg_edges <- read_csv(paste0(PTH,"neg_edges.csv"), col_types = cols(X1 = col_skip()))
pos_edges<-pos_edges <- read_csv(paste0(PTH,"pos_edges.csv"), col_types = cols(X1 = col_skip()))
neg_edges$type<-0
pos_edges$type<-1
names(pos_edges)<-names(neg_edges)
pos_neg_edges<-rbind(neg_edges,pos_edges)
rows <- sample(nrow(pos_neg_edges))
df <- pos_neg_edges[rows, ]
write.csv(df,file=paste0(PTH,"pos_neg_edges.csv"))
