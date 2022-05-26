library(readr)
library(sqldf)
library(igraph)
library(data.table)
library(foreach)
library(dplyr)

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
# STEP3: combine the gene-disease graph with gene-gene edges from Step 2
# Create a Combined graph (g_d_t.csv) with two types of edges by adding gene-gene edges in Step #2 into the gene-disease
# graph in Step #1
# =========================================================================================================================

g_d_edges<-fread(paste0(PTH,"g_d_edges.csv"),header=TRUE)
g_d_edges$V1<-NULL
g_d_edges$type<-1
g_d_edges<-data.frame(g_d_edges)


gene_gene <- fread(paste0(PTH,"gene_gene_net_edges_normal_tissue.csv"),header=TRUE)
gene_gene$V1<-NULL
names(gene_gene)<-c("from","to","weight")
g_g_t$type<-2

# from g_g_t (tissue info) edges, keep only edges that their nodes (genes) are in gene-disease (g_d_edges)
genes<-unique(g_d_edges$geneSymbol)
keep_gg_edges<-g_g_t[g_g_t$from %in% genes | g_g_t$to %in% genes,]

# add edges to create gene-disease network with gene-gene edges from gene-tissue network
names(g_d_edges)<-names(keep_gg_edges)
g_d_t<-rbind(g_d_edges,keep_gg_edges)

write.csv(g_d_t,paste0(PTH,"g_d_t.csv"))

# =========================================================================================================================
# STEP4: Create positive (`POS`) and negative (`NEG`) examples.
# =========================================================================================================================

# Positive (`POS`) Examples: sample 20% of *gene-disease* edges (type 1) and delete them from the original network (g_d_t.csv) 

g_d_t <- fread(paste0(PTH,'g_d_t.csv'))
g_d_t$V1<-NULL

gene_disease <- g_d_t[g_d_t$type==1,c(1,2)]
gene_disease<-data.frame(gene_disease)
names(gene_disease)<-c("gene","disease")

diseases <- unique(gene_disease$disease)

# sample 20% of edges as positive examples
pos_edges<-data.frame()

for(disease in diseases){
  print(disease)
  df<-gene_disease[gene_disease$disease==disease,]
  N = round(0.2*nrow(df))
  # sample gene-edges
  edges<-df[sample(1:nrow(df),N),] # sample 20% edge
  pos_edges <- rbind(pos_edges,edges)
}

# number of diseases sampled
length(unique(pos_edges$disease))

# delete `pos` edges from the graph G to create train edges for node2vec embeddings
names(pos_edges)<-c("from","to")
edges_G_tag<-anti_join(g_d_t, pos_edges, by=c("from","to"))


write.csv(edges_G_tag,file=paste0(PTH,"G_tag.csv"), row.names = FALSE)
write.csv(pos_edges,file=paste0(PTH,"pos_edges.csv"), row.names = FALSE)


# `NEG`: sample the same number as `POS` examples - randomly select pairs of nodes that have no edges in Step 3.
#(1) 50% of gene-gene edges
#(2) 50% of gene-disease edges with low GDPS

N = round(nrow(pos_edges)/2)

# (1) ---------- 50% of negative gene-gene edges ----------

genes<- unique(gene_disease$gene)

neg_gene_gene_edges<-data.frame(
                                    from = sample(genes, N, replace = TRUE), 
                                    to = sample(genes, N, replace = TRUE) 
                                )

# keep only non existing gene-gene edges
neg_gene_gene_edges <- anti_join(neg_gene_gene_edges,g_d_t, by=c("from","to"))
n = N-nrow(neg_gene_gene_edges)


# (2)---------- 50% of negative gene-disease edges ---------->>>>>>>>

edges_G_tag <- read_csv(paste0(PTH,"G_tag.csv"))

#---- make sure we sample only low GDPS gene-disease scores ----
GDPS <- as.matrix(fread(file = paste0(PTH,"disease_prediction.tsv"), sep = "\t"),rownames=1) # gene X disease matrix
GDPS1 <- GDPS[,colnames(GDPS) %in% edges_G_tag$to] # keep only diseases in G'
GDPS1 <- GDPS1[rownames(GDPS1) %in% edges_G_tag$from, ] # keep only genes in G'
rm(GDPS)

diseases <- colnames(GDPS1)

# sample 20% of edges as positive examples
neg_edges<-data.frame()
i = 1
k=1
while(i < n){
  j = i
  if( i > length(diseases)){
    j = 1
  }
  disease = diseases[j]
  ordered_genes <- rownames(data.frame(sort((GDPS1[,disease]),decreasing = FALSE))) # get genes with low GDPS score
  df<- g_d_t[g_d_t$from == ordered_genes[k] & g_d_t$to == disease, ]
  if( nrow(df) == 0 ){
    neg_edges <- rbind(
                        neg_edges,
                        c(ordered_genes[1],disease)
                      )
    i=i+1
    print(i)
  }
  else{
    k=k+1
  }
}

# make sure only non existing edges are in neg_edges
names(neg_edges) <- c("from","to")
neg_edges <- anti_join(neg_edges,g_d_t, by=c("from","to"))

#---------- combine (gene,gene) and (gene, disease) pairs ----------
neg_edges<-rbind(neg_edges,neg_gene_gene_edges)

# make sure all genes and diseases are in G'
neg_edges<-neg_edges[union(neg_edges$from,neg_edges$to) %in% union(edges_G_tag$from,edges_G_tag$to),]

# make sure negative edges are not in G (`g_d_t`)
neg_edges <- anti_join(neg_edges,g_d_t, by=c("from","to"))

write.csv(neg_edges,file=paste0(PTH,"neg_edges.csv"), row.names = FALSE)


#---------- combine pos and neg examples ----------
neg_edges$type<-0
pos_edges$type<-1
names(pos_edges)<-names(neg_edges)
pos_neg_edges<-rbind(neg_edges,pos_edges)
write.csv(pos_neg_edges,file=paste0(PTH,"pos_neg_edges.csv"), row.names = FALSE)
