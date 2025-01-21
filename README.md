# HetIG-PreDiG: A Heterogeneous Integrated Graph Model for Predicting Human Disease Genes Based on Gene Expression

## Overview
This repository contains code and links to the datasets necessary to run the HetIG-PreDiG model. 
HetIG-PreDiG is a method for predicting novel disease genes using node embeddings in heterogeneous graphs.

### Motivation
Graph analytical approaches permit identifying novel genes involved in complex diseases, but are limited by (i) inferring structural network similarity of connected genes, ignoring potentially relevant unconnected nodes; (ii) using homogeneous graphs, missing gene-disease associations’ complexity; (iii) relying on disease/gene-phenotype associations’ similarities, involving highly incomplete data; (iv) using binary classification, with gene-disease edges as positive training samples, and non-associated gene and disease nodes as negative samples that may include currently unknown disease genes; or (v) reporting predicted novel associations without systematically evaluating their accuracy. 

### Results
HetIG-PreDiG analyzes a graph with gene-gene, gene-disease, and gene-tissue associations. It allows predicting novel disease genes using low-dimensional representation of nodes accounting for network structure, and extending beyond network structure using the developed Gene-Disease Prioritization Score (GDPS) reflecting the degree of gene-disease association via gene co-expression data.

## Running the code

`data2network.R`
  * Create gene-disease graph using [Dataset 1 (DS1)](https://www.disgenet.org/downloads)
  * Create gene-tissue graph and convert it into a gene-gene graph using [Dataset 2 (DS2)](https://www.proteinatlas.org/about/download)
  * Combine the gene-disease graph with gene-gene graph.
  * Create positive and negative training examples.


`GDPS.R` 
* Calculates the Gene-Disease prioritization score (GDPS) using the gene-gene correlation matrix from [Dataset 3 (DS3)](https://maayanlab.cloud/archs4/download.html)

`HetIG-PreDiG.py` 
* Learns node embeddings using node2vec.
* Train a logistic regression model to classify pairs of nodes into a link will/not form group. 



## Miscellaneous
Please send any questions you might have about the code and/or the algorithm to alon.bartal@biu.ac.il.

## Requirements
HetIG-PreDiG is tested to work under Python 3 and R 4.1.2 with RStudio 2021.09.2+382.

## Citing
If you find HetIG-PreDiG useful for your research, please consider citing us:
```
@article{jagodnik2023hetig,
  title={HetIG-PreDiG: A heterogeneous integrated graph model for predicting human disease genes based on gene expression},
  author={Jagodnik, Kathleen M and Shvili, Yael and Bartal, Alon},
  journal={Plos one},
  volume={18},
  number={2},
  pages={e0280839},
  year={2023},
  publisher={Public Library of Science San Francisco, CA USA}
}
```

