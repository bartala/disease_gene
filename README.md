# HetIG-PreDiG: A Heterogeneous Integrated Graph Model for Predicting Human Disease Genes Based on Gene Expression

## Abstract
Graph analytical approaches permit identifying novel genes involved in complex diseases, but are limited by (i) inferring structural network similarity of connected genes, ignoring potentially relevant unconnected nodes; (ii) using homogeneous graphs, missing gene-disease associations’ complexity; (iii) relying on disease/gene-phenotype associations’ similarities, involving highly incomplete data; (iv) using binary classification, with gene-disease edges as positive training samples, and non-associated gene and disease nodes as negative samples that may include currently unknown disease genes; or (v) reporting predicted novel associations without systematically evaluating their accuracy. Addressing these limitations, we develop the Heterogeneous Integrated Graph for Disease Genes Prediction model (HIGDGP) that includes gene-gene, gene-disease, and gene-tissue associations. We predict novel disease genes using low-dimensional representation of nodes accounting for network structure, and extending beyond network structure using the developed Gene-Disease Prioritization Score (GDPS) reflecting the degree of gene-disease association via gene co-expression data. For negative training samples, we select non-associated gene and disease nodes with lower GDPS. Lastly, we evaluate the developed model’s success in predicting novel disease genes by analyzing the prediction probabilities of gene-disease associations. HIGDGP successfully predicts (Micro-F1 = 0.91) gene-disease associations, outperforming baseline models, and is validated using published literature, thus advancing our understanding of genetic diseases.

## Running the code


## Miscellaneous
Please send any questions you might have about the code and/or the algorithm to alon.bartal@biu.ac.il.

## Requirements
HetIG-PreDiG is tested to work under Python 3 and R 4.1.2 with RStudio 2021.09.2+382.

## Citing
If you find HetIG-PreDiG useful for your research, please consider citing us:
```
@article{Jagodnik2022,
  title     = {HetIG-PreDiG: A Heterogeneous Integrated Graph Model for Predicting Human Disease Genes Based on Gene Expression.},
  author    = {Jagodnik, Kathleen M. and Shvili, Yael and Bartal, Alon},
  journal   = {},
  volume    = {},
  number    = {},
  pages     = {from page– to page},
  year      = {2022}
}
```

