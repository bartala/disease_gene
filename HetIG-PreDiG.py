import numpy as np
import networkx as nx
import pandas as pd
import os
from numpy import asarray
from numpy import savetxt
from numpy import loadtxt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
import pickle
from node2vec import Node2Vec

PTH = '/path/to/data/'


# load data of the original network after edge removeal
## type=1: gene-disease
## type=2: gene-gene
df = pd.read_csv(os.path.join(PTH,'G_tag.csv'))

# create an undirected networkx graph
G_karate = nx.from_pandas_edgelist(df, 'from', 'to')
G_karate = nx.to_undirected(G_karate)

n = G_karate.number_of_nodes()
m = G_karate.number_of_edges()
print("Number of nodes :", str(n))
print("Number of edges :", str(m))
print("Number of connected components :", str(nx.number_connected_components(G_karate)))

"""## Learn node embeddings"""

# compute transition probabilities and generate walks
node2vec = Node2Vec(G_karate, dimensions=64, walk_length=5, num_walks=10, p=1, q=1, workers=4)

# Learn nodes embeddings
model = node2vec.fit(window=10, min_count=1, batch_words=4)

# Save embeddings for later use
model.wv.save_word2vec_format(os.path.join(PTH,'model/W2V_EMBEDDINGS_m4'))
# Save model for later use
model.save(os.path.join(PTH,'model/W2V_MODEL_m4'))

"""## Load the node2vec model"""

# load the node2vec saved model
from gensim.models import Word2Vec
 
model_name = 'W2V_MODEL_m4'
 
model = Word2Vec.load(os.path.join(PTH,'model',model_name))

"""
Load the positive/ negative examples and attach node embeddings for each node.
 `type` = 1 is positive example of gene-disease edge from the original graph
 `type` = 0 is negative example of edge that is not in the original graph
"""

training_testing_edges = pd.read_csv(os.path.join(PTH,'pos_neg_edges.csv'))
training_testing_edges

# get unique node names of training pos/neg edges
nodes_to = training_testing_edges['to'].tolist()
nodes_from = training_testing_edges['from'].tolist()
nodes = set(nodes_to + nodes_from)
len(nodes)

# get embeddings for each node and save to file
def get_embeddings(nodes):
  nodes_v = {}
  for node in nodes:
    try:
      nodes_v[str(node)] = model.wv.get_vector(node)
    except:
      pass
  return(nodes_v)

nodes_v = get_embeddings(nodes)

"""
For each row in `training_testing_edges` (pair of nodes --> i.e. edges) sum teir embeddings into a single vector `x`
* Recall in `training_testing_edges` each 'from' node is a gene and 'to' node is a disease 
"""

# load the gene-disease prioritization score matrix ( rows = genes, columns = diseases )
gene_dis_mat = pd.read_csv(os.path.join(PTH,"disease_prediction.tsv.zip"),compression='zip',sep="\t")

# normalize the columns of the GDPS matrix between [0,1]
x = gene_dis_mat.values
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
df = pd.DataFrame(x_scaled)
df.index = gene_dis_mat.index
df.columns = gene_dis_mat.columns

# get node embeddings
x = []
missing = []
for i in range(0,len(training_testing_edges)):
  to = training_testing_edges.iloc[i]['to'] # this is the disease
  frm = training_testing_edges.iloc[i]['from'] # this is the gene
  try:
    node_frm = nodes_v[frm]
    node_to = nodes_v[to]
    x.append(node_to + node_frm)
  except:
    missing.append(i)
 
x = pd.DataFrame(x)
print(x.shape)

modDfObj = training_testing_edges.drop(missing)
modDfObj.shape
modDfObj.to_csv(os.path.join(PTH,"modDfObj.csv")) 
x = pd.concat([modDfObj.reset_index(drop=True), x], axis=1)

attention = []
mis_g = []
for i in range(0,len(x)):
  try:
    to = x['to'].iloc[i]
    frm = x['from'].iloc[i]
    attention.append(df[to][frm]) # append gene prioritization score to embedding vector
  except:
    attention.append(0)
    
x['attention'] = attention

x.to_csv(os.path.join(PTH,"attention.csv"),index=False)

"""## Train classifiers"""

x = pd.read_csv(os.path.join(PTH,"attention.csv"))
x = x.drop(x.columns[[0, 1]], axis = 1)

def report(y_pred, y_test):
  from sklearn.metrics import confusion_matrix
  from sklearn.metrics import roc_auc_score
  
  cm = confusion_matrix(y_test, y_pred)
  print('Confusion Matrix [Actual (rows) X Predicted(columns)]: \n', cm)
  print(classification_report(y_test, y_pred))
  sensitivity1 = cm[1,1]/(cm[1,1]+cm[1,0])
  print('Sensitivity : ', sensitivity1 )
  specificity1 = cm[0,0]/(cm[0,0]+cm[0,1])
  print('Specificity : ', specificity1)
  print("AUC:",roc_auc_score(y_test, y_pred))
  print("\n")

# for CV
from sklearn.model_selection import cross_val_predict

def evaluate(X,Y,model,kfold):
  y_pred = cross_val_predict(model, X, Y, cv=10)
  print(classification_report(Y, y_pred))

# Split the file to 70% training/ 30% testing sets
modDfObj = pd.read_csv(os.path.join(PTH,"modDfObj.csv"))

y = modDfObj['type']
X_train, X_test, y_train, y_test = train_test_split(x.iloc[:,1:], y, test_size=0.3, shuffle=True)

print(X_train.shape, y_train.shape)
print(X_test.shape, y_test.shape)

# Train a logistic regression model No GDPS
logreg_nogdps = LogisticRegression(max_iter=300).fit(X_train.iloc[:,0:-1], y_train)
y_pred = logreg_nogdps.predict(X_test.iloc[:,0:-1])
report(y_pred, y_test)

# cross validation
#evaluate(X_train,y_train,logreg_nogdps,10)

# Train a logistic regression model
logreg = LogisticRegression(max_iter=300).fit(X_train, y_train)

y_pred = logreg.predict(X_test)
report(y_pred, y_test)

# save the model to disk
filename = os.path.join(PTH,'model/LR.sav')
pickle.dump(logreg, open(filename, 'wb'))

"""auc plot"""

# roc curve and auc
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
import numpy as np

# predict probabilities
lr_probs3 = model2.predict(X_test) #DFNN

lr_probs = logreg.predict_proba(X_test.iloc[:,0:]) # LR with GDPS
lr_probs = lr_probs[:, 1]


lr_probs2 = logreg_nogdps.predict_proba(X_test.iloc[:,0:-1]) # LR no GDPS
lr_probs2 = lr_probs2[:, 1]

# calculate scores
lr_auc = roc_auc_score(y_test, lr_probs)
lr_auc2 = roc_auc_score(y_test, lr_probs2)

# summarize scores
print('LR: ROC AUC=%.3f' % (lr_auc))
print('LR without GDPS: ROC AUC=%.3f' % (lr_auc2))

# calculate roc curves
lr_fpr, lr_tpr, _ = roc_curve(y_test, lr_probs)
lr_fpr2, lr_tpr2, _ = roc_curve(y_test, lr_probs2)

auc = round(metrics.roc_auc_score(y_test, lr_probs),2)
auc2 = round(metrics.roc_auc_score(y_test, lr_probs2),2)

# plot the roc curve for the model
fig1 = pyplot.gcf()
pyplot.plot(lr_fpr, lr_tpr, marker='', label='LR, AUC='+str(0.954), color='blue' )
pyplot.plot(lr_fpr2, lr_tpr2, marker='', label='LR without GDPS, AUC='+str(0.93), color='red')

# axis labels
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
 
pyplot.plot([0, 1], [0, 1], color='black', lw=1, linestyle='--')

# show the legend
pyplot.legend(loc='best')
# show the plot
pyplot.show()


fig1.savefig(os.path.join(PTH,"savefig.pdf"))

# Zoom in view of the upper left corner.

pyplot.figure(2)

fig2 = pyplot.gcf()


pyplot.xlim(0, 0.2)
pyplot.ylim(0.75, 1)
pyplot.plot([0, 1], [0, 1], 'k--')

pyplot.plot(lr_fpr, lr_tpr, marker='',  color='blue' )
pyplot.plot(lr_fpr2, lr_tpr2, marker='', color='red')

pyplot.xlabel('False positive rate')
pyplot.ylabel('True positive rate')
pyplot.show()

fig2.savefig(os.path.join(PTH,"savefig3.pdf"))
