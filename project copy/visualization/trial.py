import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter
import pickle
from visualization.data import Data

matplotlib.use('TkAgg')

with open('/Users/armanozcan/Desktop/project/visualization/model.pkl', 'rb') as f:
    model = pickle.load(f) # deserialize using load()
f.close()

with open('/Users/armanozcan/Desktop/project/visualization/ad.pkl', 'rb') as f:
    ad = pickle.load(f) # deserialize using load()
f.close()

print(model.get_hard_assignments().head())

labels,weights = model.get_soft_assignments()
SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')




# Normalize cells, log transform and compute highly variable genes
sc.pp.normalize_per_cell(SEACell_ad)
sc.pp.log1p(SEACell_ad)
sc.pp.highly_variable_genes(SEACell_ad, n_top_genes=800)

sns.distplot(SEACell_ad.to_df().var())

# Compute principal components -
# Here we use 50 components. This number may also be selected by examining variance explaint
sc.tl.pca(SEACell_ad, n_comps=10, use_highly_variable=True)

plt.plot(np.cumsum(SEACell_ad.uns['pca']['variance_ratio']))

sc.pp.neighbors(SEACell_ad, n_neighbors=5)

sc.tl.umap(SEACell_ad, n_components=2, min_dist=0.5, spread=1.0)

sc.pl.umap(SEACell_ad, edges=True, size=100)
sc.pl.scatter(SEACell_ad, basis='umap', size=100)



# TRIANGLES!!!

adjacency_matrix = (SEACell_ad.obsp["connectivities"].toarray()!=0).astype(int)

adjacency_list = []
for i in range(0,90):
    ind = np.arange(0, 90)
    adjacency_list.append(ind[adjacency_matrix[i] == 1])

count1 = 0
count2 = 0
threshold = 0.05
weights_by_cell = []
for i in range(model.A_.shape[1]):
    if(max(model.A_[:,i]) == 1):
        count2 += 1
    no = len(model.A_[:,i][model.A_[:,i]> threshold])
    if(no < 3):
        count1 +=1
        continue
    ind = np.argsort(model.A_[:,i])[::-1][:no]
    ind2 = []
    for j in range(len(ind)):
        ind2.append("{}{}".format("SEACell-", ind[j]))
    weights_by_cell.append(ind2)
print("Only one assignment: " + str(count2))
print("Only two assignments: " + str(count1 - count2))
print("Total removed assignments " + str(count1))

triangles = []
in_triangles = set()
data = [set() for x in range(90)]
for s in range(90):
    for t in adjacency_list[s]:
        if(s < t):
            for v in data[s].intersection(data[t]):
                triangles.append((v,s,t))
                in_triangles.add(v)
                in_triangles.add(s)
                in_triangles.add(t)
            data[t].add(s)
print("Number of triangles initially: " + str(len(triangles)))
print("Number of SEACELLS covered in these triangles: " + str(len(in_triangles)))
all = set(range(90))
print("Those SEACELLS were not covered: " + str(all.difference(in_triangles)))


nn_triangles = []
names = SEACell_ad.obs_names
for i in triangles:
    nn_triangles.append((names[i[0]], names[i[1]], names[i[2]]))
#print(nn_triangles)

counts = np.array([0] * len(nn_triangles))
confirmed_triangles = []
for index, triangle in enumerate(nn_triangles):
    for weights in weights_by_cell:
        if(len(set(triangle).difference(set(weights))) == 0):
            counts[index] += 1
    if(counts[index] > 3):
        confirmed_triangles.append(triangle)
#print(counts[counts > 0])
print("This is the count matrix:")
print(counts)
print("These are the triangles that were validated:")
print(confirmed_triangles)

print("Number of confirmed triangles is: " + str(len(confirmed_triangles)))
print("Number of removed triangles after comparing to weight assignments for each cell is: " + str(len(nn_triangles) - len(confirmed_triangles)))
all = set()
for i in range(90):
    all.add("{}{}".format("SEACell-", i))
confirmed_set = set()
for triangle in confirmed_triangles:
    confirmed_set.update(triangle)
print("Number of SEACells that do not exist in the remaining triangles is: " + str(len(all.difference(confirmed_set))))
print("It is this set:")
difference = all.difference(confirmed_set)
print(difference)


power = np.linalg.matrix_power(adjacency_matrix, 2)
print(np.trace(power))
power = np.linalg.matrix_power(adjacency_matrix, 3)
print(np.trace(power)/6)
power = np.linalg.matrix_power(adjacency_matrix, 4)
sum_degree_squared = np.sum(np.power(np.sum(adjacency_matrix, axis = 0), 2))
sum_degree = np.sum(np.sum(adjacency_matrix, axis = 0))
print((np.trace(power) - 2*sum_degree_squared + sum_degree)/8)


removed_triangles = set(nn_triangles).difference(set(confirmed_triangles))

data = Data(confirmed_triangles, removed_triangles, difference)
with open('data.pkl', 'wb') as f:  # open a text file
    pickle.dump(data, f) # serialize the list
f.close()


