import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter
import pickle

matplotlib.use('TkAgg')

# Load the data using
ad = sc.read('../data/cd34_multiome_rna.h5ad')

# Some plotting aesthetics

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100

sns.distplot(ad.to_df().sum(1))

# Plot cell-types for reference
sc.pl.scatter(ad, basis='umap', color='celltype', frameon=False)
# algorithms to construct convex hull

# Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
# This step should be performed after filtering
raw_ad = sc.AnnData(ad.X)
raw_ad.obs_names, raw_ad.var_names = ad.obs_names, ad.var_names
ad.raw = raw_ad

# Normalize cells, log transform and compute highly variable genes
sc.pp.normalize_per_cell(ad)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=1500)

# Compute principal components -
# Here we use 50 components. This number may also be selected by examining variance explained
sc.tl.pca(ad, n_comps=50, use_highly_variable=True)

## User defined parameters

## Core parameters
n_SEACells = 90
build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

model = SEACells.core.SEACells(ad,
                  build_kernel_on=build_kernel_on,
                  n_SEACells=n_SEACells,
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)

model.construct_kernel_matrix()
M = model.kernel_matrix

sns.clustermap(M.toarray()[:500,:500])

# Initialize archetypes
model.initialize_archetypes()

# Plot the initialization to ensure they are spread across phenotypic space
SEACells.plot.plot_initialization(ad, model)

model.fit(min_iter=10, max_iter=50)

# Check for convergence
model.plot_convergence()

plt.figure(figsize=(3,2))
sns.distplot((model.A_.T > 0.1).sum(axis=1), kde=False)
plt.title(f'Non-trivial (> 0.1) assignments per cell')
plt.xlabel('# Non-trivial SEACell Assignments')
plt.ylabel('# Cells')
plt.show()

plt.figure(figsize=(3,2))
b = np.partition(model.A_.T, -5)
sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
plt.title('Strength of top 5 strongest assignments')
plt.xlabel('$n^{th}$ strongest assignment')
plt.show()

with open('model.pkl', 'wb') as f:  # open a text file
    pickle.dump(model, f) # serialize the list
f.close()

with open('ad.pkl', 'wb') as f:  # open a text file
    pickle.dump(ad, f) # serialize the list
f.close()
