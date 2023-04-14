import scanpy as sc
import numpy as np
import SEACells
from visualization.data import Data
import networkx as nx
import plotly.graph_objects as go
import plotly

def adjacency_matrix_to_graph(adj_matrix):
    """
    Convert an adjacency matrix to a networkx graph.
    """
    graph = nx.Graph()
    graph.add_nodes_from(range(adj_matrix.shape[0]))
    for i in range(adj_matrix.shape[0]):
        for j in range(i+1, adj_matrix.shape[0]):
            if adj_matrix[i,j] == 1:
                graph.add_edge(i,j)
    return graph
def annotate_nodes(labels, graph):
    """
    Annotate nodes in a graph with a list of labels.
    """
    for i, label in enumerate(labels):
        graph.nodes[i]['label'] = label
    return graph
def add_coordinates(coords, graph):
    """
    Add coordinates to a graph.
    """
    for i, coord in enumerate(coords):
        graph.nodes[i]['x'] = coord[0]
        graph.nodes[i]['y'] = coord[1]
    return graph



def plot(adjacency_matrix, seacell_labels, UMAP_coords):
    G = adjacency_matrix_to_graph(adjacency_matrix)
    G = annotate_nodes(seacell_labels, G)
    G = add_coordinates(UMAP_coords, G)

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['x'], G.nodes[edge[0]]['y']
        x1, y1 = G.nodes[edge[1]]['x'], G.nodes[edge[1]]['y']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['x'], G.nodes[node]['y']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: ' + str(len(adjacencies[1])))

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='<br>Network graph made with Python',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        annotations=[dict(
                            text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002)],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    return plotly.io.to_html(fig)

def triangles(model, ad):

    labels, weights = model.get_soft_assignments()
    SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')

    # Normalize cells, log transform and compute highly variable genes
    sc.pp.normalize_per_cell(SEACell_ad)
    sc.pp.log1p(SEACell_ad)
    sc.pp.highly_variable_genes(SEACell_ad, n_top_genes=800)

    # Compute principal components -
    # Here we use 50 components. This number may also be selected by examining variance explaint
    sc.tl.pca(SEACell_ad, n_comps=10, use_highly_variable=True)

    sc.pp.neighbors(SEACell_ad, n_neighbors=5)

    sc.tl.umap(SEACell_ad, n_components=2, min_dist=0.5, spread=1.0)

    # TRIANGLES!!!

    adjacency_matrix = (SEACell_ad.obsp["connectivities"].toarray() != 0).astype(int)

    adjacency_list = []
    for i in range(0, 90):
        ind = np.arange(0, 90)
        adjacency_list.append(ind[adjacency_matrix[i] == 1])

    count1 = 0
    count2 = 0
    threshold = 0.05
    weights_by_cell = []
    for i in range(model.A_.shape[1]):
        if (max(model.A_[:, i]) == 1):
            count2 += 1
        no = len(model.A_[:, i][model.A_[:, i] > threshold])
        if (no < 3):
            count1 += 1
            continue
        ind = np.argsort(model.A_[:, i])[::-1][:no]
        ind2 = []
        for j in range(len(ind)):
            ind2.append("{}{}".format("SEACell-", ind[j]))
        weights_by_cell.append(ind2)

    triangles = []
    in_triangles = set()
    data = [set() for x in range(90)]
    for s in range(90):
        for t in adjacency_list[s]:
            if (s < t):
                for v in data[s].intersection(data[t]):
                    triangles.append((v, s, t))
                    in_triangles.add(v)
                    in_triangles.add(s)
                    in_triangles.add(t)
                data[t].add(s)

    nn_triangles = []
    names = SEACell_ad.obs_names
    for i in triangles:
        nn_triangles.append((names[i[0]], names[i[1]], names[i[2]]))
    # print(nn_triangles)

    counts = np.array([0] * len(nn_triangles))
    confirmed_triangles = []
    for index, triangle in enumerate(nn_triangles):
        for weights in weights_by_cell:
            if (len(set(triangle).difference(set(weights))) == 0):
                counts[index] += 1
        if (counts[index] > 3):
            confirmed_triangles.append(triangle)

    all = set()
    for i in range(90):
        all.add("{}{}".format("SEACell-", i))
    confirmed_set = set()
    for triangle in confirmed_triangles:
        confirmed_set.update(triangle)
    difference = all.difference(confirmed_set)

    removed_triangles = set(nn_triangles).difference(set(confirmed_triangles))

    data = Data(confirmed_triangles, removed_triangles, difference, adjacency_matrix, SEACell_ad.obs_names, SEACell_ad.obsm['X_umap'])

    return data
# not using now
def seacells(filename):
    ad = sc.read(filename)
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
    build_kernel_on = 'X_pca'  # key in ad.obsm to use for computing metacells
    # This would be replaced by 'X_svd' for ATAC data

    ## Additional parameters
    n_waypoint_eigs = 10  # Number of eigenvalues to consider when initializing metacells

    model = SEACells.core.SEACells(ad,
                                   build_kernel_on=build_kernel_on,
                                   n_SEACells=n_SEACells,
                                   n_waypoint_eigs=n_waypoint_eigs,
                                   convergence_epsilon=1e-5)

    model.construct_kernel_matrix()
    M = model.kernel_matrix

    # Initialize archetypes
    model.initialize_archetypes()

    model.fit(min_iter=10, max_iter=50)

    '''
    with open('/Users/armanozcan/Desktop/project/visualization/pickles/model.pkl', 'wb') as f:  # open a text file
        pickle.dump(model, f)  # serialize the list
    f.close()

    with open('/Users/armanozcan/Desktop/project/visualization/pickles/ad.pkl', 'wb') as f:  # open a text file
        pickle.dump(ad, f)  # serialize the list
    f.close()
    '''