import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from src.destot.metrics import transition_mat

def plot_slice_value(slice, value_vec, vmax=None, vmin=None):
    """
    Parameters: 
    slice: AnnData object of the slice
    value_vec: growth vector

    Returns:
    Plots the slice with each spot colored according to its value in value_vec
    """
    plt.figure()
    spatial = slice.obsm['spatial']
    sc = plt.scatter(spatial[:, 0], spatial[:, 1], c=value_vec, cmap='RdYlGn', s=50, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=20)
    plt.gca().invert_yaxis()
    plt.axis('off')

    fig = plt.gcf()
    fig_size = fig.get_size_inches()
    new_width = 20.0
    new_height = new_width * (fig_size[1] / fig_size[0])
    fig.set_size_inches(new_width, new_height)
    plt.show()
    return

def get_transition(Pi, slice1, slice2, annotation_key = 'Annotation', plotting=True, column = True):
    """
    Computes the transition matrix between cell-types with a doubly-stochastic alignment matrix.
    
    Parameters
    ----------
    Pi: np.ndarray ( n1 x n2 )
        A doubly-stochastic alignment matrix representing a probabilistic alignment from time 1 (t1) to time 2 (t2)
    slice1: AnnData
        AnnData Object at time 1 containing cell-type states
    slice2: AnnData Object
        AnnData Object at time 2  containing cell-type states
    annotation_key: str, optional (default='Annotation')
        Key in 'obs' for accessing cell-types in AnnData objects.
    plotting: bool, optional (default=True)
        Whether to visualize the transition matrix with a heatmap.
    column: bool, optional (default=True)
        Whether to return a column-normalized or a row-normalized transition matrix.
    Returns
    -------
    T: np.ndarray
        Stochastic transition matrix capturing changes in cell-type proportions from time 1 to time 2.
    labelX: list (str)
        List of cell-type labels ordered in accordance with T at time 1
    labelY: list (str)
        List of cell-type labels ordered in accordance with T at time 2
    """
    # Extract cell-types from AnnData objects
    l1, l2 = slice1.obs[annotation_key].tolist(), slice2.obs[annotation_key].tolist()
    l_merged = l1 + l2
    
    categorical_labels, celltypes = pd.factorize(l_merged) # Getting numeric indices for the labels
    mapping = dict(enumerate(celltypes))

    # Split by time
    celltypes1 = categorical_labels[:len(l1)]
    celltypes2 = categorical_labels[len(l1):]

    if column:
        # Computes transitions between cell-types assuming global labels assuming column-normalization
        T = transition_mat(Pi, celltypes1, celltypes2, column_norm = True)
    else:
        # Row-normalized
        T = transition_mat(Pi, celltypes1, celltypes2, column_norm = False)
    
    unique_celltypes = sorted(set(categorical_labels), key=list(categorical_labels).index)
    ct_labels = [mapping[ct] for ct in unique_celltypes]

    # Filtering on labels specific to each timepoint, after global transition across all labels computed
    labelX, idxX = [], []
    labelY, idxY = [], []
    
    for ct in unique_celltypes:
        if mapping[ct] in l1:
            # CT present in time 1
            labelX.append(mapping[ct])
            idxX.append(ct)
        if mapping[ct] in l2:
            # CT present in time 2
            labelY.append(mapping[ct])
            idxY.append(ct)

    # Masking out unrepresented types at each timepoint
    T = T[np.ix_(idxX, idxY)]

    # If true, plot the transitions
    if plotting:
        
        plt.rcParams['figure.dpi'] = 300
        fig, ax = plt.subplots(1,1)
        img = ax.imshow(T, vmin=0, vmax=1, cmap='Greens')

        # Setting the axis ticks & labels
        ax.set_xticks(range(len(labelY)))
        ax.set_yticks(range(len(labelX)))
        
        ax.set_xticklabels(labelY, rotation=90)
        ax.set_yticklabels(labelX)
        
        ax.tick_params(labelbottom=False,labeltop=True)
        ax.set_ylabel(r"Cell type at time $(t-1)$")
        ax.set_xlabel(r"Cell type at time $t$")
    
    return T, labelX, labelY


