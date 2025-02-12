import numpy as np
import scanpy as sc
from scipy.spatial import distance
import pandas as pd


def partial_procrustes_analysis(X, Y, pi):
    m = np.sum(pi)
    Z = (X - pi.sum(axis=1).dot(X) * (1.0 / m)).T
    W = (Y - pi.sum(axis=0).dot(Y) * (1.0 / m)).T
    H = W.dot(pi.T.dot(Z.T))
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T.dot(U.T)
    W = R.dot(W)
    return Z.T, W.T


def partial_stack_slices_pairwise(slices, pis):
    assert len(slices) == len(pis) + 1, "'slices' should have length one more than 'pis'. Please double check."
    assert len(slices) > 1, "You should have at least 2 layers."

    new_coor = []
    S1, S2 = partial_procrustes_analysis(slices[0].obsm['spatial'], slices[1].obsm['spatial'], pis[0])
    new_coor.append(S1)
    new_coor.append(S2)
    for i in range(1, len(slices) - 1):
        x, y = partial_procrustes_analysis(new_coor[i], slices[i + 1].obsm['spatial'], pis[i])
        shift = new_coor[i][0,:] - x[0,:]
        y = y + shift
        new_coor.append(y)
    new_slices = []
    for i in range(len(slices)):
        s = slices[i].copy()
        s.obsm['spatial'] = new_coor[i]
        new_slices.append(s)
    return new_slices


def transition_mat(Pi, celltypes1, celltypes2, column_norm = True):
    """
    Compute the cell type transition matrix from an alignment matrix Pi (Eq. 10 of the paper)

    Parameters:
    Pi: numpy array of shape (N1, N2)
        The alignment matrix
    celltype1: numpy array of shape N1
        An array of celltypes for each spot on slice t1
    celltype2: numpy array of shape N2
        An array of celltypes for each spot on slice t2
    column_norm: bool, optional (default=True)
        Whether to return a column-normalized or a row-normalized transition matrix.
        
    Return:
    The cell type transition matrix, a numpy array of shape (N_CT, N_CT), where N_CT is the total number of cell types
    """
    celltypes_all = np.unique(np.concatenate((celltypes1,celltypes2), axis=0))
    N_CT = celltypes_all.shape[0]
    
    T = np.zeros((N_CT, N_CT))
    for i in range(N_CT):
        for j in range(N_CT):
            ct_i = celltypes_all[i]
            ct_j = celltypes_all[j]
            mask = np.outer((celltypes1 == ct_i), (celltypes2 == ct_j))
            T[i,j] = np.sum(Pi[mask])

    if column_norm:
        # Column-normalized
        
        col_sums = T.sum(axis=0)
        non_zero_cols = col_sums != 0
        T[:, non_zero_cols] /= col_sums[non_zero_cols]
    else:
        # Row-normalized
        
        row_sums = T.sum(axis=1, keepdims=True)
        T /= row_sums
    
    return T


def growth_distortion_metric_helper(xi, celltypes1, celltypes2, T=None):
    '''
    Implementation of the growth distortion metric given the inferred growth vector xi

    Parameters:
    xi: numpy array of shape (N1, N2)
        The growth vector of a spatiotemporal alignment
    celltype1: numpy array of shape N1
        An array of celltypes for each spot on slice t1
    celltype2: numpy array of shape N2
        An array of celltypes for each spot on slice t2
    T: numpy array of shape (N1, N2)
        A cell type transition matrix to compute the growth distortion metric under.

    Returns:
    the growth distortion metric
    '''
    N1 = xi.shape[0]
    if T is not None:
        celltypes = np.unique(np.concatenate((celltypes1,celltypes2), axis=0))
    else:
        celltypes = np.intersect1d(celltypes1, celltypes2)
    
    distortion_measure = 0
    
    N_C = len(celltypes)
    m_t = np.zeros(N_C)
    m_tm1 = np.zeros(N_C)
    
    for p, celltype in np.ndenumerate(celltypes):
        m_t[p] = np.sum(celltypes2 == celltype)
        m_tm1[p] = np.sum(celltypes1 == celltype)
    
    if T is not None:
        # Transition-matrix adjusted assignment of mass-fluxes
        m_t = T @ m_t
    
    # Quantifying the accuracy of matching true baseline growth rates
    growth_rates = np.zeros(celltypes.shape[0])
    for p, celltype in np.ndenumerate(celltypes):
        dmP = (m_t[p] - m_tm1[p])/N1
        if m_tm1[p] != 0:
            growth_rates[p] = gamma_tp = (1/N1)*((m_t[p] - m_tm1[p])/m_tm1[p])
            xi_p = xi[celltypes1 == celltype]
            # Sum of squares error of individual cell growth rates,
            # relative to true cell type specific growth rate.
            distortion_measure += np.sum( (xi_p - gamma_tp * np.ones(xi_p.shape[0]) )**2 )
        else:
            # Pass if undefined.
            pass
            
    # print(f'Distortion metric value: {N1*distortion_measure}')
    return N1 * distortion_measure


def growth_distortion_metric(slice_t1, slice_t2, Pi, xi=None, annotation_key="annotation", option="infer_transition"):
    """
    Compute the growth distortion metric given two slices, an alignment matrix Pi, and a growth vector xi

    Parameters:
    slice_t1: AnnData object
        The AnnData object of the first slice, with .obs['annotation'] field storing the cell type of each spot.
    slice_t2: AnnData object
        The AnnData object of slice t2, with .obs['annotation'] field storing the cell type of each spot.
    Pi: numpy array of shape (N1, N2)
        Alignment matrix between slice t1 and slice t2
    xi: numpy array of shape (N1) or NoneType
        The growth vector from the alignment. If not input, then recomputed from Pi.
    annotation_key: String
        The key for the cell-type annotations in the AnnData object
    option: String, one of the following two options
        "no_transition": Assumes no cell type transition, and the growth distortion metric is calculcated based on the intersection of cell types in the two slices
        "infer_transition": Assumes cell type transition during development, and the growth distortion metric is calculated based on the cell type transition matrix of all cell types that minimizes the growth distortion metric for the given Pi (Section 2.3.1 of the paper)

    Returns:
    The growth distortion metric of the given Pi and xi
    """
    if xi is None:
        one_N2 = np.ones(Pi.shape[1])
        g1 = np.ones(Pi.shape[0])/Pi.shape[0]
        xi = (Pi @ one_N2 - g1)
    else:
        pass
    l1, l2 = slice_t1.obs[annotation_key].tolist(), slice_t2.obs[annotation_key].tolist()
    l_merged = l1 + l2
    categorical_labels, celltypes = pd.factorize(l_merged)
    celltypes1 = categorical_labels[:len(l1)]
    celltypes2 = categorical_labels[len(l1):]

    if option == 'no_transition':
        T = None
    elif option == 'infer_transition':
        T = transition_mat(Pi, celltypes1, celltypes2)
    else:
        raise Exception('Invalid Option')
    
    return growth_distortion_metric_helper(xi, celltypes1, celltypes2, T)


def migration_metric(slice_t1, slice_t2, Pi):
    """
    Compute the migration metric given two slices and an aligment matrix Pi

    Parameters:
    slice_t1: AnnData object
        An AnnData object of slice t1, with .obsm['spatial'] field storing the spatial coordinates.
    slice_t2: AnnData object
        An AnnData object of slice t2, with .obsm['spatial'] field storing the spatial coordinates.
    Pi: numpy array of shape (N1 x N2)
        Alignment matrix between slice t1 and slice t2

    Returns:
    The migration metric of Pi
    """
    slice_t1_newcoor_adata, slice_t2_newcoor_adata = partial_stack_slices_pairwise([slice_t1, slice_t2], [Pi])
    
    slice_t1_newcoor = slice_t1_newcoor_adata.obsm['spatial']
    slice_t2_newcoor = slice_t2_newcoor_adata.obsm['spatial']
    
    V = distance.cdist(slice_t2_newcoor, slice_t1_newcoor)
    
    # Handling zero entries in Pi @ 1_n (i.e. spot not mapped at all to next timepoint)
    ma = (Pi @ np.ones(Pi.shape[1])) > 0

    Pi_hat = Pi[ma, :] / (Pi @ np.ones(Pi.shape[1]))[ma][:, None]
    # Pi_hat = pi / (pi @ np.ones(pi.shape[1]))[:, None]
    V = V[:,ma]

    migs = np.einsum('ij,ji->i', Pi_hat, V)
    
    avg_weighted_migration_distance = np.mean(migs)
    return avg_weighted_migration_distance

