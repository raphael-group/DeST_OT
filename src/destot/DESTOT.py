from .destot_opt import LogSinkhorn_iteration

import torch
import anndata as ad
import scanpy as sc
import numpy as np
from scipy.spatial import distance
from pynvml import *


def intersect(lst1, lst2): 
    """
    param: lst1 - list
    param: lst2 - list
    
    return: list of common elements
    """
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3 

"""
This function to select a GPU is from the scSLAT package.
"""
def get_free_gpu() -> int:
    """
    Get index of GPU with least memory usage
    
    Ref
    ----------
    https://stackoverflow.com/questions/58216000/get-total-amount-of-free-gpu-memory-and-available-using-pytorch
    """
    nvmlInit()
    index = 0
    max = 0
    for i in range(torch.cuda.device_count()):
        h = nvmlDeviceGetHandleByIndex(i)
        info = nvmlDeviceGetMemoryInfo(h)
        index = i if info.free > max else index
        max = info.free if info.free > max else max
        
    # seed = np.random.randint(1000)
    # os.system(f'nvidia-smi -q -d Memory |grep Used > gpu-{str(seed)}.tmp')
    # memory_available = [int(x.split()[2]) for x in open('gpu.tmp', 'r').readlines()]
    # os.system(f'rm gpu-{str(seed)}.tmp')
    # print(memory_available)
    return index

def xi_to_growth_rate(xi, t1=0, t2=1, normalize_xi=True):
    '''
    Returns a differential growth rate given the growth vector xi.
    
    Parameters:
    xi: numpy array (N1)
        Growth vector quantifying the raw mass-flux
    t1: float
        First observation timepoint
    t2: float
        Second observation timepoint
    normalize_xi: bool
        True if xi normalized to number of cells units (this is the default output of align, using the same flag "normalize_xi"),
        False if xi directly computed from Pi without renormalization.
    '''
    N1 = xi.shape[0]
    if normalize_xi is False:
        Js = np.log(N1*xi + 1) / (t2 - t1)
    else:
        Js = np.log(xi + 1) / (t2 - t1)
    # Returning a proper growth-rate given mass-flux xi
    return Js

def align(slice_t1, slice_t2, alpha=0.2, gamma=50, \
          epsilon=1e-1, beta=0.5, max_iter=100, balanced=False, \
          use_gpu=True, normalize_xi=True, \
          check_convergence=False, spatial=True, \
         feature_key = 'X_pca', spatial_key = 'spatial'):
    """
    Run DeST-OT

    Parameters:
    slice_t1: AnnData object
        The AnnData object of the first slice, with .obsm['spatial'] field storing the spatial coordinates
    slice_t2: AnnData object
        The AnnData object of the second slice, with .obsm['spatial'] field storing the spatial coordinates
    alpha: float
        A balance parameter between the interslice feature term of the objective and the merged feature-spatial term. Default is 0.2.
    gamma: float
        A hyperparameter controlling the strength of the KL-divergence term in the objective. Default is 50.
    epsilon: float
        A hyperparameter controlling the strength of the entropic regularization in the Sinkhorn algorithm. Default is 0.1.
    max_iter: int
        The maximal number of iterations DeST-OT is run. Default is 100.
    balanced: bool
        Boolean for whether to default to a balanced OT or to use DeST-OT's semi-unbalanced routine. Default set to False.
    use_gpu: bool
        Boolean for whether to use GPU. Default is True.
    normalize_xi: bool
        Boolean for whether to normalize the growth vector xi to the unit of number of spots. If True, each entry in xi will be in the unit of number of spots, e.g. xi_i = 1 means spots i will grow into two spots in the next timepoint. Default is False.
        Note: set to False when use the output xi to compute the growth distortion metric.
    check_convergence: bool
        Boolean for whether to return the stationarity gap with respect to Pi for each iteration.

    Returns:
    Pi: numpy array of shape (N1 x N2)
        The alignment matrix
    xi: numpy array of shape (N1)
        The growth vector
    errs: list of size (max_iter)
        If check_convergence is True, return a list of stationarity gaps for each iterations. Used for checking the convergence of Pi.
    feature_key: str, optional (Default = 'X_pca')
        Key for feature data used for alignment (e.g. transcriptomics PCA vector)
    spatial_key: str, optional (Default = 'spatial')
        Key for spatial data used for alignment (e.g. spatial transcriptomics physical coordinates)
    """
    
    # subset for common genes
    common_genes = intersect(slice_t1.var.index, slice_t2.var.index)
    slice_t1 = slice_t1[:, common_genes]
    slice_t2 = slice_t2[:, common_genes]

    # Calculate gene distances
    joint_adata = ad.concat([slice_t1, slice_t2])
    sc.pp.normalize_total(joint_adata, inplace=True)
    sc.pp.log1p(joint_adata)
    sc.pp.pca(joint_adata, 30)
    
    joint_datamatrix = joint_adata.obsm[feature_key]
    
    X = joint_datamatrix[:slice_t1.shape[0], :]
    Y = joint_datamatrix[slice_t1.shape[0]:, :]
    
    C = distance.cdist(X, Y)
    C1 = distance.cdist(X, X)
    C2 = distance.cdist(Y, Y)
    
    if spatial:
        # Calculate spatial distances
        D1 = distance.cdist(slice_t1.obsm[spatial_key], slice_t1.obsm[spatial_key])
        D2 = distance.cdist(slice_t2.obsm[spatial_key], slice_t2.obsm[spatial_key])
    else:
        # Else, only use expression information
        D1, D2 = np.ones(C1.shape), np.ones(C2.shape)
    
    # Pytorch
    if use_gpu:
        gpu_index = get_free_gpu()
        device = torch.device(f'cuda:{gpu_index}' if torch.cuda.is_available() else 'cpu')
    else:
        device='cpu'
        
    C = torch.from_numpy(C).to(device)
    C1 = torch.from_numpy(C1).to(device)
    C2 = torch.from_numpy(C2).to(device)
    D1 = torch.from_numpy(D1).to(device)
    D2 = torch.from_numpy(D2).to(device)

    # Run DeST-OT
    xi, Pi, errs = LogSinkhorn_iteration(C, D1, D2, \
                                         C1, C2, alpha=alpha, \
                                         gamma=gamma, epsilon=epsilon, \
                                         beta=beta, max_iter=max_iter, \
                                         balanced=balanced, device=device)
    
    Pi = Pi.cpu().detach().numpy()
    xi = xi.cpu().detach().numpy()
    
    if normalize_xi:
        xi = slice_t1.shape[0] * xi
        
    if check_convergence:
        return Pi, xi, errs
        
    return Pi, xi