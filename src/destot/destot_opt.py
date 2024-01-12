'''
An Optimal-Transport based approach for quantifying the spatiotemporal growth of tissues.
'''
import numpy as np
import scipy
import matplotlib.pyplot as plt
import torch

def Cost(f, g, Grad, epsilon, device='cpu'):
    '''
    A matrix which is using for the broadcasted log-domain log-sum-exp trick-based updates.
    ------Parameters------
    f: torch.tensor (N1)
        First dual variable of semi-unbalanced Sinkhorn
    g: torch.tensor (N2)
        Second dual variable of semi-unbalanced Sinkhorn
    Grad: torch.tensor (N1 x N2)
        A collection of terms in our gradient for the update
    epsilon: float
        Entropic regularization for Sinkhorn
    device: 'str'
        Device tensors placed on
    '''
    return -( Grad - torch.outer(f, torch.ones(Grad.size(dim=1), device=device)) - torch.outer(torch.ones(Grad.size(dim=0), device=device), g) ) / epsilon

def normalize_M(M, C, N1, _ord='fro'):
    '''
    A function to return normalization constants for the gradient involving the merged feature-spatial matrix,
    which is normalized by a factor related to the interslice feature cost C.
    ------Parameters------
    M: torch.tensor (N x N)
        A symmetric matrix with positive entries (i.e. distance matrix, merged feature-spatial matrix)
    C: torch.tensor (N x M) or (M x N)
        Matrix of pairwise feature distances between slice 1 and 2
    N1: int
        Number of spots in first slice
    '''
    return (torch.linalg.norm(C, ord=_ord)**(1/2) / torch.linalg.norm(M, ord=_ord)) * (N1 )**(1/2)

def LogSinkhorn_iteration(C, D1, D2, C1, C2, Pi_0=None, alpha=0.2, beta=0.5, gamma=50, epsilon=1e-1, max_iter=100, balanced=False, device='cpu', dtype=torch.float64, override_EDM=False, override_FDM=False):
    '''
    Sinkhorn algorithm for the balanced/partially unbalanced case, with log-domain stabilization.
    ------Parameters------
    C: torch.tensor (N1 x N2)
        A matrix of pairwise feature distances between transcript vectors in slice 1 and slice 2 (interslice).
    D1: torch.tensor (N1 x N1)
        A matrix of pairwise Euclidean distances between points in slice 1.
    D2: torch.tensor (N2 x N2)
        A matrix of pairwise Euclidean distances between points in slice 2.
    C1: torch.tensor (N1 x N1)
        A matrix of pairwise feature distances between transcript vectors in slice 1 (intraslice).
    C2: torch.tensor (N2 x N2)
        A matrix of pairwise feature distances between transcript vectors in slice 2 (intraslice).
    Pi_0: torch.tensor (N1 x N2)
        An initialization for the alignment matrix Pi. Should respect marginals of semi-unbalanced or balanced.
    alpha: float
        A balance parameter between the interslice feature term of the objective and the merged feature-spatial term.
    beta: float
        A balance parameter between the GW (quartet) term and the triplet term of the merged feature-spatial term.
    gamma: float
        A hyperparameter controlling the strength of the KL-divergence term in the objective.
    epsilon: float
        A hyperparameter controlling the strength of the entropic regularization in the Sinkhorn algorithm.
    max_iter: int
        The maximal number of iterations DeST-OT is run.
    balanced: bool
        Boolean for whether to default to a balanced OT or to use DeST-OT's semi-unbalanced routine. Default set to False.
    device: str
        Device that torch tensors are placed on. Using GPU/'cuda' is much faster.
    dtype: torch.type
        The default datatype that the alignment and other tensors are in. Ideally torch.float64 or torch.float32
    override_EDM: bool
        Whether to override the merged feature-spatial matrix with a standard Euclidean distance matrix.
    override_FDM:
        Whether to override the merged feature-spatial matrix with an intraslice feature distance matrix.
    '''
    
    N1, N2 = C.size(dim=0), C.size(dim=1)
    
    M1 = D1 * C1
    M2 = D2 * C2
    
    if override_EDM:
        # Override fused matrix with standard EDM
        M1 = D1
        M2 = D2
    elif override_FDM:
        # Override fused matrix with feature distance matrix
        M1 = C1
        M2 = C2
    
    # Normalizing constants
    p1, p2 = normalize_M(M1, C, N1, _ord='fro'), \
                            normalize_M(M2, C, N1, _ord='fro')
    r2, r1 = torch.linalg.norm(C, ord='fro')/torch.linalg.norm(M2**2, ord='fro') * (N1**1/2), \
                            torch.linalg.norm(C, ord='fro')/torch.linalg.norm(M1**2, ord='fro') * (N1 / N2**(1/2))
    
    k = 0
    stationarity_gap = torch.inf
    
    # These are the same for this synthetic data-- change dimensions when we get there
    one_N1 = torch.ones((N1), device=device, dtype=dtype)
    one_N2 = torch.ones((N2), device=device, dtype=dtype)
    
    g1 = one_N1 / N1
    
    if balanced:
        g2 = one_N2 / N2
    else:
        g2 = one_N2 / N1
    
    log_g1 = torch.log(g1)
    log_g2 = torch.log(g2)
    
    if Pi_0 is None:
        # Will converge to partially balanced marginal
        Pi_0 = torch.outer(g1, g2).to(device)
    
    Pi_k = Pi_0
    
    f_k = torch.zeros((N1), device=device)
    g_k = torch.zeros((N2), device=device)
    
    errs = []
    
    # unbalanced coefficient
    ubc = gamma/(gamma + epsilon)
    
    grad = torch.inf
    while k < max_iter:
        
        if k % 5 == 0:
            print(f'Iteration: {k}')
        
        # Offering the user the option to choose their energy-regularization
        dist_orig = r1*(M1**2 @ Pi_k) + r2*(Pi_k @ M2**2)
        
        if balanced:
            dist_GW = -2 * ((p1*M1) @ Pi_k @ (p2*M2)) 
        else:
            dist_GW = -2*((p1*M1) @ Pi_k @ (p2*M2).T) + (p1*M1)**2 @ Pi_k @ torch.ones((N2, N2), device=device, dtype=dtype)
        
        grad = (1-alpha)*C + alpha*(beta*dist_orig + (1-beta)*dist_GW)
        
        if balanced:
            f_k = f_k + epsilon*(log_g1 - torch.logsumexp(Cost(f_k, g_k, grad, epsilon, device=device), axis=1))
            g_k = g_k + epsilon*(log_g2 - torch.logsumexp(Cost(f_k, g_k, grad, epsilon, device=device), axis=0))
        elif balanced is False:
            # Partially unbalanced coefficient used on one of the dual variables
            f_k = ubc*(f_k + epsilon*(log_g1 - torch.logsumexp(Cost(f_k, g_k, grad, epsilon, device=device), axis=1)) )
            g_k = g_k + epsilon*(log_g2 - torch.logsumexp(Cost(f_k, g_k, grad, epsilon, device=device), axis=0))
        
        Pi = torch.exp(Cost(f_k, g_k, grad, epsilon, device=device))
        # Checking whether we're approaching a fixed point
        stationarity_gap = torch.linalg.norm(Pi_k - Pi, ord='fro')
        Pi_k = Pi
        
        k+=1
        errs.append(stationarity_gap)
    
    xi = (Pi_k @ one_N2 - g1)
    
    return xi, Pi_k, errs
