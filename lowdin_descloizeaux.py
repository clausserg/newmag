import numpy as np


def des_cloizeaux(matrix=None, energies=None):
    ener_shifted = energies - np.mean(energies)
    # Compute overlap matrix using matrix multiplication
    S = matrix @ matrix.conj().T

    # SVD-based symmetric square root and its inverse
    U, s, Vh = np.linalg.svd(S)
    sqrt_s = np.sqrt(s)
    Ssq = U @ np.diag(sqrt_s) @ Vh
    Ssq_inv = Vh.conj().T @ np.diag(1 / sqrt_s) @ U.conj().T

    # Orthonormalize the matrix
    matrix_on = Ssq_inv.conj() @ matrix
    
    # Vectorized computation of HdC using Einstein summation
    HdC = np.einsum('k,ki,kj->ij', ener_shifted, matrix_on, matrix_on.conj())

    return HdC
