import numpy as np
from class_lmlsms import *
from magnetic_moments import *
import pandas as pd

def g_tensor(mag_cntr, Hcf, doublet_index):
    """
    Compute g-tensor for a selected Kramers doublet.

    Parameters
    ----------
    Hcf : (14,14) ndarray
        Crystal-field Hamiltonian in J=5/2 ⊕ 7/2 basis.
    mux, muy, muz : (14,14) ndarrays
        Magnetic moment operators (μx, μy, μz) in the same basis.
    doublet_index : int
        Index (0-based) of the *lower* member of the nearly degenerate doublet
        in ascending eigenvalue order. The partner is index+1.

    Returns
    -------
    g_tensor : (3,3) ndarray
        Effective g-tensor in lab axes.
    g_principal : (3,) ndarray
        Principal g-values (sorted descending).
    U : (3,3) ndarray
        Rotation matrix whose columns are the principal axes in lab frame.
    """
    # Diagonalize CF Hamiltonian
    evals, evecs = np.linalg.eigh(Hcf)

    # Extract doublet states (columns of evecs)
    psi_plus = evecs[:, doublet_index]
    psi_minus = evecs[:, doublet_index + 1]

    # Build projector onto doublet subspace
    P = np.column_stack([psi_plus, psi_minus])


    mu_ops = mu_operators_lmlsms(mag_cntr)  # gives a list of [mu_x, mu_y, mu_z]
    mux, muy, muz = mu_ops[0], mu_ops[1],mu_ops[2]
    
    # Project moment operators into doublet subspace (2x2 matrices)
    mux_eff = P.conj().T @ mux @ P
    muy_eff = P.conj().T @ muy @ P
    muz_eff = P.conj().T @ muz @ P

    mu_eff = [mux_eff, muy_eff, muz_eff]

    # Pauli matrices
    sigma = [
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    ]

    # Compute g tensor
    g_tensor = np.zeros((3, 3), dtype=float)
    muB = 1.0  # we assume moment operators already contain μB factor

    for alpha in range(3):
        for beta in range(3):
            g_tensor[alpha, beta] = (2 / muB) * np.trace(mu_eff[alpha] @ sigma[beta]).real

    # Principal values from SVD of g_tensor
    U, w, Vt = np.linalg.svd(g_tensor)
    #g_principal = np.sort(s)[::-1]
    #w=np.sort(w)
    g_principal={}
    g_principal['gX']=w[0]
    g_principal['gY']=w[1]
    g_principal['gZ']=w[2]
    
    return g_principal
