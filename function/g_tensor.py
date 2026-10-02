import numpy as np
import pandas as pd
from function.class_function import *
from function.helper_build_hamiltonian import *
from function.constants import *

def mu_operators_lmlsms(mag_cntr):
    sr_basis = mag_cntr.basis_lmlsms
    so_basis = mag_cntr.basis_jmj
    U = cg_lml_jmj(sr_basis, so_basis)  # CG matrix, lmlsms -> jmj
    dim = len(sr_basis)  # dimension of the basis set
    
    mux_lmlsms = np.zeros((dim,dim), dtype=complex)
    muy_lmlsms = np.zeros((dim,dim), dtype=complex)
    muz_lmlsms = np.zeros((dim,dim), dtype=complex)
    
    # Precompute Lx, Ly and Lz for each lmlsms basis ket
    lx_dict = {idx: sr_basis[idx].op_Jx() for idx in range(dim)}
    ly_dict = {idx: sr_basis[idx].op_Jy() for idx in range(dim)}
    lz_dict = {idx: sr_basis[idx].op_Jz() for idx in range(dim)}

    # Precompute Sx, Sy and Sz for each lmlsms basis ket
    sx_dict = {idx: sr_basis[idx].op_Sx() for idx in range(dim)}
    sy_dict = {idx: sr_basis[idx].op_Sy() for idx in range(dim)}
    sz_dict = {idx: sr_basis[idx].op_Sz() for idx in range(dim)}

    # Loop over all pairs of basis elements
    for idx in range(dim):
        for jdx in range(dim):
            # Get the corresponding L and S expansion for the ket
            lx_terms, sx_terms = lx_dict[jdx], sx_dict[jdx]  # precomputed above
            ly_terms, sy_terms = ly_dict[jdx], sy_dict[jdx]  # precomputed 
            lz_terms, sz_terms = lz_dict[jdx], sz_dict[jdx]  # precomputed
            
            mux_lmlsms[idx, jdx] = (-1) * muB * (inner_prod(sr_basis[idx], lx_terms) + 2.0023 * inner_prod(sr_basis[idx], sx_terms))
            muy_lmlsms[idx, jdx] = (-1) * muB * (inner_prod(sr_basis[idx], ly_terms) + 2.0023 * inner_prod(sr_basis[idx], sy_terms))
            muz_lmlsms[idx, jdx] = (-1) * muB * (inner_prod(sr_basis[idx], [lz_terms]) + 2.0023 * inner_prod(sr_basis[idx], [sz_terms]))

    # transform into the J, Mj basis
    mux_jmj = U @ mux_lmlsms @ U.conj().T
    muy_jmj = U @ muy_lmlsms @ U.conj().T
    muz_jmj = U @ muz_lmlsms @ U.conj().T
    return [mux_jmj, muy_jmj, muz_jmj]

def inner_prod(bas_vec, myvec):
    """
    Compute the inner product between bas_vec and myvec.

    Parameters:
        bas_vec: A JMJ or LMLSMS ket (treated as 'bra').
        myvec: A list of JMJ or LMLSMS kets (treated as 'ket').

    Returns:
        result: The inner product value.
    """
    result = 0

    if isinstance(bas_vec, JMJ) and isinstance(myvec, JMJ):
        if bas_vec.Ket == myvec.Ket:
            result += bas_vec.coef.conjugate() * myvec.coef
        return result
    
    if isinstance(bas_vec, JMJ) and isinstance(myvec[0], LMLSMS):
        for idx in bas_vec.lmlsms_expansion:
            for jdx in myvec:
                if idx.Ket == jdx.Ket:
                    result += idx.coef.conjugate() * jdx.coef
        return result
                    
    if isinstance(bas_vec, LMLSMS) and isinstance(myvec[0], LMLSMS):
        for idx in myvec:
            if idx.Ket == bas_vec.Ket:
                result += bas_vec.coef.conjugate() * idx.coef
        return result
        
    if isinstance(bas_vec, JMJ) and isinstance(myvec[0], JMJ):
        for idx in myvec:
            if bas_vec.Ket == idx.Ket:
                result += bas_vec.coef.conjugate() * idx.coef
        return result

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
    P = evecs[:, [doublet_index, doublet_index+1]]

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
    g_principal['g1']=w[0]
    g_principal['g2']=w[1]
    g_principal['g3']=w[2]
    
    return g_principal


