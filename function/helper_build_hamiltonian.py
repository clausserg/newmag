import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG
from collections import defaultdict

def des_cloizeaux(matrix=None, energies=None):
    '''
    Orthogonalization follwing of the 
    building Des Cloizeaux Hamiltonian
    '''
    
    ener_shifted = energies - np.mean(energies)

    # Inversion matrix S = M†M (not MM†)
    S = matrix.conj().T @ matrix

    # SVD decomposition to calculate the symmetric square root and its inverse  
    U, s, Vh = np.linalg.svd(S)
    sqrt_s = np.sqrt(s)
    Ssq = U @ np.diag(sqrt_s) @ Vh
    Ssq_inv = Vh.conj().T @ np.diag(1 / sqrt_s) @ U.conj().T

    # Orthogonalization: M_orth = M @ S^{-1/2} 
    matrix_on = matrix @ Ssq_inv

    # Construction of the effective Hamiltonian according to Des Cloizeaux:
    # HdC_ij = Σ_k (E_k * ⟨i|φ_k⟩⟨φ_k|j⟩)
    HdC = np.einsum('k,ki,kj->ij', ener_shifted, matrix_on, matrix_on.conj())

    return HdC

def harmonics_re_to_sph(basis):
    ''' 
    compute the mono-electronic matrix tranformation form RSH to CSH
    coded according to page 4 in https://doi.org/10.1016/S0166-1280(97)00185-1
    '''
    dim = len(basis)
    U = np.zeros((dim, dim), dtype=np.complex128)

    for i, m1 in enumerate(basis):
        for j, m2 in enumerate(basis):
            ml1, ml2 = m1.Ml, m2.Ml
            if abs(ml1) != abs(ml2):
                continue
            elif ml1 == 0 and ml2 == 0:
                U[i, j] = 1
            elif ml1 == ml2 and ml1 < 0:
                U[i, j] = 1j / np.sqrt(2)
            elif ml1 == ml2 and ml1 > 0:
                U[i, j] = ((-1) ** ml1) / np.sqrt(2)
            elif abs(ml1) == abs(ml2) and ml1 < 0 and ml2 > 0:
                U[i, j] = -1j * ((-1) ** ml2) / np.sqrt(2)
            elif abs(ml1) == abs(ml2) and ml1 > 0 and ml2 < 0:
                U[i, j] = 1 / np.sqrt(2)
    return U

def cg_lml_jmj(srbas, sobas):
    '''
    compute the Clebsch–Gordan matrix tranforamtion |L,ML,S,MS > to |J,MJ >
    '''
    dim = len(sobas)
    Ucg = np.zeros((dim, len(srbas)), dtype=np.complex128)

    for i, jket in enumerate(sobas):
        j, mj = float(jket.J), float(jket.Mj)
        for jdx, lket in enumerate(srbas):
            l, ml, s, ms = float(lket.L), float(lket.Ml), float(lket.S), float(lket.Ms)
            if ml + ms != mj: continue
            Ucg[i, jdx] = CG(l, ml, s, ms, j, mj).doit()

    return Ucg

def composition(eigvecs, state_index, level, mag_center, print_large):
    '''
    creation of the dict_composition[ket][%value] as a function of the basis
    '''
    vec = eigvecs[:, state_index]     
    probs = vec * vec.conj()

    composition = defaultdict(float)
    if print_large:
        if level in ['casscf-sr', 'caspt2-sr', 'nevpt2-sr']:
            kets = [f"|L={b.L}, ML={b.Ml}>" for b in mag_center.basis_lmls]
            for ket, p in zip(kets, probs):
                composition[f"{ket}"] += 100 * p                  
        elif level in ['casscf-so', 'caspt2-so', 'nevpt2-so']:
            kets = [f"|J={b.J}, MJ={b.Mj}>" for b in mag_center.basis_jmj]
            for ket, p in zip(kets, probs):
                composition[f"{ket}"] += 100 * p                  
    else:
        if level in ['casscf-sr', 'caspt2-sr', 'nevpt2-sr']:
            kets = [f"|L={b.L}, ML=±{abs(b.Ml)}>" for b in mag_center.basis_lmls]
            for ket, p in zip(kets, probs):
                composition[f"{ket}"] += 100 * p                  
        elif level in ['casscf-so', 'caspt2-so', 'nevpt2-so']:
            kets = [f"|J={b.J}, MJ=±{abs(b.Mj)}>" for b in mag_center.basis_jmj]
            for ket, p in zip(kets, probs):
                composition[f"{ket}"] += 100 * p                  
    
    composition=dict(composition)
    return composition
