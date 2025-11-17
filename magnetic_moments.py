from clebsch_gordan import *
from constants import *
from inner_product import *

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

