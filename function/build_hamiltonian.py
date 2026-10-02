import numpy as np
import pandas as pd
from fractions import Fraction
from scipy.sparse import csr_matrix, kron

from function.helper_functions import *
from function.helper_build_hamiltonian import *

def build_heff(sr_wfn, so_wfn, sr_energies, so_energies, mag_center, output, level, decimal, print_large):
    '''
    Compute the sr and so effectif hamiltonienne 
    
    INPUT
    sr_wfn : dict[sr_root][coeff], sr wavefunction
    so_wfn : dict[so_root][(str(sr_root), str(S), str(Ms))[coeff], so wavefunction 
    sr_energies : list of sr energies
    so_energies : list of so energies
    mag_center : class of functions possessing all the properties characteristic of the metallic center 
    output : molcas/orca outpute (useless here)
    level : level of calculation
    decimal : number of decimal printed in the NewMag out
    print_large : user option for large print newmag out

    OUTPUT
    sr_heff : sr effectif hamiltonienne 
    so_heff : so effectif hamiltonienne 
    '''

    # initialization of differents basis set
    basis_lmls = mag_center.basis_lmls
    basis_lmlsms = mag_center.basis_lmlsms
    basis_jmj = mag_center.basis_jmj
    dim = len(basis_jmj)

    #tranform dict[sr_root][coeff] to a matrix 
    ne_ml=mag_center.basis_ne_ml_uncoupled
    mat_wfn_sr = np.zeros((len(sr_wfn),len(ne_ml)), dtype=np.float64 )
    for roots in range(len(sr_wfn)):
        for base in range(len(sr_wfn[roots])):
            mat_wfn_sr[roots,base]=sr_wfn[roots][ne_ml[base]]

    # initialization of mono-electronic matrix orbital tranformation
    # form the RSH representation to the CSH representation
    U1 = harmonics_re_to_sph(MagneticCenter(mag_center.l, nel=1).basis_lmls)
    U1_sparse = csr_matrix(U1) # sparse of U1 to take less time 
    U = U1_sparse
    # application of U = U1 \bigotimes^{n-1} U1 (n is the number of elecron/hole)
    if mag_center.nel < (2*mag_center.l+1):
        for n in range(1,mag_center.nel):
            U = kron(U, U1_sparse, format="csr")
    elif mag_center.nel > (2*mag_center.l+1):
        for n in range(1,mag_center.nel-2*(mag_center.nel-int(2*mag_center.l+1))): 
            U = kron(U, U1_sparse, format="csr")
    mat_wfn_sr = mat_wfn_sr @ U # apply U on the sr_wfn

    # reading of the Clebsch–Gordan matrix tranformation to get the |L,ML> representation 
    # |l1,ml1,l2,ml2....,ln,mln > to |L,ML>
    if mag_center.l==3:    
        Ucg = np.genfromtxt(f"function/Ucg_mat/Ucg_f{mag_center.nel}", delimiter=" ")
    if mag_center.l==2:                             
        Ucg = np.genfromtxt(f"function/Ucg_mat/Ucg_d{mag_center.nel}", delimiter=" ")
    if mag_center.l==1:                             
        Ucg = np.genfromtxt(f"function/Ucg_mat/Ucg_p{mag_center.nel}", delimiter=" ")
    mat_wfn_sr = mat_wfn_sr @ Ucg #apply the Ucg on the sr_wfn

    #compute the des cloizeaux hamiltonienne for the sr wfn 
    sr_heff = des_cloizeaux(mat_wfn_sr, sr_energies)
    sr_level=level+'-sr'
    sr_level_print='SR-'+level.upper()
    # printing the result
    print(40*'-')
    print('SCALAR extraction part:')
    print(40*'-')                  
    print(f"\n{sr_level_print} numerical Hamiltonian:")   
    pretty_matrix(sr_heff, basis_lmls, mag_center, decimal)
    print(f"\n{sr_level_print} composition of states:")    
    pretty_composition(sr_energies, sr_heff, sr_level, mag_center, decimal, print_large)

    # if there are not so part -> finich 
    if so_wfn == None:
        return sr_heff, None
    else:
        # expression of so_wfn in |L,ML,S,MS> representation 
        comb_ml_ms = [(b.Ml, str(Fraction(b.Ms))) for b in basis_lmlsms]
        mat_wfn_so = np.zeros((dim, len(comb_ml_ms)), dtype=np.complex128)
        for state_idx, so_components in so_wfn.items():
            for (sr_idx, s_val, ms_val), so_coef in so_components.items():
                sr_wfn_tmp = mat_wfn_sr[int(sr_idx)]
                sr_wfn = list(sr_wfn_tmp.copy())
                for ml, ml_coef in zip(basis_lmls,sr_wfn):
                    row = state_idx
                    col = comb_ml_ms.index((ml.Ml, str(ms_val)))
                    mat_wfn_so[row, col] += so_coef * ml_coef

        #apply Clebsch–Gordan transformation 
        # |L,ML,S,MS> -> |J,MJ>
        Ucg_jmj = cg_lml_jmj(basis_lmlsms, basis_jmj) 
        mat_wfn_so = mat_wfn_so @ Ucg_jmj.T.conj() 

        #compute the des cloizeaux hamiltonienne for the so wfn 
        so_heff = des_cloizeaux(mat_wfn_so, so_energies)
        so_level=level+'-so'
        so_level_print='SO-'+level.upper()
        # printing the result
        print(40*'-')                 
        print('SOC extraction part:') 
        print(40*'-')                 
        print(f"\n{so_level_print} numerical Hamiltonian:")   
        pretty_matrix(so_heff, basis_jmj, mag_center, decimal)
        print(f"\n{so_level_print} composition of states:")    
        pretty_composition(so_energies, so_heff, so_level, mag_center, decimal, print_large)

    return sr_heff, so_heff

