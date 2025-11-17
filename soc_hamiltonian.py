import numpy as np
from lowdin_descloizeaux import *
from helper_functions import *
from extract_caswfn_orca import *
from extract_sociwfn_orca import *
from extract_energies_orca import *
from real_to_spherical_harmonics import *
from clebsch_gordan import *
import pandas as pd


def build_heff_so(mag_center=None, output=None, level="casscf-so"):

    level = level.lower()
    if level not in {"casscf-so", "nevpt2-so"}:
        raise ValueError("level should be 'casscf-so' or 'nevpt2-so'")
    
    basis_lmls = mag_center.basis_lmls
    basis_lmlsms = mag_center.basis_lmlsms
    basis_jmj = mag_center.basis_jmj
    dim = len(basis_jmj)

    # Preload all SOC and CASSCF wavefunctions
    all_so_wfn, all_sr_wfn, energies = {}, {}, []
    if software(output) == "orca":
        all_so_wfn = soci_wfn_orca(output, level=level)  # {state_idx: {(block, sr_idx, S, Ms): coef}}
        all_sr_wfn = casscf_wfn_orca(output)  # {sr_idx: {Ml or (Ml, Ml'): coef}}
        energies = get_energies_orca(output, level=level)
    elif software(output) == "molcas":
        return "NOT implemented!"
      
    
    # Get list of (Ml, Ms) combinations in order matching basis_lmlsms
    comb_ml_ms = [(b.Ml, str(b.Ms)) for b in basis_lmlsms]
    # Initialize wavefunction matrix
    mat_wfn = np.zeros((dim, len(comb_ml_ms)), dtype=np.complex128)
    
    for state_idx, so_components in all_so_wfn.items():
        for (sr_idx, s_val, ms_val), so_coef in so_components.items():
            sr_wfn = all_sr_wfn[int(sr_idx)]

            sr_wfn = sr_wfn.copy()
            # Sign corrections (based on known ORCA quirks)
            if software(output) == "orca":
                sr_wfn[1]  = -sr_wfn.get(1, 0)   # Column 4
                sr_wfn[-1] = -sr_wfn.get(-1, 0)  # Column 2

            for ml, ml_coef in sr_wfn.items():
                row = state_idx
                col = comb_ml_ms.index((ml, str(ms_val)))
                mat_wfn[row, col] += so_coef * ml_coef

    # Build Heff in real spherical harmonics basis
    heff_rsh = des_cloizeaux(matrix=mat_wfn, energies=energies)

    # Rotate to complex spherical harmonics (CSH)
    U = harmonics_re_to_sph(basis_lmls)  # 7x7
    Uf = np.zeros_like(heff_rsh, dtype=np.complex128)
    Uf[:7, 7:] = U
    Uf[7:, :7] = U
    
    heff_csh = Uf.conj().T @ heff_rsh @ Uf

    # Rotate from |lmlsms> to |J,MJ> basis
    Ucg = cg_lml_jmj(basis_lmlsms, basis_jmj)  # dim x dim
    heff = Ucg @ heff_csh @ Ucg.conj().T

    return heff.T
