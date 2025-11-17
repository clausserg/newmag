import numpy as np
from lowdin_descloizeaux import *
from helper_functions import *
from extract_caswfn_orca import *
from extract_energies_orca import *
from real_to_spherical_harmonics import *
from tabulate import tabulate
import pandas as pd


def build_heff_sr(mag_center, output, level="casscf-sr"):

    level = level.lower()
    if level not in {"casscf-sr", "nevpt2-sr"}:
        raise ValueError("level should be 'casscf-sr', 'nevpt2-sr'")

    basis = mag_center.basis_lmls
    dim = len(basis)

    # initialize a dict to hold wfn data
    all_wfns, energies = {}, []
    
    # Load all wavefunctions for the 7 states
    if software(output) == "orca":
        all_wfns = casscf_wfn_orca(output)  # returns a dict: {state_index: {Ml: coef, ...}}
        energies = get_energies_orca(output, level=level)
    elif software(output) == "molcas":
        return "NOT implemented!"

    # Initialize wavefunction matrix
    mat_wfn = np.zeros((len(all_wfns), dim), dtype=np.complex128)

    for i, wfn in all_wfns.items():
        mat_wfn[i, :] = [wfn.get(b.Ml, 0) for b in basis]

    # Apply ORCA phase corrections
    if software(output) == "orca":
        mat_wfn[:, 2] *= -1
        mat_wfn[:, 4] *= -1

    # Build Heff in real spherical harmonics (RSH) basis
    heff_rsh = des_cloizeaux(matrix=mat_wfn, energies=energies)

    # Transform to complex spherical harmonics (CSH) basis
    U = harmonics_re_to_sph(basis)

    heff_csh = U.conj().T @ heff_rsh @ U
    
    return heff_csh.T
