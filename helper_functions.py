import re
import numpy as np
import pandas as pd
from state_composition import *
from class_lmlsms import *
from class_jmj import *
from dictionaries import *
from aniso_parameters import *
from tabulate import tabulate

def pretty_composition(H_cf, level):
    compo_wft={}
    Energies={}
    for index_root in range(len(H_cf)):
        compo_wft[f'Root {index_root}']={}
        Energies[index_root]=[]
        Energies[index_root],compo_wft[f'Root {index_root}']=composition(H_cf, index_root, level)
        compo_wft[f'Root {index_root}']['Energies (cm**-1)']=Energies[index_root]
    
    df=pd.DataFrame(compo_wft)
    df=df.map(lambda x: f"{x.real:6.2f}")

    return df

def pretty_matrix(matrix, basis):
    """
    Pretty-print a matrix with bras and kets as row/column labels.

    Parameters
    ----------
    matrix : np.ndarray
        Square numpy matrix.
    basis : list of objects
        Each object should have attributes .l and .ml
        (e.g. basis[i].l, basis[i].ml)
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square.")

    if matrix.shape[0] != len(basis):
        raise ValueError("Matrix size and basis size do not match.")

    kets = []
    bras = []
    # Build ket and bra labels
    if type(basis[0]) == LMLSMS:
      kets = [f"|{b.L}, {b.Ml}>" for b in basis]
      bras = [f"<{b.Ml}, {b.L}|" for b in basis]
    elif type(basis[0]) == JMJ:
      kets = [f"|{b.J}, {b.Mj}>" for b in basis]
      bras = [f"<{b.Mj}, {b.J}|" for b in basis]

    # Create DataFrame
    df = pd.DataFrame(matrix, index=bras, columns=kets)
    df = df.map(lambda x: f"{x.real:6.2f}{x.imag:+6.2f}j")

    return df


def myeig(speig):
    speig = np.asarray(speig)
    return np.sort(speig - np.min(speig))


def split_real_imag(matrix):
    """
    Split a matrix into its real and imaginary parts.

    Parameters:
    matrix (sp.Matrix): The input matrix (symbolic or numeric).

    Returns:
    tuple: Two matrices - real part and imaginary part.
    """
    matrix = sp.Matrix(matrix)
    real_part, imag_part = matrix.as_real_imag()
    return real_part, imag_part


def software(filename):
    with open(filename, mode='r') as rfile:
        content = rfile.readlines()
    for line in content:
        if "molcas".lower() in line or "molcas".upper() in line:
            return "molcas"
        elif "orca".upper() in line or "orca".lower() in line:
            return "orca"


def print_cf_params(param_sr, param_so, outpute, mag):
    # parameter groups
    b2 = [B22m, B21m, B20, B21, B22]
    b4 = [B44m, B43m, B42m, B41m, B40, B41, B42, B43, B44]
    b6 = [B66m, B65m, B64m, B63m, B62m, B61m, B60, B61, B62, B63, B64, B65, B66]

    # print parameters
    dict_paramerters={}
    dict_paramerters['Spin-Free']={}
    dict_paramerters['Spin-Orbit']={}
    for group in (b2, b4, b6):
        for par in group:
                dict_paramerters['Spin-Free'][par]=param_sr[par]
                dict_paramerters['Spin-Orbit'][par]=param_so[par]
    with open(outpute, "r", encoding="utf-8") as f:
        for line in f:       
            line = line.rstrip()
            header_line = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L =  3"
            if header_line in line:
                dict_paramerters['single_ansio_L']={}
                dict_paramerters['single_ansio_L']=aniso_params(outpute,mag)    
    df=pd.DataFrame(dict_paramerters)                     
    df=df.map(lambda x: f"{x.real:6.2f}")
    return df


def verify_cfprms_soc(mag_center=None, cfp_values=None, soc_zeta=None):
    sobas = mag_center.basis_lmlsms
    hmod = sp.Matrix.zeros(len(sobas), len(sobas))  # initialize the model SO Hamiltonian
    # let us fill the Hammiltonain
    for key, val in Bkq.items():
        hmod += matrix_cf(mag_center, sobas, key) * (cfp_values[val])
    Ucg = cg_lml_jmj(mag_center.basis_lmlsms, mag_center.basis_jmj)
    newmagH = (Ucg @ hmod @ Ucg.conj().T) + np.array(matrix_soc(mag_center.basis_jmj, soc_zeta), dtype=np.complex64)
    return np.array(newmagH, dtype=np.complex64)
