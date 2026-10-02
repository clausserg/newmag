from function.helper_functions import *
from fractions import Fraction
from itertools import product
import re
from collections import defaultdict
import pandas as pd
import h5py
from copy import deepcopy
import numpy as np


def soci_wfn_molcas(filename, level, mag_center, root_sr_remove, user_so_root):
    '''
    reading the so wavefunction form molcas output
    
    INPUT
    filename : main molcas output 
    level : calculation level 
    mag_center : class of functions possessing all the properties characteristic of the metallic center
    root_sr_remove : sr root not selected in the model
    user_so_root : so root selected by the user

    OUTPUT
    so_wfn : dict[so_root][(str(sr_root), str(S), str(Ms))[coeff], so wavefunction 
    soc_root_remove : so root not selected in the model
    '''
    if level in ['casscf-so']:
        root_cas=number_of_CI_CAS_molcas(filename, mag_center)      
    elif level in ['caspt2-so']:
        root_cas=number_of_CI_PT2_molcas(filename, mag_center) 
    all_sr_root=list(range(0,root_cas+1))
    Ms_list = [-(mag_center.S - idx) for idx in range(int(2*mag_center.S + 1))]
    r_s_ms = list(product(all_sr_root,Ms_list))
    for i in range(len(r_s_ms)):
        r_s_ms[i]=list(r_s_ms[i])
        r_s_ms[i][0]=str(r_s_ms[i][0])
        r_s_ms[i][1]=str(r_s_ms[i][1])
        r_s_ms[i].insert(1, str(mag_center.S))
        r_s_ms[i]=tuple(r_s_ms[i])
    r_s_ms=tuple(r_s_ms)
    sr_root_select=sorted(list(set(all_sr_root)-set(root_sr_remove)))

    try:
        #soc extration
        marker_imag = "SOS_COEFFICIENTS_IMAG"
        marker_real = "SOS_COEFFICIENTS_REAL"
        h5_file=filename.split(".")
        del h5_file[-1]
        h5_file=str('.'.join(h5_file))
        with h5py.File(h5_file+'.rassi.h5', "r") as f:
                matrice_imag = f[marker_imag][()]  
                matrice_real = f[marker_real][()]  # Lecture de la dataset en numpy array
        matrice=matrice_real.T+matrice_imag.T*1j
        all_soc_wfn = defaultdict(dict)
        all_sfs=list(range(len(matrice)))
    
        for froot in range(len(matrice)):
                for tuple_root_ms,sfs in zip(r_s_ms, all_sfs):
                        all_soc_wfn[froot][tuple_root_ms] = matrice[sfs,froot]

    except FileNotFoundError:
        print("""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WARNING:
For OpenMolcas, to get the SO wavefunction, we recommend to use the the h5 file
(Except if there are two CASCI in two S-manifold in the output, BUT you MUST use MOLCAS_PRINT = 5)
Put the h5 in the same foler of your output we this type name: $Project.rassi.h5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!""")
        all_soc_wfn = defaultdict(dict)
        with open(filename, "r") as f:
            current_roots = []
            for line in f.read().splitlines():
                marker_SO = re.search(r"SFS\s+S\s+Ms\s+(.+)$",line.strip())
                if marker_SO:
                    roots = [int(x) for x in marker_SO.group(1).split()]
                    current_roots = roots
                    continue
                marker_SO = re.match(r"^(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(.+)$", line.strip())
                if marker_SO and current_roots:
                    froot = str(int(marker_SO.group(1))-1)
                    S = str(Fraction(float(marker_SO.group(2))).limit_denominator())
                    Ms = str(Fraction(float(marker_SO.group(3))).limit_denominator())
                    coeffs = marker_SO.group(4).split("(")
                    for root, coeff in zip(roots, coeffs[1:]):
                        real, imag = coeff.strip(") ").split(",")
                        cplx = np.complex128(np.float64(real), np.float64(imag))
                        all_soc_wfn[int(root)-1][(froot, S, Ms)]= cplx

        Ms_list = [-(mag_center.S - idx) for idx in range(int(2*mag_center.S + 1))]
        for soc_root in list(range(0,  len(all_soc_wfn), 1)):           
            for sr_root in all_sr_root:      
                for Ms in Ms_list:                                 
                    if not (str(sr_root), str(mag_center.S), str(Ms)) in all_soc_wfn[soc_root]:
                        all_soc_wfn[soc_root][(str(sr_root), str(mag_center.S),str(Ms))]=np.complex128(0)
                
    all_soc_wfn=dict(all_soc_wfn)
    if len(all_soc_wfn)==0:
        return None, [] 

    select_soc_wfn=deepcopy(all_soc_wfn)
    if user_so_root == None:
        soc_root_selected=select_soc_root(all_soc_wfn, len(all_soc_wfn), sr_root_select, mag_center,printing=True)
    else:
        print('Spin-orbit roots selected by user input')                       
        wfn_print={}
        wfn_print['User roots']={}
        for root in range(len(all_soc_wfn)):
            if root in user_so_root:
                wfn_print['User roots'][f'Root {root}']='---> Selected'
            else:
                wfn_print['User roots'][f'Root {root}']=' '
        df=pd.DataFrame(wfn_print)
        print(df.to_string(),end='\n\n')
        soc_root_selected = user_so_root

    soc_root_remove=list(set(list(range(len(all_soc_wfn))))-set(soc_root_selected))
    
    for froot in range(len(all_soc_wfn)): 
        for del_sr_root in range(len(root_sr_remove)):
            for Ms in Ms_list:
                del select_soc_wfn[froot][(str(root_sr_remove[del_sr_root]), str(mag_center.S),str(Ms))]
    
    for del_soc_root in range(len(soc_root_remove)):
        del select_soc_wfn[soc_root_remove[del_soc_root]]
    
    for key_soc, new_key in zip(soc_root_selected,range(len(mag_center.basis_jmj))):
        select_soc_wfn[new_key]=select_soc_wfn.pop(key_soc)
    
    for froot in range(len(mag_center.basis_jmj)) :
        for key_sr, new_key in zip(sr_root_select,range(root_cas)):
            for Ms in Ms_list:
                select_soc_wfn[froot][(str(new_key),str(mag_center.S),str(Ms))]=select_soc_wfn[froot].pop(((str(key_sr), str(mag_center.S),str(Ms))))
    
    so_wfn={}
    for froot in range(len(mag_center.basis_jmj)) :
        so_wfn[froot]={}
    
    for Ms in Ms_list:
        for froot in range(len(mag_center.basis_jmj)) :
            for root in range(len(sr_root_select)):
                so_wfn[froot][(str(root), str(mag_center.S),str(Ms))]=select_soc_wfn[froot][(str(root),str(mag_center.S),str(Ms))]
    
    return so_wfn, soc_root_remove

