'''
File with different util python function 
-> specific function for orca/molcas outpute
-> for printing newmag out
-> SOC analysis
-> general ERROR orca/molcas 
'''
import re
import numpy as np
import math
import pandas as pd
from fractions import Fraction
from tabulate import tabulate
from itertools import permutations
from itertools import combinations

from function.helper_build_hamiltonian import *
from function.class_function import *
#from function.constants import *
from function.crystal_field import *
from function.g_tensor import inner_prod

################################
# helper function for software #
################################

def software(filename):
    # determine the software output -> molcas or orca
    with open(filename, mode='r') as rfile:
        content = rfile.readlines()
    for line in content:
        if "molcas".lower() in line or "molcas".upper() in line:
            return "molcas"
        elif "orca".upper() in line or "orca".lower() in line:
            return "orca"

def generate_patern_particles(length_shell, number_particles):
   '''
   generation of the CSF partern for the metal

    INPUTE
    length_shell : lenght of the metal shell (f=7, d=5)
    number_particles : number of electron/hole for the metal center

    OUTPUT
    results : list of the formating CSF for the metal center 
   '''
   if number_particles < length_shell:
        results = []
        for position in combinations(range(length_shell), number_particles):
            sub = ['0'] * length_shell
            for p in position:
                sub[p] = 'u'
            results.append(''.join(sub))
        return results
 
   if number_particles > length_shell:
        results = []
        for position in combinations(range(length_shell), number_particles-length_shell):
            sub = ['u'] * length_shell
            for p in position:
                sub[p] = '2'
            results.append(''.join(sub))
        return results

def pattern_caswfn(cas_mo, mag, software):
    '''
    generation of the general CSF partern 

    INPUT
    cas_mo : formating of the CSF
    mag : class of functions possessing all the properties characteristic of the metallic center 
    software : type of software (orca or molcas)

    OUTPUT
    pattern_wfn : string with the general CSF partern
    remove_conf: list of CSFs not on the metal center
    '''
    pattern_wfn=[]
    if software == "molcas":
        pattern_wfn=[r"\s+\d+\s+"]
    pattern_center = generate_patern_particles(2*mag.l+1, mag.nel)
    list_type_cas_mo=[] 
    center_mo=0
    for mo in cas_mo:
        nb_type_mo = list(mo)
        type_mo = nb_type_mo[-1]
        del nb_type_mo[-1]
        tmp = ''.join(nb_type_mo)
        nb_mo=int(tmp)
        if  type_mo=='d': # CSF with always an occupency of 2
            for _ in range(nb_mo):
                pattern_wfn.append(r'(2)')
                list_type_cas_mo.append('d')
        elif type_mo=='v': # CSF with always an occupency of 0  
            for _ in range(nb_mo):
                pattern_wfn.append(r'(0)')
                list_type_cas_mo.append('v')
        elif type_mo=='o': # CSF with a parical occupency but not is the interst metal/shell
            pattern_wfn.append(f'[ud02]{{{nb_mo}}}')
            for _ in range(nb_mo):
                list_type_cas_mo.append('o')
        elif type_mo=='c': # CSF of interest
            pattern_wfn.append(r'(')
            for pat in pattern_center:
                p = list(pat)
                for i in range(center_mo, nb_mo+center_mo):
                    pattern_wfn.append(p[i])
                pattern_wfn.append(r'|')
            center_mo+=nb_mo
            del pattern_wfn[-1]
            pattern_wfn.append(r')')
            for _ in range(nb_mo):
                list_type_cas_mo.append('c')
    pattern_wfn=''.join(pattern_wfn)

    remove_conf=[]
    for i, n in enumerate(list_type_cas_mo):
        if n == 'd' or n == 'v' or n == 'o':
            remove_conf.append(i)
    return pattern_wfn, remove_conf

def caspt2_molcas_present(filename):
    '''
    check if caspt2 calculation in the molcas out
    '''
    with open(filename, mode='r') as rfile:
        content = rfile.readlines() 
        for line in content:            
            if "&caspt2".lower() in line or "&caspt2".upper() in line:
                return True
        return False
def nevpt2_orca_present(filename):
    '''
    check if nevpt2 calculation in the orca out
    '''
    with open(filename, mode='r') as rfile:
        content = rfile.readlines() 
        for line in content:            
            if "< NEVPT2  >" in line:
                return True
                break
            elif "TOTAL RUN TIME:"  in line:
                return False

def type_caspt2_molcas(filename):
    '''
    check the type of caspt2 calculation in the molcas out
    ss-caspt2 or ms-caspt2
    '''
    with open(filename, mode='r') as rfile:
        content = rfile.readlines() 
        for line in content:            
            if "MS-CASPT2".lower() in line or "MS-CASPT2".upper() in line:
                return 'ms-caspt2-sr'
                break
            elif "SS-CASPT2".lower() in line or "SS-CASPT2".upper() in line:
                return 'ss-caspt2-sr'
                break

def ml_order_molcas(filename, mag_center):
    '''
    reading the ml order for molcas out 
    the reading is done after the CASCI block
    '''
    ml_order = []                
    marker_inactive_OM = False   
    marker_active_OM = False     
    marker_CASCI = False  
    marker_print_OM = False       
    
    with open(filename, "r") as f:
        for line in f:
            if not marker_inactive_OM:
                if "Inactive orbitals" in line:
                    nomber_inactive_OM = int(line.split()[-1])
            if not marker_active_OM:
                if 'Active orbitals' in line:
                     nomber_active_OM = int(line.split()[-1])
            if not marker_CASCI:
                if "CASCI only, no orbital optimization will be done" in line:
                    marker_CASCI = True
                continue
            if marker_CASCI and not marker_print_OM:
                if re.search(r"^\s+\d+\s+0.0000\s+0.0000\s+", line):
                    marker_print_OM = True
                continue
            if re.search(r"^\s+"+str(nomber_inactive_OM+nomber_active_OM+1)+r"\s+0.0000\s+0.0000\s+", line):break 

            if mag_center.l==3 :
                if re.search(r"\b[A-Z]+\d*\b", line) and re.search(r"(4f|5f)", line):
                    pattern = re.compile(r'\df(\d)([ +-])')
                    match = pattern.search(line)
                    if match:
                        number = int(match.group(1)) 
                        sign = match.group(2)
                        ml = number if sign == "+" else -number
                        ml_order.append(ml)
                        if len(ml_order)==(2*mag_center.l+1):break

            if mag_center.l==2 :
                if re.search(r"\b[A-Z]+\d*\b", line) and re.search(r"(3d|4d|5d|6d)", line):
                    pattern = re.compile(r'\dd(\d)([ +-])')
                    match = pattern.search(line)
                    if match:
                        number = int(match.group(1))  
                        sign = match.group(2)        
                        ml = number if sign == "+" else -number
                        ml_order.append(ml)
                        if len(ml_order)==(2*mag_center.l+1):break

            if mag_center.l==1 :
                print('You want to study pn configuration')
                print('You need to use the --ml_order option in the parser input, specifying the correct order of the ml numbers')
                print('pz -> ml = 0')
                print('px -> ml = 1')
                print('py -> ml =-1')
                exit()
 

        if len(set(ml_order)) != 2*mag_center.l+1:
            print('ERROR: Failed to read ml order... ')
            print('You must use the --ml_order option in the parser input, specifying the correct order of the ml numbers')
            general_error_inp_molcas(mag_center)
    return ml_order

def select_sr_root_by_OA(weight_wfn, mag_center, printing=False):
    '''
    Selection of the sr root for orca/molcas out
    the selection is done using the sum of weight of each root in the CSF of interest
    if the sum of weight is upper of a criteria the root is selection as a interest root to build the Heff

    INPUTE
    weight_wfn : dict[sr_root][sum weight on interest CSF]
    mag_center : class of functions possessing all the properties characteristic of the metallic center 
    printing : printing the result

    OUTPUT
    root_select : list of selected root
    '''
    root_weight=[]
    for i in range(len(weight_wfn)):      
        root_weight.append(weight_wfn[i])
    root_select=[]
    print(f'Root selected if the total weight is upper to {np.round(100/(2*mag_center.l+1),2)}%')
    for i in range(len(weight_wfn)):
        if 1/(2*mag_center.l+1)<weight_wfn[i]:
            root_select.append(i)
    
    weight_wfn_print={}
    for root in range(len(weight_wfn)):
        if root in root_select:
            weight_wfn_print[f'Root {root}']=str(np.round(weight_wfn[root]*100,2))+' % ---> Selected'
        else:
            weight_wfn_print[f'Root {root}']=str(np.round(weight_wfn[root]*100,2))+' %'
    if printing:                                                         
        normalisation={}
        if 2*mag_center.l+1==7:
            normalisation['Scalar level: Total weight on f OA']=weight_wfn_print
        elif 2*mag_center.l+1==5:
            normalisation['Scalar level: Total weight on d OA']=weight_wfn_print
        elif 2*mag_center.l+1==3:
            normalisation['Scalar level: Total weight on p OA']=weight_wfn_print
 
        df=pd.DataFrame(normalisation)
        print(df,end='\n\n')
    return sorted(list(root_select))


def number_of_CI_CAS_molcas(filename, mag_center):
    '''
    return the number of root is the cassscf calculation for molcas out
    '''
    marker_CSFs = False     
    root_cas=-1
    with open(filename, "r") as f:                                             
        for line in f:                                                         
            if not marker_CSFs:                                                
                if 'CASCI only, no orbital optimization will be done.' in line:
                    marker_CSFs = True
            if marker_CSFs:
                if '&RASSI' in line:
                    break
                elif 'printout of CI-coefficients larger than' in line:
                    root_cas += 1# int(line.split()[-1])  
            if "Final results" in line:
                marker_CSFs = False
            
    if root_cas== -1 :
        print('ERROR: Your output is incorrect (We do not know why, but it is incorrect)')
        general_error_inp_molcas(mag_center)    
    return root_cas

def number_of_CI_CAS_orca(filename,mag_center):
    '''
    return the number of root is the cassscf calculation for orca out
    '''
    marker_CSFs = False     
    root_cas = -1
    with open(filename, "r") as f:       
        for line in f:                   
            if not marker_CSFs:         
                if "Spin-Determinant CI Printing" in line:
                    marker_CSFs = True
            if marker_CSFs:
                if 'DENSITY MATRIX' in line or 'SA-CASSCF TRANSITION ENERGIES' in line:
                    break
                elif 'ROOT' and ":  E=" in line:
                    root_cas += 1#int(line.split()[1].split(':')[0])
            if "CAS-SCF STATES FOR BLOCK" in line:
                marker_CSFs=False
    if root_cas == -1:
       general_error_inp_orca(mag_center) 
    return root_cas

def number_of_CI_PT2_molcas(filename, mag_center):
    ''' 
    return the number of root is the caspt2 calculation for molcas out 
    '''
    root_cas = -1
    with open(filename, "r") as f:       
        for line in f:                   
            if "Number of CI roots used" in line:
                root_cas += int(line.split()[-1])
    if root_cas == -1:
        print('ERROR: Your output is incorrect (We do not know why, but it is incorrect)')
        general_error_inp_molcas(mag_center)    

    return root_cas

def select_soc_root(all_so_wfn, number_soc_root, sr_root_select,mag_center, printing=False):
    '''
    Selection of the so root for orca/molcas out
    the selection is done using the projection of the coupling on the sr root 
    P = sum_k <psi_so(sr) | psi_so(sr) >
     
    INPUTE
    all_so_wfn : dict[so_root][(str(sr_root), str(S), str(Ms))[coeff], so wavefunction 
    number_soc_root : number of so root
    sr_root_select :  list of selected sr root
    mag_center : class of functions possessing all the properties characteristic of the metallic center 
    printing : printing the result
     
    OUTPUT
    soc_root_selected : list of selected so root

    '''
    projection={}
    Ms_list = [-(mag_center.S - idx) for idx in range(int(2*mag_center.S + 1))]
    for soc_root in list(range(0,  number_soc_root, 1)):
            projection[soc_root]=0
            for sr_root in sr_root_select:
                for Ms in Ms_list:
                    projection[soc_root] += (all_so_wfn[soc_root][(str(sr_root), str(mag_center.S),str(Ms))] 
                             * all_so_wfn[soc_root][(str(sr_root),  str(mag_center.S),str(Ms))].conj())
    norm=[]
    for i in range(len(projection)):
            norm.append(projection[i].real)
    soc_root_selected=[]
    for j in range(len(norm)):
            for i in range(len(projection)):
                    if len(soc_root_selected)==len(mag_center.basis_jmj):
                            break
                    elif projection[i]==max(norm):
                            soc_root_selected.append(i)
                            norm.remove(max(norm)) 
    weight_wfn_print={}    
    for root in range(len(projection)):
        if root in soc_root_selected:
            weight_wfn_print[f'Root {root}']=str(np.round(projection[root].real*100,2))+' % ---> Selected'
        else:              
            weight_wfn_print[f'Root {root}']=str(np.round(projection[root].real*100,2))+' %'
    if printing: 
        normalisation={}
        normalisation['Spin-Orbit level: coupling between scalar roots']=weight_wfn_print
        df=pd.DataFrame(normalisation)
        print(df.to_string(),end='\n\n')
    return sorted(list(set(soc_root_selected)))


def ailft_in_orca(orca_output):
    '''
    check if aiLFT block is present in the orca out
    '''
    aiLFT_marker_casscf=False
    aiLFT_marker_nevpt2=False
    ailft_level=[]
    with open(orca_output, 'r') as f:
        for line in f:
            if not aiLFT_marker_casscf:
                if 'AILFT MATRIX ELEMENTS (CASSCF)' in line:
                    aiLFT_marker_casscf=True
                    ailft_level.append('casscf')
                continue
            if not aiLFT_marker_nevpt2:
                if 'AILFT MATRIX ELEMENTS (NEVPT2)' in line:
                    aiLFT_marker_nevpt2=True
                    ailft_level.append('nevpt2')
                continue
 
        if aiLFT_marker_casscf or aiLFT_marker_nevpt2:
            return True, ailft_level
        else:
            return False, ''
                


###################
# PRINTING
################### 

def pretty_composition(list_energie, H_cf, level, mag_center, decimal, print_large):
    '''
    Pretty print wavefunction composition of specific Heff 

    INPUTE
    list_energie : list of ab-initio energies for the specific calculation level
    H_cf :  Heff at specific calculation level
    level : calculation level 
    mag_center : class of functions possessing all the properties characteristic of the metallic center 
    decimal :  number of decimal printed in the NewMag out
    print_large : user option for large print newmag out
    '''
    compo_wft={}
    eigvals, eigvecs = np.linalg.eigh(H_cf)
    Energie_model=[]

    for index_root in range(len(H_cf)):
        compo_wft[f'Root {index_root}']={}
        compo_wft[f'Root {index_root}']=composition(eigvecs, index_root, level, mag_center, print_large)
    for index_root in range(len(eigvals)):
        compo_wft[f'Root {index_root}']['Energies (cm**-1)'] = list_energie[index_root]

    df=pd.DataFrame(compo_wft)
    df=df.map(lambda x: f"{x.real:6.{decimal}f}")
    if level.split('-')[-1] == 'sr':
        print(df.to_string(), end='\n\n')

    elif level.split('-')[-1] == 'so':
        basis=mag_center.basis_jmj
        J_order = sorted(list(set([float(b.J) for b in basis])))
        if  mag_center.nel > (2*mag_center.l+1):
            J_order = list(reversed(J_order))

        for J in J_order:
            J_values = [Fraction(j) for j in J_order]
            J_target = Fraction(J)

            block_sizes = [int(2*j + 1) for j in J_values]
            indices = np.cumsum([0] + block_sizes)
            i = J_values.index(J_target)
            block = df.iloc[:,indices[i]:indices[i+1]]
            print(30*'-')
            print(f'Composition of roots mainly on |J={J_target}>', end='\n\n')
            print(block.to_string(), end='\n\n')

    return

def pretty_matrix(matrix, basis, mag_center, decimal):
    """
    Pretty print of Heff in a specific basis |L,ML> or |J,MJ>

    INPUTE
    matrix :  Heff at specific calculation level
    basis : specific basis define with mag_center 
    mag_center : class of functions possessing all the properties characteristic of the metallic center 
    decimal :  number of decimal printed in the NewMag out

    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square.")

    if matrix.shape[0] != len(basis):
        raise ValueError("Matrix size and basis size do not match.")

    kets = []
    bras = []
    nb_J=[]
    # Build ket and bra labels
    if type(basis[0]) == LMLSMS:
        kets = [f"|{b.L}, {b.Ml}>" for b in basis]
        bras = [f"<{b.L}, {b.Ml}|" for b in basis]
    elif type(basis[0]) == JMJ:
        kets = [f"|{b.J}, {b.Mj}>" for b in basis]
        bras = [f"<{b.J}, {b.Mj}|" for b in basis]
        J_order = sorted(list(set([float(b.J) for b in basis])))
        if  mag_center.nel > (2*mag_center.l+1):
            J_order = list(reversed(J_order))

    # Create DataFrame
    df = pd.DataFrame(matrix, index=bras, columns=kets)
    df = df.map(lambda x: f"{x.real:6.{decimal}f}{x.imag:+6.{decimal}f}j")

    if type(basis[0]) == LMLSMS:
        print(df.to_string(), end='\n\n')

    elif type(basis[0]) == JMJ:   
        print('The matrix block composition in <J_1 | J_2> notation:') 
        block_matrix=[]
        for J_i in J_order:
            row = []
            for J_j in J_order:
                row.append(f"<{Fraction(J_j)} | {Fraction(J_i)}>")
            block_matrix.append(row)
        print(tabulate(block_matrix, tablefmt='fancy_grid'), end='\n\n')

        for J1 in J_order:
            for J2 in J_order:
                J_values = [Fraction(j) for j in J_order]
                J1 = Fraction(J1)
                J2 = Fraction(J2)
                block_sizes = [int(2*j + 1) for j in J_values]

                # Compute cumulative indices
                indices = np.cumsum([0] + block_sizes)
            
                i = J_values.index(J2)
                j = J_values.index(J1)
            
                block = df.iloc[
                    indices[i]:indices[i+1],
                    indices[j]:indices[j+1]
                ]
                print(30*'-')
                print(f'BLOC: <J={J1}, M_J | J={J2}, M_J>', end='\n\n')
                print(block.to_string(),end='\n\n')
    return

#def myeig(speig):
#    speig = np.asarray(speig)
#    return np.sort(speig - np.min(speig))
#
#def split_real_imag(matrix):
#    """
#    Split a matrix into its real and imaginary parts.
#
#    Parameters:
#    matrix (sp.Matrix): The input matrix (symbolic or numeric).
#
#    Returns:
#    tuple: Two matrices - real part and imaginary part.
#    """
#    matrix = sp.Matrix(matrix)
#    real_part, imag_part = matrix.as_real_imag()
#    return real_part, imag_part

##############################
# SOC analyses
##############################

def soc_analysies(Heff, cfp, mag, decimal, soc_xyz, print_large):
    '''
    main function of the soc analyse 
    1- compute the lamda/zeta constant using the spherical approx. 
    2- reconstruction of the Hamiltonian of crystal field and spin-orbit
    3- evaluation of the model Hamiltonian 
    
    INPUT
    Heff : dict[level][sr or so][array Hamiltonian], the all matrix Hamiltonian compute with newmag
    cfp : dict[level][sr or so][Bkq][value], the all CFPs compute with newmag
    mag : class of functions possessing all the properties characteristic of the metallic center
    decimal : number of decimal printed in the NewMag out
    print_large : user option for large print newmag out 
    '''
    b2 = [B22m, B21m, B20, B21, B22]
    b4 = [B44m, B43m, B42m, B41m, B40, B41, B42, B43, B44]
    b6 = [B66m, B65m, B64m, B63m, B62m, B61m, B60, B61, B62, B63, B64, B65, B66]
    b8 = [B88m, B87m, B86m, B85m, B84m, B83m, B82m, B81m, B80, B81, B82, B83, B84, B85, B86, B87, B88]
    b10 = [B1010m, B109m, B108m, B107m, B106m, B105m, B104m, B103m, B102m, B101m, B100, B101, B102, B103, B104, B105, B106, B107, B108, B109, B1010]
    b12 = [B1212m, B1211m, B1210m, B129m, B128m, B127m, B126m, B125m, B124m, B123m, B122m, B121m, B120, B121, B122, B123, B124, B125, B126, B127, B128, B129, B1210, B1211, B1212]
    Bkq = {         
    (2, 0): B20, (2, 1): B21, (2, 2): B22, (2,-1): B21m, (2,-2): B22m,
                
    (4, 0): B40, (4, 1): B41, (4, 2): B42, (4, 3): B43, (4, 4): B44, 
    (4,-1): B41m, (4,-2): B42m, (4,-3): B43m, (4,-4): B44m,
                
    (6, 0): B60, (6, 1): B61, (6, 2): B62, (6, 3): B63, (6, 4): B64, (6, 5): B65, (6, 6): B66, 
    (6,-1): B61m, (6,-2): B62m, (6,-3): B63m, (6,-4): B64m, (6,-5): B65m, (6,-6): B66m,
                
    (8, 0): B80, 
    (8,  1): B81, (8,  2): B82, (8,  3): B83, (8,  4): B84, (8,  5): B85, (8,  6): B86,  (8, 7): B87,  (8, 8): B88,
    (8, -1): B81m, (8, -2): B82m, (8, -3): B83m, (8, -4): B84m, (8, -5): B85m, (8, -6): B86m, (8,-7): B87m, (8,-8): B88m,
                
    (10, 0): B100,
    (10,  1): B101, (10,  2): B102, (10,  3): B103, (10,  4): B104, (10,  5): B105, (10,  6): B106,  (10, 7): B107,  (10, 8): B108, (10, 9): B109,(10, 10): B1010,
    (10,  -1): B101m, (10,  -2): B102m, (10,  -3): B103m, (10,  -4): B104m, (10,  -5): B105m, (10,  -6): B106m,  (10, -7): B107m,  (10, -8): B108m, (10, -9): B109m,(10, -10): B1010m,
                
    (12, 0): B120,
    (12,  1): B121, (12,  2): B122, (12,  3): B123, (12,  4): B124, (12,  5): B125, (12,  6): B126,  (12, 7): B127,  (12, 8): B128, (12, 9): B129,(12, 10): B1210,(12, 11): B1211, (12, 12): B1212,
    (12,  -1): B121m, (12,  -2): B122m, (12,  -3): B123m, (12,  -4): B124m, (12,  -5): B125m, (12,  -6): B126m,  (12, -7): B127m,  (12, -8): B128m, (12, -9): B129m,(12, -10): B1210m,(12, -11): B1211m,  
    (12, -12): B1212m,
                
    }               

    for level in Heff.keys():
        if level.split('-')[-1] == 'sr' and Heff[level] is not None:
            sr_level_print='SR-'+str(level.split("-")[0]).upper()
            print('------------------------------------------')
            print(f'Reconstruction at the {sr_level_print} level:')
            print('------------------------------------------')
            print("(Please be patient, this could take a little time...)",end='\n\n')
            H_ab = Heff[level]
            E_ab, wf_ab = np.linalg.eigh(H_ab)
            E_ab = sorted(E_ab - min(E_ab))
            dico_Energie={}
            dico_Energie[f"Ab initio"] = {}                         
            for i in range(len(E_ab)):
                dico_Energie[f"Ab initio"][f'root {i}'] = np.round(float(E_ab[i]),decimal)
            dico_Energie[f"Ab initio"]['------'] = '------' 
            dico_Energie[f"Ab initio"]['RMSD'] = '/'
            dico_Energie[f"Ab initio"]['MAE'] = '/'
            dico_Energie[f"Ab initio"]['MAEER'] = '/'
            H_cf = np.zeros((len(mag.basis_lmls),len(mag.basis_lmls)), dtype=np.complex128)
            O_k='O'
            for group in (b2, b4, b6, b8, b10, b12):
                keys = [key for key, val in Bkq.items() if val in group]
                O_k += str(keys[0][0])
                H_cf += cf_hamiltonian(mag_center=mag, bkq_values=cfp[level], bkq_group=keys, ls_constants=False)
                H_mod = H_cf
                E_mod, wf_mod = np.linalg.eigh(H_mod)                   
                E_mod = sorted(E_mod - min(E_mod))                            
                if print_large:
                    print(f'======Construction of H_CF({O_k})======')
                    pretty_matrix(H_mod, mag.basis_lmls, mag, decimal)
                    print('\n') 
                    print(f"======Composition of H_CF({O_k})======")
                    pretty_composition(E_mod, H_mod, level, mag, decimal, print_large) 
                dico_Energie[f"{O_k}"] = {}
                for i in range(len(E_mod)):
                    dico_Energie[f"{O_k}"][f'root {i}'] = np.round(float(E_mod[i]),decimal)
                dico_Energie[f"{O_k}"]['------'] = '------' 
                dE=[]
                dEabs=[]
                dEsr=[]
                for E_ref, E_Ok in zip(E_ab, E_mod):
                    dE.append(E_ref-E_Ok)
                    dEabs.append(np.abs(E_ref-E_Ok))
                    dEsr.append((E_ref-E_Ok)**2)
                Guihery_Nathalie_Metric= np.round(np.abs(sum(dE))/(len(H_mod)*(max(E_ab)))*100, decimal)
                RMSD_E = np.round(np.sqrt(sum(dEsr)/len(H_mod)), decimal)
                MAE_E  = np.round(sum(dEabs)/len(H_mod), decimal)
                dico_Energie[f"{O_k}"]['RMSD'] = RMSD_E
                dico_Energie[f"{O_k}"]['MAE'] = MAE_E
                dico_Energie[f"{O_k}"]['MAEER'] = Guihery_Nathalie_Metric 
                O_k += ' + O'
            if print_large:
                print(f'=================Rest: H_ab - H_mod=================')
                print(f'=====H_mod = H_CF(O2 + O4 + O6 + O8 + O10 + O12)=====')
                H_rest = H_ab - H_mod
                pretty_matrix(H_rest, mag.basis_lmls, mag, decimal)
                print('\n') 
                 
            df=pd.DataFrame(dico_Energie)
            print('\n------Reconstruction of the Energy Spectrum------',end='\n\n')
            print(tabulate(df, headers='keys',tablefmt='pipe'),end='\n\n')
            print('RMSD: Root Mean Square Deviation on energies in cm**-1' )
            print('MAE: Mean Absolute Error on energies in cm**-1')
            print('MAEER: Mean Absolute Error on Energies Relative to the ab initio spectral width in %' ,end='\n\n')

            
        if level.split('-')[-1] == 'so' and Heff[level] is not None:
            so_level_print='SO-'+str(level.split("-")[0]).upper()
            print('------------------------------------------')
            print(f'Reconstruction at the {so_level_print} level:')
            print('With H_abinito = H_SOC(Lambda) + H_CF(B_k^q * O_k^q)')
            print('------------------------------------------')
            H_ab = Heff[level]
            E_ab, wf_ab = np.linalg.eigh(H_ab)
            E_ab = sorted(E_ab - min(E_ab))
            dico_Energie={}
            dico_Energie[f"Ab initio"] = {}                         
            for i in range(len(E_ab)):
                dico_Energie[f"Ab initio"][f'root {i}'] = np.round(float(E_ab[i]),decimal)
            dico_Energie[f"Ab initio"]['------'] = '------' 
            dico_Energie[f"Ab initio"]['RMSD'] = '/'
            dico_Energie[f"Ab initio"]['MAE'] = '/'
            dico_Energie[f"Ab initio"]['MAEER'] = '/'
            if soc_xyz==False:
                Zeta, Lambda = extract_soc_iso(mag_center=mag, heff=H_ab)
                H_so = matrix_ls_jmj_iso(mag) * Lambda
                print('LAMBDA = '+str(np.round(Lambda,decimal))+' cm**-1',end='\n')
                print('ZETA = '+str(np.round(Zeta,decimal))+' cm**-1',end='\n\n')
                print("(Please be patient, this could take a little time...)")
            if soc_xyz==True:
                Zeta_xyz, Lambda_xyz = extract_soc_aniso(mag_center=mag, heff=H_ab)
                H_so_x = matrix_ls_jmj_aniso(mag,'x') * Lambda_xyz[0]
                H_so_y = matrix_ls_jmj_aniso(mag,'y') * Lambda_xyz[1]
                H_so_z = matrix_ls_jmj_aniso(mag,'z') * Lambda_xyz[2]
                H_so = H_so_x + H_so_y + H_so_z
                print('LAMBDA_x = '+str(np.round(Lambda_xyz[0],decimal))+' cm**-1',end='\n')
                print('LAMBDA_y = '+str(np.round(Lambda_xyz[1],decimal))+' cm**-1',end='\n')
                print('LAMBDA_z = '+str(np.round(Lambda_xyz[2],decimal))+' cm**-1',end='\n')
                print('LAMBDA_xyz = '+str(np.round((Lambda_xyz[0] + Lambda_xyz[1] + Lambda_xyz[2])/3,decimal))+' cm**-1',end='\n\n')

                print('ZETA_x = '+str(np.round(Zeta_xyz[0],decimal))+' cm**-1',end='\n')
                print('ZETA_y = '+str(np.round(Zeta_xyz[1],decimal))+' cm**-1',end='\n')
                print('ZETA_z = '+str(np.round(Zeta_xyz[2],decimal))+' cm**-1',end='\n')
                print('LAMBDA_xyz = '+str(np.round((Zeta_xyz[0] + Zeta_xyz[1] + Zeta_xyz[2])/3,decimal))+' cm**-1',end='\n\n')
                print("(Please be patient, this could take a little time...)")
                
            if print_large:
                print(f'======Construction of H_SO======')
                pretty_matrix(H_so, mag.basis_jmj, mag, decimal)
                print('\n') 
                E_so, wf_so = np.linalg.eigh(H_so)
                E_so = sorted(E_so - min(E_so))
                print(f"======Composition of H_SO======")
                pretty_composition(E_so, H_so, level, mag, decimal, print_large)

            E_mod, wf_mod = np.linalg.eigh(H_so)
            E_mod = sorted(E_mod - min(E_mod))
            dico_Energie["H_SO"] = {}
            for i in range(len(E_mod)):
                dico_Energie["H_SO"][f'root {i}'] = np.round(float(E_mod[i]),decimal)
            dico_Energie["H_SO"]['------'] = '------' 
            dE=[]
            dEabs=[]
            dEsr=[]
            for E_ref, E_Ok in zip(E_ab, E_mod):
                dE.append(E_ref-E_Ok)
                dEabs.append(np.abs(E_ref-E_Ok))
                dEsr.append((E_ref-E_Ok)**2)
            Guihery_Nathalie_Metric= np.round(np.abs(sum(dE))/(len(dE)*(max(E_ab)))*100, decimal)
            RMSD_E = np.round(np.sqrt(sum(dEsr)/len(dEsr)), decimal)
            MAE_E  = np.round(sum(dEabs)/len(dEabs), decimal)
            dico_Energie["H_SO"]['RMSD'] = RMSD_E
            dico_Energie["H_SO"]['MAE'] = MAE_E                                                    
            dico_Energie["H_SO"]['MAEER'] = Guihery_Nathalie_Metric 

            
            H_cf = np.zeros((len(mag.basis_jmj),len(mag.basis_jmj)), dtype=np.complex128)
            O_k='O'
            for group in (b2, b4, b6, b8, b10, b12):
                keys = [key for key, val in Bkq.items() if val in group]
                O_k += str(keys[0][0])
                H_cf += cf_hamiltonian(mag_center=mag, bkq_values=cfp[level], bkq_group=keys, ls_constants=True)
                H_mod = H_cf + H_so
                E_mod, wf_mod = np.linalg.eigh(H_mod)                   
                E_mod = sorted(E_mod - min(E_mod))
                if print_large:
                    print(f'======Construction of H_SOC + H_CF({O_k})======')
                    pretty_matrix(H_mod, mag.basis_jmj, mag, decimal)
                    print('\n') 
                    print(f"======Composition of H_SOC + H_CF({O_k})======")
                    pretty_composition(E_mod, H_mod, level, mag, decimal, print_large) 
                dico_Energie[f"H_SO + {O_k}"] = {}
                for i in range(len(E_mod)):
                    dico_Energie[f"H_SO + {O_k}"][f'root {i}'] = np.round(float(E_mod[i]),decimal)
                dico_Energie[f"H_SO + {O_k}"]['------'] = '------' 
                dE=[]
                dEabs=[]
                dEsr=[]
                for E_ref, E_Ok in zip(E_ab, E_mod):
                    dE.append(E_ref-E_Ok)
                    dEabs.append(np.abs(E_ref-E_Ok))
                    dEsr.append((E_ref-E_Ok)**2)
                Guihery_Nathalie_Metric= np.round(np.abs(sum(dE))/(len(dE)*(max(E_ab)))*100, decimal)
                RMSD_E = np.round(np.sqrt(sum(dEsr)/len(dEsr)), decimal)
                MAE_E  = np.round(sum(dEabs)/len(dEabs), decimal)
                dico_Energie[f"H_SO + {O_k}"]['RMSD'] = RMSD_E
                dico_Energie[f"H_SO + {O_k}"]['MAE'] = MAE_E                                                    
                dico_Energie[f"H_SO + {O_k}"]['MAEER'] = Guihery_Nathalie_Metric 
                O_k += ' + O'
            if print_large:
                print(f'===================Rest: H_ab - H_mod===================')
                print(f'===H_mod = H_SO + H_CF(O2 + O4 + O6 + O8 + O10 + O12)===')
                H_rest = H_ab - H_mod
                pretty_matrix(H_rest, mag.basis_jmj, mag, decimal)
                print('\n') 
                
            df=pd.DataFrame(dico_Energie)
            print('\n------Reconstruction of the Energy Spectrum------',end='\n\n')
            print(tabulate(df, headers='keys',tablefmt='pipe'),end='\n\n')
            print('RMSD: Root Mean Square Deviation on energies in cm**-1' )
            print('MAE: Mean Absolute Error on energies in cm**-1')
            print('MAEER: Mean Absolute Error on Energies Relative to the ab initio spectral width in %' ,end='\n\n')
    return



def op_ls_jmj_iso(jmj_ket):
    """
    Compute the spin-orbit coupling terms for a given |JMJ > 
    """
    def lz_sz(term):
        """Compute the Lz * Sz term."""
        return term.op_Sz().op_Jz()

    def lp_sm(term):
        """Compute the 1/2 (L+ * S-) term."""
        return term.op_Sm().op_Jp().times_cst(sp.Rational(1, 2))
        #a = term.op_Sm().op_Jp()
        #a.coef *= sp.Rational(1, 2)
        #return a

    def lm_sp(term):
        """Compute the 1/2 (L- * S+) term."""
        return term.op_Sp().op_Jm().times_cst(sp.Rational(1, 2))
        #a = term.op_Sp().op_Jm()
        #a.coef *= sp.Rational(1, 2)
        #return a

    # Generate all terms and filter out zero-coefficient results in a single step
    res = []

    for op in (lz_sz, lp_sm, lm_sp):
        for t in jmj_ket.lmlsms_expansion:
            term = op(t)
            if term.coef != 0:
                res.append(term)
 
#    res = [
#        term for op in (lz_sz, lp_sm, lm_sp) 
#        for term in (op(t) for t in jmj_ket.lmlsms_expansion) 
#        if term.coef != 0
#    ]
    return res 

def matrix_ls_jmj_iso(mag_cntr):
    '''
    Compute the SO matrix to do a ITO method between this and the numerical one
    '''
    basis = mag_cntr.basis_jmj
    # Initialize the SO matrix
    SO_mat = sp.zeros(len(basis), len(basis))

    # Precompute op_soc for each bket in the basis
    soc_dict = {idx: op_ls_jmj_iso(basis[idx]) for idx in range(len(basis))}

    # Loop over all pairs of basis elements
    for idx in range(len(basis)):
        for jdx in range(len(basis)):
            # Get the corresponding SOC expansion for the ket
            soc = soc_dict[jdx]  # Precomputed SOC terms for bket
            SO_mat[idx, jdx] = inner_prod(basis[idx], soc)  # Compute inner product

    # Return the matrix (scaled by zeta if necessary)
    return np.array(SO_mat, dtype=np.float64)

def extract_soc_iso(mag_center=None, heff=None):
    '''
    extract the lambda/zeta constant by the ITO method
    '''
    soc_extracted = {}

    # Convert heff to numeric, this heff is in the JMJ basis
    heff = np.array(heff).astype(np.complex128)
    # Rotate to JMJ coupled basis
    bas_lmlsms = mag_center.basis_lmlsms
    bas_jmj = mag_center.basis_jmj

    ls_mod = matrix_ls_jmj_iso(mag_center)
    num = np.trace(heff @ ls_mod)
    denom = np.trace(ls_mod @ ls_mod)
    Lambda = np.real(num / denom)
    if mag_center.nel < int(2*mag_center.l+1):
        Zeta=Lambda*2*float(mag_center.S)
    elif mag_center.nel > int(2*mag_center.l+1):
        Zeta=-Lambda*2*float(mag_center.S)
    return Zeta, Lambda

def op_ls_jmj_aniso(jmj_ket,xyz):
    """
    Compute the xyz spin-orbit coupling terms for a given |JMJ > 
    """
    def p_lp_sp(term):
        """Compute the 1/4 (L+ * S+ ) term."""
        return term.op_Sp().op_Jp().times_cst(sp.Rational(1, 4)) 

    def p_lm_sm(term):
        """Compute the 1/4 (L- * S- ) term."""
        return term.op_Sm().op_Jm().times_cst(sp.Rational(1, 4)) 

    def m_lp_sp(term):
        """Compute the -1/4 (L+ * S+ ) term."""
        return term.op_Sp().op_Jp().times_cst(sp.Rational(-1, 4)) 

    def m_lm_sm(term):
        """Compute the -1/4 (L- * S- ) term."""
        return term.op_Sm().op_Jm().times_cst(sp.Rational(-1, 4)) 

    def p_lp_sm(term):
        """Compute the 1/4 (L+ * S-) term."""
        return term.op_Sm().op_Jp().times_cst(sp.Rational(1, 4)) 

    def p_lm_sp(term):
        """Compute the 1/4 (L- * S+ ) term."""
        return term.op_Sp().op_Jm().times_cst(sp.Rational(1, 4)) 

    def lz_sz(term):
        """Compute the Lz * Sz term."""
        return term.op_Sz().op_Jz()


    # Generate all terms and filter out zero-coefficient results in a single step
    res = []

    if xyz=='x':
        for op in (p_lp_sp, p_lp_sm, p_lm_sp, p_lm_sm):
            for t in jmj_ket.lmlsms_expansion:
                term = op(t)
                if term.coef != 0:
                    res.append(term)
    if xyz=='y':
        for op in (m_lp_sp, p_lp_sm, p_lm_sp, m_lm_sm):
            for t in jmj_ket.lmlsms_expansion:
                term = op(t)
                if term.coef != 0:
                    res.append(term)
    if xyz=='z':
            for t in jmj_ket.lmlsms_expansion:
                term = lz_sz(t)
                if term.coef != 0:
                    res.append(term)
 
    return res 

def matrix_ls_jmj_aniso(mag_cntr,xyz):
    '''
    Compute the SO matrix to do a ITO method between this and the numerical one
    '''
    basis = mag_cntr.basis_jmj
    # Initialize the SO matrix
    SO_mat = sp.zeros(len(basis), len(basis))

    # Precompute op_soc for each bket in the basis
    soc_dict = {idx: op_ls_jmj_aniso(basis[idx],xyz) for idx in range(len(basis))}

    # Loop over all pairs of basis elements
    for idx in range(len(basis)):
        for jdx in range(len(basis)):
            # Get the corresponding SOC expansion for the ket
            soc = soc_dict[jdx]  # Precomputed SOC terms for bket
            SO_mat[idx, jdx] = inner_prod(basis[idx], soc)  # Compute inner product

    # Return the matrix (scaled by zeta if necessary)
    return np.array(SO_mat, dtype=np.float64)

def extract_soc_aniso(mag_center=None, heff=None):
    '''
    extract the lambda/zeta constant by the ITO method
    '''
    soc_extracted = {}

    # Convert heff to numeric, this heff is in the JMJ basis
    heff = np.array(heff).astype(np.complex128)
    # Rotate to JMJ coupled basis
    bas_lmlsms = mag_center.basis_lmlsms
    bas_jmj = mag_center.basis_jmj

    list_Lambda=[]
    list_Zeta=[]
    for xyz in ['x','y','z']:
        ls_mod = matrix_ls_jmj_aniso(mag_center, xyz)
        num = np.trace(heff @ ls_mod)
        denom = np.trace(ls_mod @ ls_mod)
        Lambda = np.real(num / denom)
        if mag_center.nel < int(2*mag_center.l+1):
            Zeta=Lambda*2*float(mag_center.S)
        elif mag_center.nel > int(2*mag_center.l+1):
            Zeta=-Lambda*2*float(mag_center.S)
        list_Lambda.append(Lambda)
        list_Zeta.append(Zeta)
    return list_Zeta, list_Lambda



def cf_hamiltonian(mag_center=None, bkq_values=None, bkq_group=None, ls_constants=None):
    '''
    reconstruction of the H_CF to determine the H_SO matrix 
    H_ab = H_SOC + H_CF <=> H_SOC = H_ab - H_CF
    '''
    bas_lmlsms = mag_center.basis_lmlsms
    bas_lmls = mag_center.basis_lmls
    bas_jmj = mag_center.basis_jmj
    dim = len(bas_lmlsms)
#    hcf_lmlsms = sp.Matrix.zeros(len(bas_lmlsms), len(bas_lmlsms))  # initialize the model SR/SO Hamiltonians
    hcf_lmls, hcf_lmlsms = sp.Matrix.zeros(len(bas_lmls), len(bas_lmls)), sp.Matrix.zeros(len(bas_lmlsms), len(bas_lmlsms))  # initialize the model SR/SO Hamiltonians
    ls_matrix_jmj = np.zeros((len(bas_jmj), len(bas_jmj)))
    # let us fill the CF part
    for key, val in Bkq.items():
        if key in bkq_group:
            if ls_constants == False:
                hcf_lmls += matrix_cf(mag_center, bas_lmls, key) * (bkq_values[val])
            else:
                hcf_lmlsms += matrix_cf(mag_center, bas_lmlsms, key) * (bkq_values[val])

    # deal with SOC part
#    if ls_constants != None:
#        ls_matrix_jmj = matrix_ls_jmj(mag_center) * ls_constants
    
#    hcf_lmls, hcf_jmj = np.array(hcf_lmls, dtype=np.complex128), np.array(hcf_jmj, dtype=np.complex128)
    if ls_constants == False:
        hcf_lmls = np.array(hcf_lmls.evalf() , dtype=np.complex128)
        return hcf_lmls
    else :
        ucg = cg_lml_jmj(bas_lmlsms, bas_jmj)
        hcf_jmj = ucg @ hcf_lmlsms @ ucg.T
        hcf_jmj = np.array(hcf_jmj, dtype=np.complex128)
        return hcf_jmj

######################
# GENERAL ERROR
####################

def general_error_inp_molcas(mag_center,end=True):
    print('exemple OpenMolcas_v25.06 input:')
    if mag_center.l==3: OA='f'
    if mag_center.l==2: OA='d'
    if mag_center.l==1: OA='p'
    n_cas = (2*mag_center.S+1)/(2*mag_center.l+2)
    n_cas *= math.comb(2*mag_center.l+2, int(mag_center.nel/2-mag_center.S))
    n_cas *= math.comb(2*mag_center.l+2, int(mag_center.nel/2+mag_center.S+1))
    print(f"""

>>> EXPORT MOLCAS_PRINT = 3
&SEWARD  
.        
.        
.        
End of input
         
&SCF     
.        
.        
.        
End of input
         
*Compute CAS({mag_center.nel},{2*mag_center.l+1})SCF
&RASSCF  
LUMORB   
Spin={2*mag_center.S+1}
Symmetry=1
nActEl={mag_center.nel} 0 0
Inactive= -> To define by yourself 
RAS2 = {2*mag_center.l+1}
CIRoots={n_cas} {n_cas} 1
Iter=200 100
ORBListing=all
ORBAppear=compact
PRWF=0   
End of input
         
*Orbital localisation of the {2*mag_center.l+1} 
*active orbital, using Cholesky methode
*(can work with other methode, PAO, PM,...)
&LOCALISATION
NFrozen= -> To define by yourself 
NORbitals= {2*mag_center.l+1}
CHOLesky 
end of input
         
*Compute CAS({mag_center.nel},{2*mag_center.l+1})CI
*with the localised Orb file to get each CFS in the single OA basis 
&RASSCF  
LUMORB   
cionly   
Spin={2*mag_center.S+1}
Symmetry=1
nActEl={mag_center.nel} 0 0
Inactive= -> To define by yourself
RAS2 = {2*mag_center.l+1} 
CIRoots={n_cas} {n_cas} 1
Iter=200 100                                                                                
ORBListing=all
ORBAppear=compact
PRWF=0   
End of input

*if you want a CASPT2 level
&CASPT2 
Multistate=all
CONVergence=1.0d-08
PRWF=0
NoMult
End of input                        
 
>>COPY $Project.JobMix JOB001

*Compute SOCI procedure on the all scalar roots
>>> EXPORT MOLCAS_PRINT = 5
&RASSI
EJOB -> if you want a CASPT2 level
NROF=1 all
SPINORBIT
THRS=0
SOCOupling=0
end of input

""") 
    if end:
        print('EXIT Program')
        exit()
    return
 
def general_error_inp_orca(mag_center,end=True):
    print('exemple ORCA_v6.1 input for the %casscf block:')                                                                         
    if mag_center.l==3: OA='f'
    if mag_center.l==2: OA='d'
    if mag_center.l==1: OA='p'
    n_cas = (2*mag_center.S+1)/(2*mag_center.l+2)
    n_cas *= math.comb(2*mag_center.l+2, int(mag_center.nel/2-mag_center.S))
    n_cas *= math.comb(2*mag_center.l+2, int(mag_center.nel/2+mag_center.S+1))
    print(f"""

! largeprint
%casscf
nel     {mag_center.nel}
norb    {2*mag_center.l+1}
mult    {2*mag_center.S+1}
nroots  {n_cas}
printwf det
actorbs {OA}orbs
ci
    TPrintWF   0
end
rel
    dosoc true
    tprint 0.0
end
end

""")
 
    if end:
        print('EXIT Program')
        exit()
    return

