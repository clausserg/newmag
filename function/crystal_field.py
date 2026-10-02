import sympy as sp
import numpy as np
import re
from copy import deepcopy
from tabulate import tabulate
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from function.class_function import *
#from function.helper_functions import *
from function.helper_build_hamiltonian import *

def cfps_calc(Heff, mag, out, decimal):
    '''
    main function of the ITO method to extract the Bkq parameters
    
    INPUT
    Heff : dict[level][sr or so][array heff], all effectif hamiltonienne
    mag : class of functions possessing all the properties characteristic of the metallic center 
    out : main orca/molcas output to extract single_aniso block
    decimal :  number of decimal printed in the NewMag out 

    OUTPUT
    cfp : dict[level][Bkq][value], dictionary with all bkq values as a function of the calculation level
    '''
    cfp={}
    for level in Heff.keys():
        if Heff[level] is not None:
            if level.split('-')[-1] == 'sr':
                cfp[level] = extract_cf_prms(mag_center=mag, heff=Heff[level], soc=False)
            elif level.split('-')[-1] == 'so':
                cfp[level] = extract_cf_prms(mag_center=mag, heff=Heff[level],soc=True)
            elif level.split('-')[-1] == 'ailft':
                cfp[level] = extract_cf_prms(mag_center=MagneticCenter(l=mag.l, nel=1),heff=Heff[level],soc=False)
     
    print("All B_k^q values are expressed in cm**-1 (using Steven's non normalized operators)")
    print("All S_k^q values are expressed in cm**-1 (after correcting the Steven's parameters to ensure normalization)")
    print('(x_aniso: From the x^th iteration of the single_aniso program in the output)')
    print()
    
    #print the Bkq DataFrame with tabulate python package
    print(tabulate(print_cf_params(cfp, cfp.keys(), out, mag, decimal), headers='keys',tablefmt='pipe'))
    return cfp

def print_cf_params(param, list_calc_level, outpute, mag, decimal):
    '''
    print the bkq value as a DataFrame and calculation of S_k / S^q / S values
    
    INPUT
    param : dict[level][Bkq][value], dictionary with all bkq values as a function of the calculation level
    list_calc_level : list of calculation level 
    outpute : main orca/molcas output to extract single_aniso block
    mag : class of functions possessing all the properties characteristic of the metallic center 
    decimal : number of decimal printed in the NewMag out 
    
    OUTPUT
    df : DataFrame of the Bkq
    '''
    # parameter groups                                                                                                                                           
    b2 = [B22m, B21m, B20, B21, B22]
    b4 = [B44m, B43m, B42m, B41m, B40, B41, B42, B43, B44]
    b6 = [B66m, B65m, B64m, B63m, B62m, B61m, B60, B61, B62, B63, B64, B65, B66]
    b8 = [B88m, B87m, B86m, B85m, B84m, B83m, B82m, B81m, B80, B81, B82, B83, B84, B85, B86, B87, B88]
    b10 = [B1010m, B109m, B108m, B107m, B106m, B105m, B104m, B103m, B102m, B101m, B100, B101, B102, B103, B104, B105, B106, B107, B108, B109, B1010]
    b12 = [B1212m, B1211m, B1210m, B129m, B128m, B127m, B126m, B125m, B124m, B123m, B122m, B121m, B120, B121, B122, B123, B124, B125, B126, B127, B128, B129, B1210, B1211, B1212]
     
    # print parameters
    dict_paramerters={}
    list_calc_level=sorted(list_calc_level)
    for level in list_calc_level:
        dict_paramerters[level]={}
        for group in (b2, b4, b6, b8, b10, b12):
            for par in group:
                dict_paramerters[level][par]=param[level][par]
    aniso_J=False
    aniso_L=False
    iter_J=1
    iter_L=1
    with open(outpute, "r", encoding="utf-8") as f:
        for line in f:       
            line = line.rstrip()
            header_line_L = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L ="
            header_line_J = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET J ="
            if header_line_J in line:
                J_value=str(line.split()[-1])
                dict_paramerters[f'{iter_J}_aniso_J={J_value}']={}
                for group in (b2, b4, b6, b8, b10, b12):
                    for par in group:
                        dict_paramerters[f'{iter_J}_aniso_J={J_value}'][par]=0
                iter_J+=1
                aniso_J=True
            elif header_line_L in line:
                L_value=str(line.split()[-1])
                dict_paramerters[f'{iter_L}_aniso_L={L_value}']={}
                for group in (b2, b4, b6, b8, b10, b12):
                    for par in group:
                        dict_paramerters[f'{iter_L}_aniso_L={L_value}'][par]=0
                iter_L+=1
                aniso_L=True
     
    if aniso_J != aniso_L:
        i=0
        for iteration in range(1,iter_J):
            i+=1
            L_value='None'
            dict_paramerters[f'{iteration}_aniso_L={L_value}']={}
            for group in (b2, b4, b6, b8, b10, b12):
                for par in group:    
                    dict_paramerters[f'{iteration}_aniso_L={L_value}'][par]=0
        for iteration in range(1,iter_L):
            i+=1
            J_value='None'
            dict_paramerters[f'{iteration}_aniso_J={J_value}']={} 
            for group in (b2, b4, b6, b8, b10, b12):
                for par in group:    
                    dict_paramerters[rf'{iteration}_aniso_J={J_value}'][par]=0


    if iter_J != iter_L:
        inter = i
    else :
        inter = iter_J
    for iteration in range(1,inter):
        dict_paramerters[f'{iteration}_aniso_L={L_value}'],dict_paramerters[f'{iteration}_aniso_J={J_value}']=aniso_params(outpute,mag,iteration)    

    dict_Normalization_Steven={2:{2:np.sqrt(3), 1:2*np.sqrt(3), 0:1}, 
                                  4:{4:np.sqrt(35), 3:2*np.sqrt(70), 2:2*np.sqrt(5), 1:2*np.sqrt(10), 0:1}, 
                                  6:{6:np.sqrt(231/2), 5:3*np.sqrt(154), 4:3*np.sqrt(7), 3:np.sqrt(210), 2:np.sqrt(105/2), 1:2*np.sqrt(21), 0:1}}

    Normalized_Steven_dict_paramerters = deepcopy(dict_paramerters)
    for level in Normalized_Steven_dict_paramerters.keys():    
        for k, group in zip([2,4,6], (b2,b4,b6)):
            for par in group:
                q=int(list(str(par))[-1])
                Normalized_Steven_dict_paramerters[level][par] = Normalized_Steven_dict_paramerters[level][par]/(dict_Normalization_Steven[k][q])

    for level in dict_paramerters.keys():
        dict_paramerters[level]['---S_k---']='---'
        for k, group in zip([2,4,6], (b2, b4, b6)):
            sum_bkq=0
            for par in group:
                q=int(list(str(par))[-1])
                sum_bkq += Normalized_Steven_dict_paramerters[level][par]**2
            dict_paramerters[level][f'S_{k}']=np.sqrt(1/(2*k+1)*sum_bkq)

    for level in dict_paramerters.keys():
        dict_paramerters[level]['---S^q---']='---'
        for q in range(0,7):
            sum_bkq=0
            for k, group in zip([2,4,6], (b2, b4, b6)):
                for par in group:  
                    if int(list(str(par))[-1])==q:
                        sum_bkq += 1/(2*k+1)*Normalized_Steven_dict_paramerters[level][par]**2
            dict_paramerters[level][f'S^{q}']=np.sqrt(sum_bkq)

    for level in dict_paramerters.keys():
        dict_paramerters[level]['---S---']='---'
        for k in [2,4,6]:
            dict_paramerters[level][f'S form S_k']=np.sqrt(1/3*(dict_paramerters[level][f'S_{2}']**2
                                                               +dict_paramerters[level][f'S_{4}']**2
                                                               +dict_paramerters[level][f'S_{6}']**2))
        for q in range(0,7):
            dict_paramerters[level][f'S form S_q']=np.sqrt(1/3*(dict_paramerters[level][f'S^{0}']**2 
                                                               +dict_paramerters[level][f'S^{1}']**2
                                                               +dict_paramerters[level][f'S^{2}']**2
                                                               +dict_paramerters[level][f'S^{3}']**2
                                                               +dict_paramerters[level][f'S^{4}']**2
                                                               +dict_paramerters[level][f'S^{5}']**2
                                                               +dict_paramerters[level][f'S^{6}']**2))
            
    
    df=pd.DataFrame(dict_paramerters)                     
    df = df.map(lambda x: f"{x.real:6.{decimal}f}" if not isinstance(x, str) else x)
    return df


def matrix_cf(magcent, basis, sto):
    '''
    building the Okq matrix to apply the ITO method
    for the extraction of the Bkq
    
    INPUT 
    magcent :class of functions possessing all the properties characteristic of the metallic center
    basis : specific basis define with magcent
    sto : Steven operator 


    OUTPUT
    Okq_mat : Okq matrix 
    '''
    Okq_mat = sp.zeros(len(basis), len(basis))
    for idx in range(len(basis)):
        for jdx in range(len(basis)):
            aket, bket = basis[idx], basis[jdx]
            expansion = Okq[sto](bket)
            for term in expansion:
                if sto[0] == 2 and isinstance(bket, JMJ):
                    term.coef *= magcent.factor_abg(bket)[0]
                if sto[0] == 2 and isinstance(bket, LMLSMS):
                    term.coef *= magcent.factor_abg()[0]
                if sto[0] == 4 and isinstance(bket, JMJ):
                    term.coef *= magcent.factor_abg(bket)[1]
                if sto[0] == 4 and isinstance(bket, LMLSMS):
                    term.coef *= magcent.factor_abg()[1]
                if sto[0] == 6 and isinstance(bket, JMJ):
                    term.coef *= magcent.factor_abg(bket)[2]
                if sto[0] == 6 and isinstance(bket, LMLSMS):
                    term.coef *= magcent.factor_abg()[2]
                if sto[0] not in [2,4,6]: 
                    term.coef *= 1    
                if term.Ket == aket.Ket:
                    Okq_mat[idx, jdx] += term.coef
    return Okq_mat

def extract_cf_prms(mag_center=None, heff=None, soc=False):
    '''
    ITO procedure to extract the Bkq parameters

    INPUT
    mag_center : class of functions possessing all the properties characteristic of the metallic center 
    heff : numerical effectif hamiltonienne (ab-initio)
    soc : bool True/False

    OUTPUT
    cfp_extracted : dict[bkq][value], all bkq values form the heff
    '''
    cfp_extracted = {}

    # Convert heff to numeric
    heff = np.array(heff).astype(np.complex128)
    if soc:
        basis = mag_center.basis_lmlsms
        bas_jmj = mag_center.basis_jmj
        ucg = cg_lml_jmj(basis, bas_jmj)
        ucg_conj_t = ucg.conj().T
        heff_eff = ucg_conj_t @ heff @ ucg  # rotated heff
    else:
        basis = mag_center.basis_lmls
        heff_eff = heff  # unrotated

    # Shared matrix_cf cache
    cf_cache = {}

    def process_key_val(key_val):
        key, val = key_val
        if key not in cf_cache:
            cf_cache[key] = np.array(matrix_cf(mag_center, basis, key)).astype(np.complex128)
        cf_mat = cf_cache[key]
        num = np.trace(heff_eff @ cf_mat)
        denom = np.trace(cf_mat @ cf_mat)
        if denom!=np.complex128(0):
            cfp_value = np.real(num / denom)
        else:cfp_value=0
        return val, cfp_value

    # Parallel execution
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_key_val, kv) for kv in Bkq.items()]
        for future in as_completed(futures):
            val, result = future.result()
            cfp_extracted[val] = result

    return cfp_extracted

def aniso_params(file_path: str, mag, iteration):

    '''
    reading the bkq form sinlge_aniso program form molcas/orca out

    INPUT
    file_path : name of main output molcas/orca
    mag : class of functions possessing all the properties characteristic of the metallic center 
    iteration : number of out sinlge_aniso in the main output molcas/orca   

    OUTPUT
    dict_single_aniso_L : dict[bkq][value] for L representation
    dict_single_aniso_J : dict[bkq][value] for J_ground_state representation 
    '''
    b2 = [B22m, B21m, B20, B21, B22]
    b4 = [B44m, B43m, B42m, B41m, B40, B41, B42, B43, B44]     
    b6 = [B66m, B65m, B64m, B63m, B62m, B61m, B60, B61, B62, B63, B64, B65, B66]
    b8 = [B88m, B87m, B86m, B85m, B84m, B83m, B82m, B81m, B80, B81, B82, B83, B84, B85, B86, B87, B88] 
    b10 = [B1010m, B109m, B108m, B107m, B106m, B105m, B104m, B103m, B102m, B101m, B100, B101, B102, B103, B104, B105, B106, B107, B108, B109, B1010]
    b12 = [B1212m, B1211m, B1210m, B129m, B128m, B127m, B126m, B125m, B124m, B123m, B122m, B121m, B120, B121, B122, B123, B124, B125, B126, B127, B128, B129, B1210, B1211, B1212]

    dict_single_aniso_L={}
    dict_single_aniso_J={}
    for group in (b2, b4, b6, b8, b10, b12):      
        for par in group:           
            dict_single_aniso_L[par]=0
            dict_single_aniso_J[par]=0

    header_line_L = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L ="
    header_line_J = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET J ="
    table_header = "k |  q  |    (K)^2    |         B(k,q)        |"
    stop_marker = "********************************************************************************"

    a,b,g=mag.Stevens_coeff()
    #except TypeError :
    #    print('Error Single_ansio reading, pass')
    #    return dict_single_aniso_L, dict_single_aniso_J
    Bkq = {}
    in_section = False
    in_table = False
    count=1
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()

            # Step 1: wait until the main header is found
            if not in_section:
                if header_line_J in line:
                    if count==iteration:
                        in_section = True
                    else: count+=1
                continue

            # Step 2: wait for the right table header
            if in_section and not in_table:
                if table_header in line:
                    in_table = True
                continue

            # Step 3: parse the table lines until stop marker
            if in_table:
                # Stop if we reach the stop marker
                if stop_marker in line:
                    break

                # Skip separators
                if line.startswith("----") or line.startswith("------------------------------------------------|"):
                    continue

                # Match table entries
                match = re.match(
                    r"^\s*(\d+)\s*\|\s*(-?\d+)\s*\|\s*[\d.E+-]+\s*\|\s*([-\d.E+]+)\s*\|",
                    line,
                )
                if match:
                    k = int(match.group(1))
                    q = int(match.group(2))
                    B_value = np.float64(match.group(3))
                    if k==2  and a!=0:
                        if q==-2:
                            dict_single_aniso_J[B22m]=B_value*a
                        elif q==-1:
                            dict_single_aniso_J[B21m]=B_value*a
                        elif q==0:
                            dict_single_aniso_J[B20]=B_value*a
                        elif q==1:
                            dict_single_aniso_J[B21]=B_value*a
                        elif q==2:
                            dict_single_aniso_J[B22]=B_value*a
                    elif k==4  and b!=0:
                        if q==-4:
                            dict_single_aniso_J[B44m]=B_value*b
                        elif q==-3:
                            dict_single_aniso_J[B43m]=B_value*b
                        elif q==-2:
                            dict_single_aniso_J[B42m]=B_value*b
                        elif q==-1:
                            dict_single_aniso_J[B41m]=B_value*b
                        elif q==0:
                            dict_single_aniso_J[B40]=B_value*b
                        elif q==1:
                            dict_single_aniso_J[B41]=B_value*b
                        elif q==2:
                            dict_single_aniso_J[B42]=B_value*b
                        elif q==3:
                            dict_single_aniso_J[B43]=B_value*b
                        elif q==4:
                            dict_single_aniso_J[B44]=B_value*b
                    elif k==6  and g!=0:
                        if q==-6:
                            dict_single_aniso_J[B66m]=B_value*g
                        elif q==-5:
                            dict_single_aniso_J[B65m]=B_value*g
                        elif q==-4:
                            dict_single_aniso_J[B64m]=B_value*g
                        elif q==-3:
                            dict_single_aniso_J[B63m]=B_value*g
                        elif q==-2:
                            dict_single_aniso_J[B62m]=B_value*g
                        elif q==-1:
                            dict_single_aniso_J[B61m]=B_value*g
                        elif q==0:
                            dict_single_aniso_J[B60]=B_value*g
                        elif q==1:
                            dict_single_aniso_J[B61]=B_value*g
                        elif q==2:
                            dict_single_aniso_J[B62]=B_value*g
                        elif q==3:
                            dict_single_aniso_J[B63]=B_value*g
                        elif q==4:
                            dict_single_aniso_J[B64]=B_value*g
                        elif q==5:
                            dict_single_aniso_J[B65]=B_value*g
                        elif q==6:
                            dict_single_aniso_J[B66]=B_value*g

                    elif k==8:
                        if q==-8:
                            dict_single_aniso_J[B88m]=B_value
                        elif q==-7: 
                            dict_single_aniso_J[B87m]=B_value
                        elif q==-6:
                            dict_single_aniso_J[B86m]=B_value
                        elif q==-5:
                            dict_single_aniso_J[B85m]=B_value
                        elif q==-4:
                            dict_single_aniso_J[B84m]=B_value
                        elif q==-3:
                            dict_single_aniso_J[B83m]=B_value
                        elif q==-2:
                            dict_single_aniso_J[B82m]=B_value
                        elif q==-1:
                            dict_single_aniso_J[B81m]=B_value
                        elif q==0:
                            dict_single_aniso_J[B80]=B_value
                        elif q==1:
                            dict_single_aniso_J[B81]=B_value
                        elif q==2:
                            dict_single_aniso_J[B82]=B_value
                        elif q==3:
                            dict_single_aniso_J[B83]=B_value
                        elif q==4:
                            dict_single_aniso_J[B84]=B_value
                        elif q==5:
                            dict_single_aniso_J[B85]=B_value
                        elif q==6:
                            dict_single_aniso_J[B86]=B_value
                        elif q==7:
                            dict_single_aniso_J[B87]=B_value
                        elif q==8:
                            dict_single_aniso_J[B88]=B_value
                    elif k==10:
                        if q==-10:
                            dict_single_aniso_J[B1010m]=B_value
                        elif q==-9: 
                            dict_single_aniso_J[B109m]=B_value
                        elif q==-8: 
                            dict_single_aniso_J[B108m]=B_value
                        elif q==-7:
                            dict_single_aniso_J[B107m]=B_value
                        elif q==-6:  
                            dict_single_aniso_J[B106m]=B_value
                        elif q==-5:
                            dict_single_aniso_J[B105m]=B_value
                        elif q==-4:
                            dict_single_aniso_J[B104m]=B_value
                        elif q==-3:
                            dict_single_aniso_J[B103m]=B_value
                        elif q==-2:
                            dict_single_aniso_J[B102m]=B_value
                        elif q==-1:
                            dict_single_aniso_J[B101m]=B_value
                        elif q==0:
                            dict_single_aniso_J[B100]=B_value
                        elif q==1:
                            dict_single_aniso_J[B101]=B_value
                        elif q==2:
                            dict_single_aniso_J[B102]=B_value
                        elif q==3:
                            dict_single_aniso_J[B103]=B_value
                        elif q==4:
                            dict_single_aniso_J[B104]=B_value
                        elif q==5:
                            dict_single_aniso_J[B105]=B_value
                        elif q==6:
                            dict_single_aniso_J[B106]=B_value
                        elif q==7:
                            dict_single_aniso_J[B107]=B_value
                        elif q==8:
                            dict_single_aniso_J[B108]=B_value
                        elif q==9:
                            dict_single_aniso_J[B109]=B_value
                        elif q==10:
                            dict_single_aniso_J[B1010]=B_value

                    elif k==12:
                        if q==-12:
                            dict_single_aniso_J[B1212m]=B_value
                        elif q==-11: 
                            dict_single_aniso_J[B1211m]=B_value
                        elif q==-10: 
                            dict_single_aniso_J[B1010m]=B_value
                        elif q==-9: 
                            dict_single_aniso_J[B129m]=B_value
                        elif q==-8: 
                            dict_single_aniso_J[B128m]=B_value
                        elif q==-7:
                            dict_single_aniso_J[B127m]=B_value
                        elif q==-6:  
                            dict_single_aniso_J[B126m]=B_value
                        elif q==-5:
                            dict_single_aniso_J[B125m]=B_value
                        elif q==-4:
                            dict_single_aniso_J[B124m]=B_value
                        elif q==-3:
                            dict_single_aniso_J[B123m]=B_value
                        elif q==-2:
                            dict_single_aniso_J[B122m]=B_value
                        elif q==-1:
                            dict_single_aniso_J[B121m]=B_value
                        elif q==0:
                            dict_single_aniso_J[B120]=B_value
                        elif q==1:
                            dict_single_aniso_J[B121]=B_value
                        elif q==2:
                            dict_single_aniso_J[B122]=B_value
                        elif q==3:
                            dict_single_aniso_J[B123]=B_value
                        elif q==4:
                            dict_single_aniso_J[B124]=B_value
                        elif q==5:
                            dict_single_aniso_J[B125]=B_value
                        elif q==6:
                            dict_single_aniso_J[B126]=B_value
                        elif q==7:
                            dict_single_aniso_J[B127]=B_value
                        elif q==8:
                            dict_single_aniso_J[B128]=B_value
                        elif q==9:
                            dict_single_aniso_J[B129]=B_value
                        elif q==10:
                            dict_single_aniso_J[B1210]=B_value
                        elif q==11:
                            dict_single_aniso_J[B1211]=B_value
                        elif q==12:
                            dict_single_aniso_J[B1212]=B_value




    Bkq = {}
    in_section = False
    in_table = False
    count=1 
    a, b, g = mag.factor_abg()
    a=np.float64(a)
    b=np.float64(b)
    g=np.float64(g)
    with open(file_path, "r", encoding="utf-8") as f:
        i=0
        for line in f:
            line = line.rstrip()

            # Step 1: wait until the main header is found
            if not in_section:
                if header_line_L in line:
                    if count==iteration:
                        in_section = True
                    else: count+=1
                continue

            # Step 2: wait for the right table header
            if in_section and not in_table:
                if table_header in line:
                    in_table = True
                continue

            # Step 3: parse the table lines until stop marker
            if in_table:
                # Stop if we reach the stop marker
                if stop_marker in line:
                    break

                # Skip separators
                if line.startswith("----") or line.startswith("------------------------------------------------|"):
                    continue

                # Match table entries
                match = re.match(
                    r"^\s*(\d+)\s*\|\s*(-?\d+)\s*\|\s*[\d.E+-]+\s*\|\s*([-\d.E+]+)\s*\|",
                    line,
                )
                if match:
                    k = int(match.group(1))
                    q = int(match.group(2))
                    B_value = np.float64(match.group(3))
                    if k==2  and a!=0:
                        if q==-2:
                            dict_single_aniso_L[B22m]=B_value/a
                        elif q==-1:
                            dict_single_aniso_L[B21m]=B_value/a
                        elif q==0:
                            dict_single_aniso_L[B20]=B_value/a
                        elif q==1:
                            dict_single_aniso_L[B21]=B_value/a
                        elif q==2:
                            dict_single_aniso_L[B22]=B_value/a
                    elif k==4  and b!=0:
                        if q==-4:
                            dict_single_aniso_L[B44m]=B_value/b
                        elif q==-3:
                            dict_single_aniso_L[B43m]=B_value/b
                        elif q==-2:
                            dict_single_aniso_L[B42m]=B_value/b
                        elif q==-1:
                            dict_single_aniso_L[B41m]=B_value/b
                        elif q==0:
                            dict_single_aniso_L[B40]=B_value/b
                        elif q==1:
                            dict_single_aniso_L[B41]=B_value/b
                        elif q==2:
                            dict_single_aniso_L[B42]=B_value/b
                        elif q==3:
                            dict_single_aniso_L[B43]=B_value/b
                        elif q==4:
                            dict_single_aniso_L[B44]=B_value/b
                    elif k==6  and g!=0:
                        if q==-6:
                            dict_single_aniso_L[B66m]=B_value/g
                        elif q==-5:
                            dict_single_aniso_L[B65m]=B_value/g
                        elif q==-4:
                            dict_single_aniso_L[B64m]=B_value/g
                        elif q==-3:
                            dict_single_aniso_L[B63m]=B_value/g
                        elif q==-2:
                            dict_single_aniso_L[B62m]=B_value/g
                        elif q==-1:
                            dict_single_aniso_L[B61m]=B_value/g
                        elif q==0:
                            dict_single_aniso_L[B60]=B_value/g
                        elif q==1:
                            dict_single_aniso_L[B61]=B_value/g
                        elif q==2:
                            dict_single_aniso_L[B62]=B_value/g
                        elif q==3:
                            dict_single_aniso_L[B63]=B_value/g
                        elif q==4:
                            dict_single_aniso_L[B64]=B_value/g
                        elif q==5:
                            dict_single_aniso_L[B65]=B_value/g
                        elif q==6:
                            dict_single_aniso_L[B66]=B_value/g
                    elif k==8:
                        if q==-8:
                            dict_single_aniso_L[B88m]=B_value
                        elif q==-7: 
                            dict_single_aniso_L[B87m]=B_value
                        elif q==-6:
                            dict_single_aniso_L[B86m]=B_value
                        elif q==-5:
                            dict_single_aniso_L[B85m]=B_value
                        elif q==-4:
                            dict_single_aniso_L[B84m]=B_value
                        elif q==-3:
                            dict_single_aniso_L[B83m]=B_value
                        elif q==-2:
                            dict_single_aniso_L[B82m]=B_value
                        elif q==-1:
                            dict_single_aniso_L[B81m]=B_value
                        elif q==0:
                            dict_single_aniso_L[B80]=B_value
                        elif q==1:
                            dict_single_aniso_L[B81]=B_value
                        elif q==2:
                            dict_single_aniso_L[B82]=B_value
                        elif q==3:
                            dict_single_aniso_L[B83]=B_value
                        elif q==4:
                            dict_single_aniso_L[B84]=B_value
                        elif q==5:
                            dict_single_aniso_L[B85]=B_value
                        elif q==6:
                            dict_single_aniso_L[B86]=B_value
                        elif q==7:
                            dict_single_aniso_L[B87]=B_value
                        elif q==8:
                            dict_single_aniso_L[B88]=B_value
                    elif k==10:
                        if q==-10:
                            dict_single_aniso_L[B1010m]=B_value
                        elif q==-9: 
                            dict_single_aniso_L[B109m]=B_value
                        elif q==-8: 
                            dict_single_aniso_L[B108m]=B_value
                        elif q==-7:
                            dict_single_aniso_L[B107m]=B_value
                        elif q==-6:  
                            dict_single_aniso_L[B106m]=B_value
                        elif q==-5:
                            dict_single_aniso_L[B105m]=B_value
                        elif q==-4:
                            dict_single_aniso_L[B104m]=B_value
                        elif q==-3:
                            dict_single_aniso_L[B103m]=B_value
                        elif q==-2:
                            dict_single_aniso_L[B102m]=B_value
                        elif q==-1:
                            dict_single_aniso_L[B101m]=B_value
                        elif q==0:
                            dict_single_aniso_L[B100]=B_value
                        elif q==1:
                            dict_single_aniso_L[B101]=B_value
                        elif q==2:
                            dict_single_aniso_L[B102]=B_value
                        elif q==3:
                            dict_single_aniso_L[B103]=B_value
                        elif q==4:
                            dict_single_aniso_L[B104]=B_value
                        elif q==5:
                            dict_single_aniso_L[B105]=B_value
                        elif q==6:
                            dict_single_aniso_L[B106]=B_value
                        elif q==7:
                            dict_single_aniso_L[B107]=B_value
                        elif q==8:
                            dict_single_aniso_L[B108]=B_value
                        elif q==9:
                            dict_single_aniso_L[B109]=B_value
                        elif q==10:
                            dict_single_aniso_L[B1010]=B_value
                    elif k==12:
                        if q==-12:
                            dict_single_aniso_J[B1212m]=B_value
                        elif q==-11: 
                            dict_single_aniso_J[B1211m]=B_value
                        elif q==-10: 
                            dict_single_aniso_J[B1010m]=B_value
                        elif q==-9: 
                            dict_single_aniso_J[B129m]=B_value
                        elif q==-8: 
                            dict_single_aniso_J[B128m]=B_value
                        elif q==-7:
                            dict_single_aniso_J[B127m]=B_value
                        elif q==-6:  
                            dict_single_aniso_J[B126m]=B_value
                        elif q==-5:
                            dict_single_aniso_J[B125m]=B_value
                        elif q==-4:
                            dict_single_aniso_J[B124m]=B_value
                        elif q==-3:
                            dict_single_aniso_J[B123m]=B_value
                        elif q==-2:
                            dict_single_aniso_J[B122m]=B_value
                        elif q==-1:
                            dict_single_aniso_J[B121m]=B_value
                        elif q==0:
                            dict_single_aniso_J[B120]=B_value
                        elif q==1:
                            dict_single_aniso_J[B121]=B_value
                        elif q==2:
                            dict_single_aniso_J[B122]=B_value
                        elif q==3:
                            dict_single_aniso_J[B123]=B_value
                        elif q==4:
                            dict_single_aniso_J[B124]=B_value
                        elif q==5:
                            dict_single_aniso_J[B125]=B_value
                        elif q==6:
                            dict_single_aniso_J[B126]=B_value
                        elif q==7:
                            dict_single_aniso_J[B127]=B_value
                        elif q==8:
                            dict_single_aniso_J[B128]=B_value
                        elif q==9:
                            dict_single_aniso_J[B129]=B_value
                        elif q==10:
                            dict_single_aniso_J[B1210]=B_value
                        elif q==11:
                            dict_single_aniso_J[B1211]=B_value
                        elif q==12:
                            dict_single_aniso_J[B1212]=B_value


    return dict_single_aniso_L, dict_single_aniso_J

def sto_O20(myket):
    mj = myket.Mj
    x = myket.J * (myket.J + 1)
    newket = deepcopy(myket)
    newket.coef = (3 * mj*mj - x)
    return [newket]

def sto_O21(myket):
    # 1/4 [JzJp + JzJm + JpJz + JmJz]
    a = myket.op_Jp().op_Jz().times_cst(sp.Rational(1,4))
    b = myket.op_Jm().op_Jz().times_cst(sp.Rational(1,4))
    c = myket.op_Jz().op_Jp().times_cst(sp.Rational(1,4))
    d = myket.op_Jz().op_Jm().times_cst(sp.Rational(1,4))
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O21m(myket):
    # -i/4 [JzJp - JzJm + JpJz - JmJz]
    a = myket.op_Jp().op_Jz().times_cst(-sp.I/4)
    b = myket.op_Jm().op_Jz().times_cst(sp.I/4)
    c = myket.op_Jz().op_Jp().times_cst(-sp.I/4)
    d = myket.op_Jz().op_Jm().times_cst(sp.I/4)
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O22(myket):
    # 1/2 [JpJp + JmJm]
    a = myket.op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O22m(myket):
    # -i/2 [JpJp - JmJm]
    a = myket.op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O40(myket):
    mj = myket.Mj
    x = myket.J * (myket.J + 1)
    new_coef = 35 * mj**4 - (30*x -25)*mj*mj + 3*x*x - 6*x
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]

def sto_O41(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**3 - (3*x+1)*mj) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**3 - (3*x+1)*mj) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**3 - (3*x+1)*mj) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O41m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**3 - (3*x+1)*mj) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(-1*cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**3 - (3*x+1)*mj) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**3 - (3*x+1)*mj) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O42(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**2 - x -5) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**2 - x -5) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**2 - x -5) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O42m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**2 - x -5) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**2 - x -5) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**2 - x -5) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O43(myket):
    # 1st term
    mj = myket.Mj
    cst = mj* sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * sp.Rational(1,4)
    mj = d.Mj
    d.coef *= mj * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O43m(myket):
    # 1st term
    mj = myket.Mj
    cst = mj * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * (-sp.I/4)
    mj = d.Mj
    d.coef *= mj * (sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O44(myket):
    # 1st term
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O44m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O60(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    new_coef = 231*mj**6 - (315*x-735)*mj**4 + (105*x*x-525*x+294)*mj*mj - 5*x**3 + 40*x*x - 60*x
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]

def sto_O61(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O61m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(-1*cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O62(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O62m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O63(myket):
    # 1st term
    mj, x = myket.Mj, myket.J*(myket.J+1)
    cst = (11*mj**3 - (3*x+59)*mj) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J*(c.J+1)
    c.coef *= (11*mj**3 - (3*x+59)*mj) * sp.Rational(1,4)
    mj, x = d.Mj, d.J*(d.J+1)
    d.coef *= (11*mj**3 - (3*x+59)*mj) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O63m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J*(myket.J+1)
    cst = (11*mj**3 - (3*x+59)*mj) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J*(c.J+1)
    c.coef *= (11*mj**3 - (3*x+59)*mj) * (-sp.I/4)
    mj, x = d.Mj, d.J*(d.J+1)
    d.coef *= (11*mj**3 - (3*x+59)*mj) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O64(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (11*mj**2 -x - 38) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (11*mj**2 -x - 38) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (11*mj**2 -x - 38) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O64m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (11*mj**2 -x - 38) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (11*mj**2 -x - 38) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (11*mj**2 -x - 38) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O65(myket):
    # 1st term
    mj = myket.Mj
    cst = mj* sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * sp.Rational(1,4)
    mj = d.Mj
    d.coef *= mj * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O65m(myket):
    # 1st term
    mj = myket.Mj
    cst = mj* (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * (-sp.I/4)
    mj = d.Mj
    d.coef *= mj * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O66(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O66m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O80(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    new_coef = sp.Rational(1,2)*(12870*mj**8 - 12012*(-9+2*x)*mj**6 + 2310*(81-56*x+6*x**2)*mj**4 +70*x*(-144+108*x-20*x**2+x**3) - 12*mj**2*(-4566+9898*x-3045*x**2+210*x**3))
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]

def sto_O81(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 5005*mj**6 + 1430*mj**7 - 5005*mj**4*(-7+x) - 1001*mj**5*(-19+2*x) + 385*mj**3*(131-36*x+2*x**2) + 385*mj**2*(112-41*x+3*x**2)+mj*(22356-12488*x+1785*x**2-70*x**3)-35*(-144+108*x-20*x**2+x**3))  * sp.Rational(1,4)
    cst2 = (-5005*mj**6 + 1430*mj**7 + 5005*mj**4*(-7+x) - 1001*mj**5*(-19+2*x) + 385*mj**3*(131-36*x+2*x**2) - 385*mj**2*(112-41*x+3*x**2)+mj*(22356-12488*x+1785*x**2-70*x**3)+35*(-144+108*x-20*x**2+x**3))  * sp.Rational(1,4)
    a = myket.op_Jp().times_cst(cst1)
    b = myket.op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O81m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 5005*mj**6 + 1430*mj**7 - 5005*mj**4*(-7+x) - 1001*mj**5*(-19+2*x) + 385*mj**3*(131-36*x+2*x**2) + 385*mj**2*(112-41*x+3*x**2)+mj*(22356-12488*x+1785*x**2-70*x**3)-35*(-144+108*x-20*x**2+x**3))  * -sp.I/4
    cst2 = (-5005*mj**6 + 1430*mj**7 + 5005*mj**4*(-7+x) - 1001*mj**5*(-19+2*x) + 385*mj**3*(131-36*x+2*x**2) - 385*mj**2*(112-41*x+3*x**2)+mj*(22356-12488*x+1785*x**2-70*x**3)+35*(-144+108*x-20*x**2+x**3))  * -sp.I/4
    a = myket.op_Jp().times_cst(cst1)
    b = myket.op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O82(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 858*mj**5 + 143*mj**6 - 143*mj**4*(-22+x) - 572*mj**3*(-12+x) + 11*mj**2*(853 - 119*x+3*x**2) + 22*mj*(333 -67*x+3*x**2) - (-2520 + 702*x-53*x**2+x**3))  * sp.Rational(1,2)
    cst2 = (-858*mj**5 + 143*mj**6 - 143*mj**4*(-22+x) + 572*mj**3*(-12+x) + 11*mj**2*(853 - 119*x+3*x**2) - 22*mj*(333 -67*x+3*x**2) - (-2520 + 702*x-53*x**2+x**3))  * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O82m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 858*mj**5 + 143*mj**6 - 143*mj**4*(-22+x) - 572*mj**3*(-12+x) + 11*mj**2*(853 - 119*x+3*x**2) + 22*mj*(333 -67*x+3*x**2) - (-2520 + 702*x-53*x**2+x**3))  *  -sp.I/2
    cst2 = (-858*mj**5 + 143*mj**6 - 143*mj**4*(-22+x) + 572*mj**3*(-12+x) + 11*mj**2*(853 - 119*x+3*x**2) - 22*mj*(333 -67*x+3*x**2) - (-2520 + 702*x-53*x**2+x**3))  *  -sp.I/2
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O83(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = (( 3+2*mj)*(840 + 234*mj**3 + 39*mj**4 -78*(-15+x)*mj -106*x +3*x**2 -13*mj**2*(-57+2*x))) * sp.Rational(1,4)
    cst2 = ((-3+2*mj)*(840 - 234*mj**3 + 39*mj**4 +78*(-15+x)*mj -106*x +3*x**2 -13*mj**2*(-57+2*x))) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O83m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = (( 3+2*mj)*(840 + 234*mj**3 + 39*mj**4 -78*(-15+x)*mj -106*x +3*x**2 -13*mj**2*(-57+2*x))) * -sp.I/4
    cst2 = ((-3+2*mj)*(840 - 234*mj**3 + 39*mj**4 +78*(-15+x)*mj -106*x +3*x**2 -13*mj**2*(-57+2*x))) * -sp.I/4
    a = myket.op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O84(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 520*mj**3 + 65*mj**4 - 13*(-139 + 2*x)*mj**2 -52*(-59+2*x)*mj + 2100 -122*x +x**2) * sp.Rational(1,2)
    cst2 = (-520*mj**3 + 65*mj**4 - 13*(-139 + 2*x)*mj**2 +52*(-59+2*x)*mj + 2100 -122*x +x**2) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O84m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 520*mj**3 + 65*mj**4 - 13*(-139 + 2*x)*mj**2 -52*(-59+2*x)*mj + 2100 -122*x +x**2) * -sp.I/2
    cst2 = (-520*mj**3 + 65*mj**4 - 13*(-139 + 2*x)*mj**2 +52*(-59+2*x)*mj + 2100 -122*x +x**2) * -sp.I/2
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O85(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = (( 5 + 2*mj)*(42 + 25*mj + 5*mj**2 -x)) * sp.Rational(1,4)
    cst2 = ((-5 + 2*mj)*(42 - 25*mj + 5*mj**2 -x)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O85m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = (( 5 + 2*mj)*(42 + 25*mj + 5*mj**2 -x)) * -sp.I/4
    cst2 = ((-5 + 2*mj)*(42 - 25*mj + 5*mj**2 -x)) * -sp.I/4
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O86(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 90*mj + 15*mj**2 -(-147+x)) * sp.Rational(1,2)
    cst2 = (-90*mj + 15*mj**2 -(-147+x)) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O86m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst1 = ( 90*mj + 15*mj**2 -(-147+x)) * (-sp.I/2)
    cst2 = (-90*mj + 15*mj**2 -(-147+x)) * (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O87(myket):
    mj = myket.Mj
    cst1 = ( 7+2*mj) * sp.Rational(1,4)
    cst2 = (-7+2*mj) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O87m(myket):
    mj = myket.Mj
    cst1 = ( 7+2*mj)* (-sp.I/4)
    cst2 = (-7+2*mj)* (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O88(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O88m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O100(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    new_coef = 46189*mj**10 - 36465*mj**8*(-22+3*x) + 3003*mj**6*(1199-450*x+30*x**2)-715*mj**4*(-6248+5481*x-966*x**2+42*x**3) - 63*x*(2880-2304*x+508*x**2-40*x**3+x**4) +33*mj**2*(32208 - 78900*x+29680*x**2-3290*x**3+105*x**4)
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]

def sto_O101(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 37791*mj**8 +8398*mj**9 -27846*mj**6*(-21+2*x) -5304*mj**7*(-41+3*x) +546*mj**5*(2603-474*x+18*x**2) +819*mj**4*(2661-620*x+30*x**2) -234*mj**2*(-8186+3640*x-427*x**2+14*x**3) -52*mj**3*(-49073+17073*x-1596*x**2+42*x**3) +63*(2880-2304*x+508*x**2-40*x**3+x**4) +6*mj*(146904-90504*x+15946*x**2-1022*x**3+21*x**4))* sp.Rational(1,4)
    cst2 = (-37791*mj**8 +8398*mj**9 +27846*mj**6*(-21+2*x) -5304*mj**7*(-41+3*x) +546*mj**5*(2603-474*x+18*x**2) -819*mj**4*(2661-620*x+30*x**2) +234*mj**2*(-8186+3640*x-427*x**2+14*x**3) -52*mj**3*(-49073+17073*x-1596*x**2+42*x**3) -63*(2880-2304*x+508*x**2-40*x**3+x**4) +6*mj*(146904-90504*x+15946*x**2-1022*x**3+21*x**4))* sp.Rational(1,4)
    a = myket.op_Jp().times_cst(cst1)
    b = myket.op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O101m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 37791*mj**8 +8398*mj**9 -27846*mj**6*(-21+2*x) -5304*mj**7*(-41+3*x) +546*mj**5*(2603-474*x+18*x**2) +819*mj**4*(2661-620*x+30*x**2) -234*mj**2*(-8186+3640*x-427*x**2+14*x**3) -52*mj**3*(-49073+17073*x-1596*x**2+42*x** 3) +63*(2880-2304*x+508*x**2-40*x**3+x**4) +6*mj*(146904-90504*x+15946*x**2-1022*x**3+21*x**4))* (-sp.I/4)
    cst2 = (-37791*mj**8 +8398*mj**9 +27846*mj**6*(-21+2*x) -5304*mj**7*(-41+3*x) +546*mj**5*(2603-474*x+18*x**2) -819*mj**4*(2661-620*x+30*x**2) +234*mj**2*(-8186+3640*x-427*x**2+14*x**3) -52*mj**3*(-49073+17073*x-1596*x**2+42*x** 3) -63*(2880-2304*x+508*x**2-40*x**3+x**4) +6*mj*(146904-90504*x+15946*x**2-1022*x**3+21*x**4))* (-sp.I/4)
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O102(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 33592*mj**7 +4199*mj**8 -3094*mj**6*(-59+2*x) -6188*mj**5*(-101+6*x) +91*mj**4*(16601-1640*x+30*x**2) +364*mj**3*(6877-960*x+30*x**2) -26*mj**2*(-106074+20230*x-1057*x**2+14*x**3) -52*mj*(-34956+8694*x-637*x**2+14*x**3) +7*(77760-24768*x+2484*x**2-92*x**3+x**4))* sp.Rational(1,2)
    cst2 = (-33592*mj**7 +4199*mj**8 -3094*mj**6*(-59+2*x) +6188*mj**5*(-101+6*x) +91*mj**4*(16601-1640*x+30*x**2) -364*mj**3*(6877-960*x+30*x**2) -26*mj**2*(-106074+20230*x-1057*x**2+14*x**3) +52*mj*(-34956+8694*x-637*x**2+14*x**3) +7*(77760-24768*x+2484*x**2-92*x**3+x**4))* sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O102m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 33592*mj**7 +4199*mj**8 -3094*mj**6*(-59+2*x) -6188*mj**5*(-101+6*x) +91*mj**4*(16601-1640*x+30*x**2) +364*mj**3*(6877-960*x+30*x**2) -26*mj**2*(-106074+20230*x-1057*x**2+14*x**3) -52*mj*(-34956+8694*x-637*x**2+14*x** 3) +7*(77760-24768*x+2484*x**2-92*x**3+x**4))* (-sp.I/2)
    cst2 = (-33592*mj**7 +4199*mj**8 -3094*mj**6*(-59+2*x) +6188*mj**5*(-101+6*x) +91*mj**4*(16601-1640*x+30*x**2) -364*mj**3*(6877-960*x+30*x**2) -26*mj**2*(-106074+20230*x-1057*x**2+14*x**3) +52*mj*(-34956+8694*x-637*x**2+14*x** 3) +7*(77760-24768*x+2484*x**2-92*x**3+x**4))* (-sp.I/2)
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O103(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 6783*mj**6 +646*mj**7 -1785*mj**4*(-79+3*x) -119*mj**5*(-329+6*x) +63*mj**2*(8054-735*x+15*x**2) +7*mj**3*(47677-3000*x+30*x**2) +mj*(453024-56154*x+1897*x**2-14*x**3) -21*(-8640+1392*x-68*x**2+x**3))* sp.Rational(1,4)
    cst2 = (-6783*mj**6 +646*mj**7 +1785*mj**4*(-79+3*x) -119*mj**5*(-329+6*x) -63*mj**2*(8054-735*x+15*x**2) +7*mj**3*(47677-3000*x+30*x**2) +mj*(453024-56154*x+1897*x**2-14*x**3) +21*(-8640+1392*x-68*x**2+x**3))* sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O103m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 6783*mj**6 +646*mj**7 -1785*mj**4*(-79+3*x) -119*mj**5*(-329+6*x) +63*mj**2*(8054-735*x+15*x**2) +7*mj**3*(47677-3000*x+30*x**2) +mj*(453024-56154*x+1897*x**2-14*x**3) -21*(-8640+1392*x-68*x**2+x**3))* (-sp.I/4)
    cst2 = (-6783*mj**6 +646*mj**7 +1785*mj**4*(-79+3*x) -119*mj**5*(-329+6*x) -63*mj**2*(8054-735*x+15*x**2) +7*mj**3*(47677-3000*x+30*x**2) +mj*(453024-56154*x+1897*x**2-14*x**3) +21*(-8640+1392*x-68*x**2+x**3))* (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O104(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 3876*mj**5 + 323*mj**6 -2040*mj**3*(-39+x) -85*mj**4*(-269+3*x) +12*mj*(16792-1075*x+15*x**2) +mj**2*(168152-7305*x+45*x**2) -(-105840+9252*x-218*x**2+x**3)) * sp.Rational(1,2)
    cst2 = (-3876*mj**5 + 323*mj**6 +2040*mj**3*(-39+x) -85*mj**4*(-269+3*x) -12*mj*(16792-1075*x+15*x**2) +mj**2*(168152-7305*x+45*x**2) -(-105840+9252*x-218*x**2+x**3)) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O104m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 3876*mj**5 + 323*mj**6 -2040*mj**3*(-39+x) -85*mj**4*(-269+3*x) +12*mj*(16792-1075*x+15*x**2) +mj**2*(168152-7305*x+45*x**2) -(-105840+9252*x-218*x**2+x**3))* (-sp.I/2)
    cst2 = (-3876*mj**5 + 323*mj**6 +2040*mj**3*(-39+x) -85*mj**4*(-269+3*x) -12*mj*(16792-1075*x+15*x**2) +mj**2*(168152-7305*x+45*x**2) -(-105840+9252*x-218*x**2+x**3))* (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O105(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 8075*mj**4 +646*mj**5 -340*mj**3*(-134+x) -425*mj**2*(-329+6*x) +15*(10584-500*x+5*x**2) +2*mj*(114627-3625*x+15*x**2)) * sp.Rational(1,4)
    cst2 = (-8075*mj**4 +646*mj**5 -340*mj**3*(-134+x) +425*mj**2*(-329+6*x) -15*(10584-500*x+5*x**2) +2*mj*(114627-3625*x+15*x**2)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O105m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 8075*mj**4 +646*mj**5 -340*mj**3*(-134+x) -425*mj**2*(-329+6*x) +15*(10584-500*x+5*x**2) +2*mj*(114627-3625*x+15*x**2))* (-sp.I/4)
    cst2 = (-8075*mj**4 +646*mj**5 -340*mj**3*(-134+x) +425*mj**2*(-329+6*x) -15*(10584-500*x+5*x**2) +2*mj*(114627-3625*x+15*x**2))* (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]


def sto_O106(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 3876*mj**3 +323*mj**4 -17*mj**2*(-1127+6*x) -102*mj*(-443+6*x) +3*(14112-338*x+x**2)) * sp.Rational(1,2)
    cst2 = (-3876*mj**3 +323*mj**4 -17*mj**2*(-1127+6*x) +102*mj*(-443+6*x) +3*(14112-338*x+x**2)) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O106m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 3876*mj**3 +323*mj**4 -17*mj**2*(-1127+6*x) -102*mj*(-443+6*x) +3*(14112-338*x+x**2))* (-sp.I/2)
    cst2 = (-3876*mj**3 +323*mj**4 -17*mj**2*(-1127+6*x) +102*mj*(-443+6*x) +3*(14112-338*x+x**2))* (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O107(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = (( 7+2*mj)*(288 + 133*mj + 19*mj**2 - 3*x)) * sp.Rational(1,4)
    cst2 = ((-7+2*mj)*(288 - 133*mj + 19*mj**2 - 3*x)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O107m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = (( 7+2*mj)*(288 + 133*mj + 19*mj**2 - 3*x))* (-sp.I/4)
    cst2 = ((-7+2*mj)*(288 - 133*mj + 19*mj**2 - 3*x))* (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O108(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 324+152*mj+19*mj**2-x) * sp.Rational(1,2)
    cst2 = ( 324-152*mj+19*mj**2-x) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O108m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 324+152*mj+19*mj**2-x)* (-sp.I/2)
    cst2 = ( 324-152*mj+19*mj**2-x)* (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O109(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 9+2*mj) * sp.Rational(1,4)
    cst2 = (-9+2*mj) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O109m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 9+2*mj)* (-sp.I/4)
    cst2 = (-9+2*mj)* (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1010(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1010m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1212(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1212m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1211(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 11+2*mj) * sp.Rational(1,4)
    cst2 = (-11+2*mj) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1211m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 11+2*mj)* (-sp.I/4)
    cst2 = (-11+2*mj)* (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1210(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = (605+230*mj+23*mj**2-x) * sp.Rational(1,2)
    cst2 = (605-230*mj+23*mj**2-x) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O1210m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = (605+230*mj+23*mj**2-x) * (-sp.I/2)
    cst2 = (605-230*mj+23*mj**2-x) * (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O129(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = (( 9+2*mj)*(550+207*mj+23*mj**2-3*x)) * sp.Rational(1,4)
    cst2 = ((-9+2*mj)*(550-207*mj+23*mj**2-3*x)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O129m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = (( 9+2*mj)*(550+207*mj+23*mj**2-3*x)) * (-sp.I/4)
    cst2 = ((-9+2*mj)*(550-207*mj+23*mj**2-3*x)) * (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O128(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 2576*mj**3 +161*mj**4 -7*mj**2*(-2365+6*x) -56*mj*(-893+6*x) +(59400-722*x+x**2)) * sp.Rational(1,2)
    cst2 = (-2576*mj**3 +161*mj**4 -7*mj**2*(-2365+6*x) +56*mj*(-893+6*x) +(59400-722*x+x**2)) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O128m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 2576*mj**3 +161*mj**4 -7*mj**2*(-2365+6*x) -56*mj*(-893+6*x) +(59400-722*x+x**2)) * (-sp.I/2)
    cst2 = (-2576*mj**3 +161*mj**4 -7*mj**2*(-2365+6*x) +56*mj*(-893+6*x) +(59400-722*x+x**2)) * (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O127(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 5635*mj**4 +322*mj**5 -140*mj**3*(-306+x) -245*mj**2*(-709+6*x) +35*(9504-218*x+x**2) +2*mj*(185749-2805*x+5*x**2)) * sp.Rational(1,4)
    cst2 = (-5635*mj**4 +322*mj**5 -140*mj**3*(-306+x) +245*mj**2*(-709+6*x) -35*(9504-218*x+x**2) +2*mj*(185749-2805*x+5*x**2)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O127m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 =( 5635*mj**4 +322*mj**5 -140*mj**3*(-306+x) -245*mj**2*(-709+6*x) +35*(9504-218*x+x**2) +2*mj*(185749-2805*x+5*x**2)) * (-sp.I/4)
    cst2 =(-5635*mj**4 +322*mj**5 -140*mj**3*(-306+x) +245*mj**2*(-709+6*x) -35*(9504-218*x+x**2) +2*mj*(185749-2805*x+5*x**2)) * (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O126(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 55062*mj**5 +3059*mj**6 -665*mj**4*(-688+3*x) -7980*mj**3*(-274+3*x) +19*mj**2*(328739-6315*x+15*x**2) +114*mj*(87827-2535*x+15*x**2) -5*(-1397088+55578*x-575*x**2+x**3)) * sp.Rational(1,2)
    cst2 = (-55062*mj**5 +3059*mj**6 -665*mj**4*(-688+3*x) +7980*mj**3*(-274+3*x) +19*mj**2*(328739-6315*x+15*x**2) -114*mj*(87827-2535*x+15*x**2) -5*(-1397088+55578*x-575*x**2+x**3)) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O126m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 55062*mj**5 +3059*mj**6 -665*mj**4*(-688+3*x) -7980*mj**3*(-274+3*x) +19*mj**2*(328739-6315*x+15*x**2) +114*mj*(87827-2535*x+15*x**2) -5*(-1397088+55578*x-575*x**2+x**3)) * (-sp.I/2)
    cst2 = (-55062*mj**5 +3059*mj**6 -665*mj**4*(-688+3*x) +7980*mj**3*(-274+3*x) +19*mj**2*(328739-6315*x+15*x**2) -114*mj*(87827-2535*x+15*x**2) -5*(-1397088+55578*x-575*x**2+x**3)) * (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O125(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 15295*mj**6 +874*mj**7 -3325*mj**4*(-205+3*x) -133*mj**5*(-985+6*x) +19*mj**3*(120079-3020*x+10*x**2) +95*mj**2*(51056-1905*x+15*x**2) +mj*(6010260-306652*x+4135*x**2-10*x**3) -5*(-665280+44076*x-880*x**2+5*x**3)) * sp.Rational(1,4)
    cst2 = (-15295*mj**6 +874*mj**7 +3325*mj**4*(-205+3*x) -133*mj**5*(-985+6*x) +19*mj**3*(120079-3020*x+10*x**2) -95*mj**2*(51056-1905*x+15*x**2) +mj*(6010260-306652*x+4135*x**2-10*x**3) +5*(-665280+44076*x-880*x**2+5*x**3)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O125m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 15295*mj**6 +874*mj**7 -3325*mj**4*(-205+3*x) -133*mj**5*(-985+6*x) +19*mj**3*(120079-3020*x+10*x**2) +95*mj**2*(51056-1905*x+15*x**2) +mj*(6010260-306652*x+4135*x**2-10*x**3) -5*(-665280+44076*x-880*x**2+5*x**3)) * (-sp.I/4)
    cst2 = (-15295*mj**6 +874*mj**7 +3325*mj**4*(-205+3*x) -133*mj**5*(-985+6*x) +19*mj**3*(120079-3020*x+10*x**2) -95*mj**2*(51056-1905*x+15*x**2) +mj*(6010260-306652*x+4135*x**2-10*x**3) +5*(-665280+44076*x-880*x**2+5*x**3)) * (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O124(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 118864*mj**7 +7429*mj**8 -4522*mj**6*(-221+2*x) -18088*mj**5*(-295+6*x) +323*mj**4*(59767-2040*x+10*x**2) +2584*mj**3*(18439-920*x+10*x**2) -34*mj**2*(-2279002+154234*x-2805*x**2+10*x**3) -136*mj*(-553650+48442*x-1285*x**2+10*x**3) +5*(6652800-728400*x+26188*x**2-340*x**3+x**4)) * sp.Rational(1,2)
    cst2 = (-118864*mj**7 +7429*mj**8 -4522*mj**6*(-221+2*x) +18088*mj**5*(-295+6*x) +323*mj**4*(59767-2040*x+10*x**2) -2584*mj**3*(18439-920*x+10*x**2) -34*mj**2*(-2279002+154234*x-2805*x**2+10*x**3) +136*mj*(-553650+48442*x-1285*x**2+10*x**3) +5*(6652800-728400*x+26188*x**2-340*x**3+x**4)) * sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O124m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 118864*mj**7 +7429*mj**8 -4522*mj**6*(-221+2*x) -18088*mj**5*(-295+6*x) +323*mj**4*(59767-2040*x+10*x**2) +2584*mj**3*(18439-920*x+10*x**2) -34*mj**2*(-2279002+154234*x-2805*x**2+10*x**3) -136*mj*(-553650+48442*x-1285*x**2+10*x**3) +5*(6652800-728400*x+26188*x**2-340*x**3+x**4)) * (-sp.I/2)
    cst2 = (-118864*mj**7 +7429*mj**8 -4522*mj**6*(-221+2*x) +18088*mj**5*(-295+6*x) +323*mj**4*(59767-2040*x+10*x**2) -2584*mj**3*(18439-920*x+10*x**2) -34*mj**2*(-2279002+154234*x-2805*x**2+10*x**3) +136*mj*(-553650+48442*x-1285*x**2+10*x**3) +5*(6652800-728400*x+26188*x**2-340*x**3+x**4)) * (-sp.I/2)
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O123(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 200583*mj**8 +14858*mj**9 -7752*mj**7*(-205+3*x) -40698*mj**6*(-203+6*x) +1938*mj**5*(15689-762*x+6*x**2) +2907*mj**4*(27541-1920*x+30*x**2) -68*mj**3*(-2190575+205803*x-5280*x**2+30*x**3) -306*mj**2*(-610706+73962*x-2715*x**2+30*x**3) +45*(1108800-206400*x+13148*x**2-340*x**3+3*x**4) +6*mj*(23730600-3600564*x+178042*x**2-3230*x**3+15*x**4))  * sp.Rational(1,4)
    cst2 = (-200583*mj**8 +14858*mj**9 -7752*mj**7*(-205+3*x) +40698*mj**6*(-203+6*x) +1938*mj**5*(15689-762*x+6*x**2) -2907*mj**4*(27541-1920*x+30*x**2) -68*mj**3*(-2190575+205803*x-5280*x**2+30*x**3) +306*mj**2*(-610706+73962*x- 2715*x**2+30*x**3) -45*(1108800-206400*x+13148*x**2-340*x**3+3*x**4) +6*mj*(23730600-3600564*x+178042*x**2-3230*x**3+15*x**4)) * sp.Rational(1,4)
    a = myket.op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O123m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 200583*mj**8 +14858*mj**9 -7752*mj**7*(-205+3*x) -40698*mj**6*(-203+6*x) +1938*mj**5*(15689-762*x+6*x**2) +2907*mj**4*(27541-1920*x+30*x**2) -68*mj**3*(-2190575+205803*x-5280*x**2+30*x**3) -306*mj**2*(-610706+73962*x- 2715*x**2+30*x**3) +45*(1108800-206400*x+13148*x**2-340*x**3+3*x**4) +6*mj*(23730600-3600564*x+178042*x**2-3230*x**3+15*x**4))  * (-sp.I/4)
    cst2 = (-200583*mj**8 +14858*mj**9 -7752*mj**7*(-205+3*x) +40698*mj**6*(-203+6*x) +1938*mj**5*(15689-762*x+6*x**2) -2907*mj**4*(27541-1920*x+30*x**2) -68*mj**3*(-2190575+205803*x-5280*x**2+30*x**3) +306*mj**2*(-610706+73962*x- 2715*x**2+30*x**3) -45*(1108800-206400*x+13148*x**2-340*x**3+3*x**4) +6*mj*(23730600-3600564*x+178042*x**2-3230*x**3+15*x**4)) * (-sp.I/4)
    a = myket.op_Jp().op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O122(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 74290*mj**9 +7429*mj**10 -4845*mj**8*(-113+3*x) -38760*mj**7*(-67+3*x) +9690*mj**5*(2501-258*x+6*x**2) +969*mj**6*(9563-710*x+10*x**2) -85*mj**4*(-554219+77439*x-3000*x**2+30*x**3) -340*mj**3*(-193637+34803*x-1860*x**2+30*x**3) +30*mj*(1228584-346356*x+32182*x**2-1190*x**3+15*x**4) +3*mj**2*(20983308-4771720*x+345870*x**2-9350*x**3+75*x**4)-3*(-3326400+1152000*x-135444*x**2+6808*x**3-145*x**4+x**5)) * sp.Rational(1,2)
    cst2 = (-74290*mj**9 +7429*mj**10 -4845*mj**8*(-113+3*x) +38760*mj**7*(-67+3*x) -9690*mj**5*(2501-258*x+6*x**2) +969*mj**6*(9563-710*x+10*x**2) -85*mj**4*(-554219+77439*x-3000*x**2+30*x**3) +340*mj**3*(-193637+34803*x-1860*x**2+30*x**3) -30*mj*(1228584-346356*x+32182*x**2-1190*x**3+15*x**4) +3*mj**2*(20983308-4771720*x+345870*x**2-9350*x**3+75*x**4)-3*(-3326400+1152000*x-135444*x**2+6808*x**3-145*x**4+x**5))* sp.Rational(1,2)
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O122m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 74290*mj**9 +7429*mj**10 -4845*mj**8*(-113+3*x) -38760*mj**7*(-67+3*x) +9690*mj**5*(2501-258*x+6*x**2) +969*mj**6*(9563-710*x+10*x**2) -85*mj**4*(-554219+77439*x-3000*x**2+30*x**3) -340*mj**3*(-193637+34803*x-1860*x**2+30*x**3) +30*mj*(1228584-346356*x+32182*x**2-1190*x**3+15*x**4) +3*mj**2*(20983308-4771720*x+345870*x**2-9350*x**3+75*x**4)-3*(-3326400+1152000*x-135444*x**2+6808*x**3-145*x**4+x**5))* (-sp.I/2)
    cst2 = (-74290*mj**9 +7429*mj**10 -4845*mj**8*(-113+3*x) +38760*mj**7*(-67+3*x) -9690*mj**5*(2501-258*x+6*x**2) +969*mj**6*(9563-710*x+10*x**2) -85*mj**4*(-554219+77439*x-3000*x**2+30*x**3) +340*mj**3*(-193637+34803*x-1860*x**2+30*x**3) -30*mj*(1228584-346356*x+32182*x**2-1190*x**3+15*x**4) +3*mj**2*(20983308-4771720*x+345870*x**2-9350*x**3+75*x**4)-3*(-3326400+1152000*x-135444*x**2+6808*x**3-145*x**4+x**5)) * (-sp.I/2)
    a = myket.op_Jp().op_Jp().times_cst(cst1)
    b = myket.op_Jm().op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O121(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 572033*mj**10 +104006*mj**11 -124355*mj**9*(-37+2*x) -373065*mj**8*(-44+3*x) +42636*mj**7*(1383-180*x+5*x**2) +74613*mj**6*(1793-290*x+10*x**2) -6545*mj**5*(-39755+9498*x-630*x**2+12*x**3) -6545*mj**4*(-53702+15879*x-1290*x**2+30*x**3) +154*mj**3*(2341414-937880*x+113055*x**2-5100*x**3+75*x**4) +231*mj**2*(1065768-531580*x+78120*x**2-4250*x**3+75*x**4) -231*(-86400+72000*x-17544*x**2+1708*x**3-70*x**4+x**5) -3*mj*(-34637280+22740960*x-4531076*x**2+367752*x**3-12705*x**4+154*x**5)) * sp.Rational(1,4)

    cst2 = (-572033*mj**10 +104006*mj**11 -124355*mj**9*(-37+2*x) +373065*mj**8*(-44+3*x) +42636*mj**7*(1383-180*x+5*x**2) -74613*mj**6*(1793-290*x+10*x**2) -6545*mj**5*(-39755+9498*x-630*x**2+12*x**3) +6545*mj**4*(-53702+15879*x- 1290*x**2+30*x**3) +154*mj**3*(2341414-937880*x+113055*x**2-5100*x**3+75*x**4) -231*mj**2*(1065768-531580*x+78120*x**2-4250*x**3+75*x**4) +231*(-86400+72000*x-17544*x**2+1708*x**3-70*x**4+x**5) -3*mj*(-34637280+22740960*x-4531076*x**2+367752*x**3-12705*x**4+154*x**5)) * sp.Rational(1,4)
    a = myket.op_Jp().times_cst(cst1)
    b = myket.op_Jm().times_cst(cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O121m(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1) 
    cst1 = ( 572033*mj**10 +104006*mj**11 -124355*mj**9*(-37+2*x) -373065*mj**8*(-44+3*x) +42636*mj**7*(1383-180*x+5*x**2) +74613*mj**6*(1793-290*x+10*x**2) -6545*mj**5*(-39755+9498*x-630*x**2+12*x**3) -6545*mj**4*(-53702+15879*x- 1290*x**2+30*x**3) +154*mj**3*(2341414-937880*x+113055*x**2-5100*x**3+75*x**4) +231*mj**2*(1065768-531580*x+78120*x**2-4250*x**3+75*x**4) -231*(-86400+72000*x-17544*x**2+1708*x**3-70*x**4+x**5) -3*mj*(-34637280+22740960*x-4531076*x**2+367752*x**3-12705*x**4+154*x**5)) * (-sp.I/4)
    cst2 = (-572033*mj**10 +104006*mj**11 -124355*mj**9*(-37+2*x) +373065*mj**8*(-44+3*x) +42636*mj**7*(1383-180*x+5*x**2) -74613*mj**6*(1793-290*x+10*x**2) -6545*mj**5*(-39755+9498*x-630*x**2+12*x**3) +6545*mj**4*(-53702+15879*x- 1290*x**2+30*x**3) +154*mj**3*(2341414-937880*x+113055*x**2-5100*x**3+75*x**4) -231*mj**2*(1065768-531580*x+78120*x**2-4250*x**3+75*x**4) +231*(-86400+72000*x-17544*x**2+1708*x**3-70*x**4+x**5) -3*mj*(-34637280+22740960*x-4531076*x**2+367752*x**3-12705*x**4+154*x**5)) * (-sp.I/4)
    a = myket.op_Jp().times_cst(cst1)
    b = myket.op_Jm().times_cst(-1*cst2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O120(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    new_coef = 676039*mj**12 -323323*mj**10*(-65+6*x) +138567*mj**8*(1391-330*x+15*x**2)-17017*mj**6*(-35945+17622*x-2010*x**2+60*x**3) +1001*mj**4*(606164-618090*x+139245*x**2-10200*x**3+225*x**4) +231*x*(-86400+72000*x-17544*x**2+1708*x**3-70*x**4+x**5) -39*mj**2*(-3176160+8488392*x-3601048*x**2+501116*x**3-26565*x**4+462*x**5)
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]


"""
Sympy symbols used in the construction of the Stevens CF
and Zeeman matrices across the NewMag project
"""

# SO constant
Z = sp.symbols("Z",real=True)

# Magnetic couplings
J0, J1, J2, J3 = sp.symbols("J_0 J_1 J_2 J_3",real=True)

# Stevens CF parameters
B20, B21, B22 = sp.symbols("B_2^0 B_2^1 B_2^2",real=True)
B21m, B22m = sp.symbols("B_2^-1 B_2^-2",real=True)

B40, B41, B42, B43, B44 = sp.symbols("B_4^0 B_4^1 B_4^2 B_4^3 B_4^4", real=True)
B41m, B42m, B43m, B44m = sp.symbols("B_4^-1 B_4^-2, B_4^-3 B_4^-4", real=True)

B50, B51, B52, B53, B54, B55 = sp.symbols("B_5^0 B_5^1 B_5^2 B_5^3 B_5^4 B_5^5", real=True)
B51m, B52m, B53m, B54m, B55m = sp.symbols("B_5^-1 B_5^-2, B_5^-3 B_5^-4 B_5^-5", real=True)

B60, B61, B62, B63, B64, B65, B66 = sp.symbols("B_6^0 B_6^1 B_6^2 B_6^3 B_6^4 B_6^5 B_6^6", real=True)
B61m, B62m, B63m, B64m, B65m, B66m = sp.symbols("B_6^-1 B_6^-2, B_6^-3 B_6^-4 B_6^-5 B_6^-6", real=True)

B80, B81, B82, B83, B84, B85, B86, B87, B88 = sp.symbols("B_8^0 B_8^1 B_8^2 B_8^3 B_8^4 B_8^5 B_8^6 B_8^7 B_8^8", real=True)
B81m, B82m, B83m, B84m, B85m, B86m, B87m, B88m = sp.symbols("B_8^-1 B_8^-2, B_8^-3 B_8^-4 B_8^-5 B_8^-6 B_8^-7 B_8^-8", real=True)

B100, B101, B102, B103, B104, B105, B106, B107, B108, B109, B1010 = sp.symbols("B_10^0 B_10^1 B_10^2 B_10^3 B_10^4 B_10^5 B_10^6 B_10^7 B_10^8 B_10^9 B_10^10", real=True)
B101m, B102m, B103m, B104m, B105m, B106m, B107m, B108m, B109m, B1010m = sp.symbols("B_10^-1 B_10^-2, B_10^-3 B_10^-4 B_10^-5 B_10^-6 B_10^-7 B_10^-8 B_10^-9 B_10^-10", real=True)

B120, B121, B122, B123, B124, B125, B126, B127, B128, B129, B1210, B1211, B1212 = sp.symbols("B_12^0 B_12^1 B_12^2 B_12^3 B_12^4 B_12^5 B_12^6 B_12^7 B_12^8 B_12^9 B_12^10 B_12^11 B_12^12", real=True)
B121m, B122m, B123m, B124m, B125m, B126m, B127m, B128m, B129m, B1210m, B1211m, B1212m = sp.symbols("B_12^-1 B_12^-2, B_12^-3 B_12^-4 B_12^-5 B_12^-6 B_12^-7 B_12^-8 B_12^-9 B_12^-10 B_12^-11 B_12^-12", real=True)


# g-factors
gxx, gyy, gzz, gxy, gyx, gxz, gzx, gyz, gzy = sp.symbols("g_xx g_yy g_zz g_xy g_yx g_xz g_zx g_yz g_zy", real=True)

# Magnetic field in x, y, z directions
Bx, By, Bz = sp.symbols("B_x B_y B_z", real=True)
 
# Lande g-factor
gJ = sp.symbols("g_J", real=True)
 
# Bohr magneton
BM = sp.symbols("\\mu_B", real=True)


# Gather tuples for the rank and projection of the operator
Okq = {
    (2, 0): sto_O20, (2, 1): sto_O21, (2, 2): sto_O22, (2,-1): sto_O21m, (2,-2): sto_O22m,

    (4, 0): sto_O40, (4, 1): sto_O41, (4, 2): sto_O42, (4, 3): sto_O43, (4, 4): sto_O44, 
    (4,-1): sto_O41m, (4,-2): sto_O42m, (4,-3): sto_O43m, (4,-4): sto_O44m,

    (6, 0): sto_O60, (6, 1): sto_O61, (6, 2): sto_O62, (6, 3): sto_O63, (6, 4): sto_O64, (6, 5): sto_O65, (6, 6): sto_O66, 
    (6,-1): sto_O61m, (6,-2): sto_O62m, (6,-3): sto_O63m, (6,-4): sto_O64m, (6,-5): sto_O65m, (6,-6): sto_O66m,
    
    (8, 0): sto_O80, 
    (8, 1): sto_O81, (8,  2): sto_O82, (8,  3): sto_O83, (8,  4): sto_O84, (8,  5): sto_O85, (8,  6): sto_O86, (8,  7): sto_O87, (8,  8): sto_O88, 
    (8, -1): sto_O81m, (8, -2): sto_O82m,(8, -3): sto_O83m, (8, -4): sto_O84m, (8, -5): sto_O85m,(8, -6): sto_O86m, (8,-7): sto_O87m, (8,-8): sto_O88m,

    (10, 0): sto_O100,
    (10, 1): sto_O101, (10,  2): sto_O102, (10,  3): sto_O103, (10,  4): sto_O104, (10,  5): sto_O105, (10,  6): sto_O106, (10,  7): sto_O107, (10,  8): sto_O108,(10,  9): sto_O109,(10, 10): sto_O1010,
    (10,-1): sto_O101m, (10, -2): sto_O102m, (10,  -3): sto_O103m, (10,  -4): sto_O104m, (10,  -5): sto_O105m, (10,  -6): sto_O106m, (10, -7): sto_O107m, (10,  -8): sto_O108m,(10,  -9): sto_O109m,(10, -10): sto_O1010m,

    (12, 0): sto_O120,
    (12, 1): sto_O121, (12,  2): sto_O122, (12,  3): sto_O123, (12,  4): sto_O124, (12,  5): sto_O125, (12,  6): sto_O126, (12,  7): sto_O127, (12,  8): sto_O128,(12,  9): sto_O129,(12, 10): sto_O1210, (12, 11): sto_O1211, (12, 12): sto_O1212,
    (12,-1): sto_O121m, (12, -2): sto_O122m, (12,  -3): sto_O123m, (12,  -4): sto_O124m, (12,  -5): sto_O125m, (12,  -6): sto_O126m, (12, -7): sto_O127m, (12,  -8): sto_O128m,(12,  -9): sto_O129m,(12, -10): sto_O1210m, (12, -11): sto_O1211m, (12, -12): sto_O1212m

}

# Gather the CFPs in a dictionary
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
    (12,  -1): B121m, (12,  -2): B122m, (12,  -3): B123m, (12,  -4): B124m, (12,  -5): B125m, (12,  -6): B126m,  (12, -7): B127m,  (12, -8): B128m, (12, -9): B129m,(12, -10): B1210m,(12, -11): B1211m,(12, -12): B1212m,

}

# Make a dictionary relating the STOs with the proportionality constants used in the ITOs of SINGLE_ANISO
ito_cst = {
    (2, 0): 1, (2, 1): sp.sqrt(6), (2, 2): sp.sqrt(3)/sp.sqrt(2), (2,-1): sp.sqrt(6), (2,-2): sp.sqrt(3)/sp.sqrt(2),
    
    (4, 0): 1, (4, 1): 2*sp.sqrt(5), (4, 2): sp.sqrt(10), (4, 3): 2*sp.sqrt(35), (4, 4): sp.sqrt(35)/sp.sqrt(2), 
    (4,-1): 2*sp.sqrt(5), (4,-2): sp.sqrt(10), (4,-3): 2*sp.sqrt(35), (4,-4): sp.sqrt(35)/sp.sqrt(2),
    
    (6, 0): 1, (6, 1): sp.sqrt(42), (6, 2): sp.sqrt(105)/2, (6, 3): sp.sqrt(105), (6, 4): 3*sp.sqrt(7)/sp.sqrt(2), (6, 5): 3*sp.sqrt(77), (6, 6): sp.sqrt(231)/2, 
    (6,-1): sp.sqrt(42), (6,-2): sp.sqrt(105)/2, (6,-3): sp.sqrt(105), (6,-4): 3*sp.sqrt(7)/sp.sqrt(2), (6,-5): 3*sp.sqrt(77), (6,-6): sp.sqrt(231)/2
}


