from function.helper_functions import *
from itertools import permutations
from copy import deepcopy  
import numpy as np 
import pandas as pd        
import h5py

def casscf_wfn_orca(orca_output, mag_center, pattern_wfn, ml_order, remove_conf, ion_wft, user_sr_root):
    '''
    reading the sr wavefunction form orca output
    
    INPUT
    filename : main orca output 
    mag_center : class of functions possessing all the properties characteristic of the metallic center
    pattern_wfn : pattern of the CSF in the main molcas output
    ml_order : list of the ml make up the CSF of the metal center
    remove_conf : list of CSF not on the metal center
    ion_wft : ion orca output  
    user_root : sr root selected by the user

    OUTPUT
    sr_wfn : dict[ml][coeff], sr wavefunction
    root_sr_remove : sr root not selected in the model
    '''
    sr_wfn = {}  # {state_index: {config: coeff}}
    weight_wfn = {}
    marker_CI=False
    root_cas=number_of_CI_CAS_orca(orca_output,mag_center)
    len_CFS=len(remove_conf)+2*mag_center.l+1
    root=-1
    with open(orca_output, 'r') as f:
        for line in f:
            if not marker_CI:
                if "Spin-Determinant CI Printing" in line:
                    marker_CI=True
                continue
            if "ROOT" and ":  E=" in line:
                root+=1#int(line.split()[1].split(':')[0])
                sr_wfn[root]={}
                weight_wfn[root]=0
                for base in mag_center.basis_ne_ml_uncoupled :
                    sr_wfn[root][base]=[] 
                    if len(set(base)) != len(base):
                        sr_wfn[root][base] = np.float64(0)
            if "CAS-SCF STATES FOR BLOCK" in line:
                marker_CI=False
            
            if re.search(pattern_wfn,line):   
                conf=line.split()[0].split('[')[1].split(']')[0]
                if len(list(conf)) != len_CFS:
                    if len(list(conf)) != len_CFS:
                        print('ERROR: the lenght of the CFS are different of your date input')
                        print('use the cas_mo input')
                        print('There are at most 4 categories of MOs in the active space and the sum of (d+c+o+v) must match the number of active orbitals:')
                        print('d = doubly occupied -> Orbitals with a strict occupancy of 2 ("inactive")')
                        print(f'c = centered orbitals -> Orbitals at the magnetic center used to build the model space (metal-centered p, d or f orbitals), the total number of c must be equal to {2*mag.l+1}')
                        print('o = other occupied -> Orbitals belonging to the active space but not to the model one (other metal-centered orbitals or ligand-centered orbitals)')
                        print('v = vacant orbitals -> Orbitals with a strict occupancy of 0 ("virtual")')

                        general_error_inp_orca(mag_center)

                norm = set()
                for i in remove_conf:
                    if -len(conf) <= i < len(conf):
                        norm.add(i if i >= 0 else len(conf) + i)
                conf=[elem for idx, elem in enumerate(conf) if idx not in norm]
                conf=''.join(conf)
                ml_idx=[m.span()[0] for m in re.finditer(re.compile('u'), conf)]
                weight_wfn[root]+=np.float64(line.split()[1])**2
                
                ml_value=[]
                for i in range(len(ml_idx)):
                    ml_value.append((ml_order[ml_idx[i]]))
                
                dico_ml_sign = {}
                index_map = {val: i for i, val in enumerate(ml_value)}
                for perm in permutations(ml_value):
                    inv = 0
                    for i in range(len(perm)):
                        for j in range(i + 1, len(perm)):
                            if index_map[perm[i]] > index_map[perm[j]]:
                                inv += 1
                    sign = 1 if inv % 2 == 0 else -1
                    dico_ml_sign[perm] = sign
                for i,conf in enumerate(dico_ml_sign):
                    sr_wfn[root][conf] = np.float64(line.split()[1]) * dico_ml_sign[conf]

    if user_sr_root == None:
        root_select_by_OA=select_sr_root_by_OA(weight_wfn,mag_center,printing=True)
    else : 
        print('Scalar roots selected by user input')
        wfn_print={}
        wfn_print['User roots']={}
        for root in range(len(sr_wfn)):                                        
            if root in user_sr_root:                                                
                wfn_print['User roots'][f'Root {root}']='---> Selected'
            else:
                wfn_print['User roots'][f'Root {root}']=' '
        df=pd.DataFrame(wfn_print)                                                                   
        print(df.to_string(),end='\n\n')                                                   
        root_select_by_OA = user_sr_root

    all_sr_root=list(range(0,root_cas+1))
    root_sr_remove_by_OA=sorted(list(set(all_sr_root)-set(root_select_by_OA)))
    root_sr_remove = root_sr_remove_by_OA
    
    for root_idx in range(len(sr_wfn)):
        for base in mag_center.basis_ne_ml_uncoupled:
            if len(set(base)) != len(base):
                sr_wfn[root_idx][base] = int(0)
    for root in root_sr_remove:
            del sr_wfn[root]
    for key_sr, new_key in zip(root_select_by_OA,range(len(mag_center.basis_lmls))):
        sr_wfn[new_key]=sr_wfn.pop(key_sr)

    ne_ml=mag_center.basis_ne_ml_uncoupled
    if mag_center.l==3:                
        for root in range(len(sr_wfn)):
            for base in range(len(sr_wfn[root])):
                if (3 in ne_ml[base]) or (-3 in ne_ml[base]):
                    if (3 in ne_ml[base]) and (-3 in ne_ml[base]):
                        continue      
                    sr_wfn[root][ne_ml[base]] *= -1              
#    elif mag_center.l==2:                
#        for root in range(len(sr_wfn)):
#            for base in range(len(sr_wfn[root])):
#                if (2 in ne_ml[base]) or (-2 in ne_ml[base]):
#                    if (2 in ne_ml[base]) and (-2 in ne_ml[base]):
#                        continue      
#                    sr_wfn[root][ne_ml[base]] *= -1              
#

    if user_sr_root==None:
        if len(sr_wfn) < len(mag_center.basis_lmls):
            print('There are some roots missing from your output.')
            print(f'You indicated that your magnetic center has a dimension of 2L+1 = {len(mag_center.basis_lmls)}')
            print(f'Therefore, your output should contain at least {len(mag_center.basis_lmls)} roots')
            print('EXIT Program') 
            exit()
        elif len(sr_wfn) > len(mag_center.basis_lmls):
            print(f"The wavefunction has more than 2L+1 = {len(mag_center.basis_lmls)} roots to represente the fundamental term.")
            if ion_wft==None:
                print(f"We apply the pseudo-L approximation (by selecting the first 2L+1={len(mag_center.basis_lmls)} roots) \n\n")
                for root_idx in range(len(sr_wfn)):
                    if root_idx >= len(mag_center.basis_lmls):
                        del sr_wfn[root_idx] 
                        root_sr_remove.append(root_idx)
    
    if ion_wft!=None:
        if len(sr_wfn) != len(ion_wft):
            print('There are too few or too many roots in the output of the ion.')
            print(f'Your molecular output has {len(sr_wfn)} roots; therefore, your ionic output must have the same number.')
            print(f"We apply the pseudo-L approximation (by selecting the first 2L+1={len(mag_center.basis_lmls)} roots). \n\n")
            for root_idx in range(len(sr_wfn)):            
                if root_idx >= len(mag_center.basis_lmls): 
                    del sr_wfn[root_idx]                   
                    root_sr_remove.append(root_idx)        
        elif len(sr_wfn) == len(ion_wft):
            print("selection form the ion output")
            for _ in [0]:
                Projection_root={}
                for root_mol in range(len(sr_wfn)):
                    Projection_root[root_mol]=[]
                    for root_ion in range(len(ion_wft)):
                        S_mol_ion=0
                        for (conf_mol,conf_ion) in zip(list(sr_wfn[root_mol].keys()),list(ion_wft[root_ion].keys())):
                                S_mol_ion+=sr_wfn[root_mol][conf_mol] * ion_wft[root_ion][conf_ion]
                        S_mol_ion = S_mol_ion**2
                        Projection_root[root_mol].append(S_mol_ion)
                weight_root={}
                for root in range(len(Projection_root)):
                    weight_root[root]=[]
                    total = sum(Projection_root[root])
                    total_GS=0
                    for GS in range(len(mag_center.basis_lmls)):
                        total_GS += Projection_root[root][GS] 
                    weight_root[root]=total_GS/total
                root_select_by_term=[]
                print(f'Root selected if the projection is upper to {np.round(200/(2*mag_center.l+1),2)}%')
                for root in range(len(weight_root)):
                    if weight_root[root]>2/(2*mag_center.l+1):
                        root_select_by_term.append(root)
    
                if len(root_select_by_term)!=len(mag_center.basis_lmls):
                    print('Internal error during the projection')
                    print(f"We apply the pseudo-L approximation (by selecting the first 2L+1={len(mag_center.basis_lmls)} roots). \n\n")
                    for root_idx in range(len(sr_wfn)):           
                        if root_idx >= len(mag_center.basis_lmls):
                            del sr_wfn[root_idx]                  
                            root_sr_remove.append(root_idx)       
                    break
                all_sr_root=list(range(0,len(sr_wfn)))
                root_sr_remove_term=sorted(list(set(all_sr_root)-set(root_select_by_term)))
                printing=True
                weight_wfn_print={} 
                for root in range(len(weight_root)):
                    if root in root_select_by_term:
                        weight_wfn_print[f'Root {root}']=str(np.round(weight_root[root]*100,2))+' % ---> Selected'
                    else:
                        weight_wfn_print[f'Root {root}']=str(np.round(weight_root[root]*100,2))+' %'
                if printing:                                                         
                    normalisation={}
                    normalisation['Scalar projection']=weight_wfn_print
                    print('Scalar: Projection between ion and molcular wavefunction')
                    print("Projection = Sum_{root}^{2L+1} | Sum_{k} <psi_{molecule}^{k}|psi_{ion}^{k}> |_{root}^{2} ")
                    df=pd.DataFrame(normalisation)
                    print(df,end='\n\n')
                for root in root_sr_remove_term:
                    del sr_wfn[root]
                for key_sr, new_key in zip(root_select_by_term,range(len(mag_center.basis_lmls))):
                    sr_wfn[new_key]=sr_wfn.pop(key_sr)   
                root_sr_remove = root_sr_remove + root_sr_remove_term

    return sr_wfn, root_sr_remove
