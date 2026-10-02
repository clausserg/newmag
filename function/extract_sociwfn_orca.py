from function.helper_functions import *
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction


def soci_wfn_orca(output, level, mag_center, root_sr_remove, user_so_root):
    '''                                                                             
     reading the so wavefunction form orca output                                  
                                                                                     
     INPUT                                                                           
     filename : main orca output                                                   
     level : calculation level                                                       
     mag_center : class of functions possessing all the properties characteristic of the metallic center
     root_sr_remove : sr root not selected in the model                              
     user_so_root : so root selected by the user                                     
                                                                                     
     OUTPUT                                                                          
     so_wfn : dict[so_root][(str(sr_root), str(S), str(Ms))[coeff], so wavefunction                               
     soc_root_remove : so root not selected in the model                             
    '''                                                                             

    level = level.lower()
    if level not in {"casscf-so", "nevpt2-so"}:
        raise ValueError("level should be 'casscf-so' or 'nevpt2-so'")
        
    with open(output, "r") as file:
        content = file.readlines()

    wavefunctions = {}
    state_index = None
    current_wfn = defaultdict(complex)

    # Map level to keyword in trigger line
    level_map = {
        "casscf-so": "CASSCF",
        "nevpt2-so": "NEVPT2",
    }

    keyword = level_map.get(level.lower(), "CASSCF")  # default fallback
    trigger = f"QDPT WITH {keyword} DIAGONAL ENERGIES"
    
    parsing = False
    
    for line in content:
        # Wait until trigger line appears
        if not parsing:
            if trigger in line:
                parsing = True
            continue  # skip lines until trigger is found
        
        if " STATE" in line and ":" in line:
            # When a new state starts, save the previous one
            if state_index is not None and current_wfn:
                wavefunctions[state_index] = dict(current_wfn)
                current_wfn.clear()

            # Extract state index
            try:
                state_index = int(line.split()[1].strip(":"))
            except (IndexError, ValueError):
                state_index = None

        elif state_index is not None:
            row = line.split()
            if len(row) == 8 and row[4] == "0":
                key = (row[5], row[6], row[7])
                coef = float(row[1]) + float(row[2]) * 1j
                coef = np.complex128(np.float64(row[1])+np.float64(row[2])* 1j)
                current_wfn[key] += np.complex128(coef)

        # Stop when passing known end blocks
        if "Center of nuclear charge           = (" in line or "COMPUTING QDPT PROPERTIES" in line:
            if state_index is not None and current_wfn:
                wavefunctions[state_index] = dict(current_wfn)
            break
        
    all_soc_wfn={}
    if len(wavefunctions)==0:
        return None, [] 

    for root, compo in wavefunctions.items():
        all_soc_wfn[root] = dict(
            sorted(
                compo.items(),
                key=lambda item: float(Fraction(item[0][2]))
            )
        )
    Ms_list = [-(mag_center.S - idx) for idx in range(int(2*mag_center.S + 1))]
    root_cas=number_of_CI_CAS_orca(output, mag_center)
    all_sr_root=list(set(range(0,root_cas+1)))
    sr_root_select=sorted(list(set(all_sr_root)-set(root_sr_remove)))

    for soc_root in list(range(0,  len(all_soc_wfn), 1)):           
        for sr_root in all_sr_root:      
            for Ms in Ms_list:                                 
                if not (str(sr_root), str(mag_center.S), str(Ms)) in all_soc_wfn[soc_root]:
                    all_soc_wfn[soc_root][(str(sr_root), str(mag_center.S),str(Ms))]=np.complex128(0) 
    
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
