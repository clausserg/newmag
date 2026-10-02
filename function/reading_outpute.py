from function.extract_caswfn_molcas import *
from function.extract_caswfn_orca import *
from function.extract_sociwfn_molcas import *
from function.extract_sociwfn_orca import *
from function.extract_energies_molcas import *
from function.extract_energies_orca import *
from function.extract_caswfn_ion import *
from function.extract_aiLFT_orca import *
import copy

def reading_outpute(out, mag, ion_out, cas_mo, user_sr_root, user_so_root, user_ml_order):
    '''
    read the all importante information in the outpute file given by the user
    
    INPUTE
    - out : main outpute file given by the user
    - mag : class of functions possessing all the properties characteristic of the metallic center
    - ion_out : ion outpute file given by the user
    - cas_mo : formating of the CSF
    - user_sr_root : list of sr root selected by user
    - user_so_root : list of so root selected by user
    - user_ml_order : list of ml in the main output that make up the CSF
    
    OUTPUE
    - wavefunction : dict[level][ml][coeff], dictionary with the wavefunction
    - energies : dict[level][energies], dictionary with the energies of the wavefunction
    - list_calc_level : list of the calculation level
    - Heff : dict[level][array of Heff], dictionary with each effectif hamiltonienne
    '''
    if ion_out!=None:
        ion_wft = casscf_wfn_ion(ion_out, mag)
    else: ion_wft=None
            
    wavefunction, energies = {}, {}
    Heff={} 
    list_calc_level = []
    pattern_wfn=[]
    ml_order=[] 
    if software(out) == "orca":
        pattern_wfn, remove_conf=pattern_caswfn(cas_mo, mag, "orca")
        if user_ml_order == None:
            if mag.l==3: ml_order = [0, 1, -1, 2, -2, 3, -3]
            elif mag.l==2: ml_order = [0, 1, -1, 2, -2]
            elif mag.l==1: ml_order = [0, 1, -1]
        else : ml_order = user_ml_order 
        print('--- CASSCF LEVEL ---')
        wavefunction['casscf-sr'], root_sr_remove = casscf_wfn_orca(out, mag, pattern_wfn, ml_order, remove_conf, ion_wft, user_sr_root)
        wavefunction['casscf-so'], soc_root_remove = soci_wfn_orca(out, 'casscf-so', mag, root_sr_remove, user_so_root) 
        energies['casscf-sr'] = get_energies_orca(out, 'casscf-sr', root_sr_remove, soc_root_remove)
        energies['casscf-so'] = get_energies_orca(out, 'casscf-so', root_sr_remove, soc_root_remove)
        list_calc_level.append('casscf')
            
        if nevpt2_orca_present(out):
            print('--- NEVPT2 LEVEL ---')
            wavefunction['nevpt2-sr'], root_sr_remove = casscf_wfn_orca(out, mag, pattern_wfn, ml_order, remove_conf, ion_wft, user_sr_root)
            wavefunction['nevpt2-so'], soc_root_remove = soci_wfn_orca(out, 'nevpt2-so', mag, root_sr_remove, user_so_root) 
            energies['nevpt2-sr'] = get_energies_orca(out, 'nevpt2-sr', root_sr_remove, soc_root_remove)
            energies['nevpt2-so'] = get_energies_orca(out, 'nevpt2-so', root_sr_remove, soc_root_remove)
            list_calc_level.append('nevpt2')
            
        ailft_present, ailft_level = ailft_in_orca(out)
        if ailft_present:
            if 'casscf' in ailft_level:
                Heff['casscf-ailft'] = H_ailft_orca(out,mag, 'casscf') 
            if 'nevpt2' in ailft_level:
                Heff['nevpt2-ailft'] = H_ailft_orca(out,mag, 'nevpt2') 
            
    elif software(out) == "molcas":
        pattern_wfn, remove_conf=pattern_caswfn(cas_mo, mag, "molcas")   
        if user_ml_order == None:
            ml_order = ml_order_molcas(out, mag)
        else : ml_order = user_ml_order 
        if caspt2_molcas_present(out):
            if type_caspt2_molcas(out)=='ms-caspt2-sr':
                print('MS-CASPT2 is not implemented...')
                print('And never will be :)')
                general_error_inp_molcas(mag)

            if type_caspt2_molcas(out)=='ss-caspt2-sr':
                print('--- CASSCF LEVEL ---')
                wavefunction['casscf-sr'], root_sr_remove = casscf_wfn_molcas(out, mag, pattern_wfn, ml_order, remove_conf, ion_wft, user_sr_root)
                print('--- CASPT2 LEVEL ---')
                wavefunction['caspt2-sr'], root_sr_remove = casscf_wfn_molcas(out, mag, pattern_wfn, ml_order, remove_conf, ion_wft, user_sr_root)
                wavefunction['casscf-so']=None
                wavefunction['caspt2-so'], soc_root_remove = soci_wfn_molcas(out, 'caspt2-so', mag, root_sr_remove, user_so_root)
                energies['casscf-sr'] = get_energies_molcas(out, 'casscf-sr', root_sr_remove, soc_root_remove)                        
                energies['caspt2-sr'] = get_energies_molcas(out, 'caspt2-sr', root_sr_remove, soc_root_remove)
                energies['caspt2-so'] = get_energies_molcas(out, 'caspt2-so', root_sr_remove, soc_root_remove)
                energies['casscf-so'] = None
                list_calc_level.append('casscf')
                list_calc_level.append('caspt2')
            else:
                 print("There are one or more ‘CASPT2’ iterations in the molcas output file, but no actual CASPT2 calculation...")
                 print("We will stop the program here to be sure we do not produce bullshit :) ")
                 general_error_inp_molcas(mag)
        else:
             print('--- CASSCF LEVEL ---')
             wavefunction['casscf-sr'], root_sr_remove = casscf_wfn_molcas(out, mag, pattern_wfn, ml_order, remove_conf, ion_wft, user_sr_root)
             wavefunction['casscf-so'], soc_root_remove = soci_wfn_molcas(out, 'casscf-so', mag, root_sr_remove, user_so_root)
             energies['casscf-sr'] = get_energies_molcas(out, 'casscf-sr', root_sr_remove, soc_root_remove)
             energies['casscf-so'] = get_energies_molcas(out, 'casscf-so', root_sr_remove, soc_root_remove)
             list_calc_level.append('casscf')

    return wavefunction, energies, list_calc_level, Heff
