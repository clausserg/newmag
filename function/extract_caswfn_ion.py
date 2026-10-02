from copy import deepcopy 
import re
import numpy as np        
import math
import pandas as pd       
import h5py               
from function.helper_functions import *

def casscf_wfn_ion(ion_out, mag_center):
    '''
    reading the sr wavefunction form molcas/orca ion output using a minimal CAS
    
    INPUT
    ion_out : ion molcas/orca output 
    mag_center : class of functions possessing all the properties characteristic of the metallic center

    OUTPUT
    sr_wfn : dict[ml][coeff], sr ion wavefunction
    '''
    if software(ion_out) == "molcas":
        pattern_wfn_tmp = generate_patern_particles(2*mag_center.l+1, mag_center.nel)                       
        tmp='|'.join(pattern_wfn_tmp)
        pattern_wfn = ['(', tmp ,')']
        pattern_wfn=''.join(pattern_wfn)

        ml_order = ml_order_molcas(ion_out, mag_center) 
        sr_wfn={}           
        conf_wfn={}         
        sign_conf={}        
        marker_CASCI = False
        with open(ion_out, "r") as f:                                              
            for line in f:
                if not marker_CASCI:
                    if "CASCI only, no orbital optimization will be done" in line:
                        marker_CASCI = True
                    continue
                if "printout of CI-coefficients larger than" in line:
                    root=int(line.split()[-1])-1
                    sr_wfn[root]={}
                    conf_wfn[root]={}
                    sign_conf[root]={}
                    for base in mag_center.basis_ne_ml_uncoupled :
                        sr_wfn[root][base]=[]
                        conf_wfn[root][base]=[]
                        sign_conf[root][base]=[]
                        if len(set(base)) != len(base):
                            sr_wfn[root][base] = np.float64(0)

                if re.search(pattern_wfn,line):   
                    conf=line.split()[1]
                    if len(list(conf)) != 2*mag_center.l+1:
                        error_molcas(mag_center)


                    ml_idx=[ml.span()[0] for ml in re.finditer(re.compile('u'), conf)]
                    
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
                    if mag_center.l == 3:
                        dico_ml_sign = phase_rassi_molcas_f(mag_center.nel, int(line.split()[0]), dico_ml_sign)
                    if mag_center.l == 2:                                                                       
                        dico_ml_sign = phase_rassi_molcas_d(mag_center.nel, int(line.split()[0]), dico_ml_sign)
                    for i,conf in enumerate(dico_ml_sign):
                        conf_wfn[root][conf].append(int(line.split()[0])-1)
                        sign_conf[root][conf].append(dico_ml_sign[conf])
                        sr_wfn[root][conf] = float(line.split()[2]) * dico_ml_sign[conf]


        try:
            marker="CI_VECTORS"           
            h5_file=ion_out.split(".")   
            del h5_file[-1]               
            h5_file=str('.'.join(h5_file))
            with h5py.File(h5_file+'.rasscf.h5', "r") as f:             
                wfn_cas = f[marker][()]
                wfn_cas = wfn_cas.T
                sr_wfn = {}
                for root_idx in range(len(wfn_cas.T)):
                    sr_wfn[root_idx]={}
                    for base in mag_center.basis_ne_ml_uncoupled:
                        sr_wfn[root_idx][base] = wfn_cas[conf_wfn[root_idx][base],root_idx]
                        sr_wfn[root_idx][base] *= sign_conf[root_idx][base]
                        if len(set(base)) != len(base):   
                            sr_wfn[root_idx][base] = int(0)
                        sr_wfn[root_idx][base] = float(sr_wfn[root_idx][base])

        except FileNotFoundError: 
            for root_idx in range(len(sr_wfn)):       
                for base in mag_center.basis_ne_ml_uncoupled:
                    if len(set(base)) != len(base):   
                        sr_wfn[root_idx][base] = int(0)

            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('WARNING:') 
            print('For OpenMolcas, to get the CAS wavefunction, we recommend to use the the h5 file') 
            print('Put the h5 in the same folder of your output we this type name: $Project.rasscf.h5')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')  
    
    if software(ion_out) == "orca":
        pattern_wfn_tmp = generate_patern_particles(2*mag_center.l+1, mag_center.nel)                       
        tmp='|'.join(pattern_wfn_tmp)
        pattern_wfn = ['(', tmp ,')']
        pattern_wfn=''.join(pattern_wfn)

        if mag_center.l==3: ml_order = [0, 1, -1, 2, -2, 3, -3]
        elif mag_center.l==2: ml_order = [0, 1, -1, 2, -2]

        sr_wfn = {} 
        marker_CI=False
        with open(ion_out, 'r') as f:
            for line in f:
                if not marker_CI:
                    if "Spin-Determinant CI Printing" in line:
                        marker_CI=True
                    continue
                if " ROOT" and ":  E=" in line:
                    root=int(line.split()[1].split(':')[0])
                    sr_wfn[root]={}
                    for base in mag_center.basis_ne_ml_uncoupled :
                        sr_wfn[root][base]=[] 
                        if len(set(base)) != len(base):
                            sr_wfn[root][base] = np.float64(0)
                    
                if re.search(pattern_wfn,line):   
                    conf=line.split()[0].split('[')[1].split(']')[0]
                    if len(list(conf)) != 2*mag_center.l+1:
                        print('ERROR in the reading ion output')
                        print(f'This output must be a minimal CAS -> CAS({mag_center.nel},{2*mag_center.l+1})')
                        error_orca(mag_center) 
                    ml_idx=[m.span()[0] for m in re.finditer(re.compile('u'), conf)]
                    
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

        for root_idx in range(len(sr_wfn)):
            for base in mag_center.basis_ne_ml_uncoupled:
                if len(set(base)) != len(base):
                    sr_wfn[root_idx][base] = int(0)       

        ne_ml=mag_center.basis_ne_ml_uncoupled                                           
        if mag_center.l==3:                
            for root in range(len(sr_wfn)):
                for base in range(len(sr_wfn[root])):
                    if (3 in ne_ml[base]) or (-3 in ne_ml[base]):
                        if (3 in ne_ml[base]) and (-3 in ne_ml[base]):
                            continue      
                        sr_wfn[root][ne_ml[base]] *= -1              

    return sr_wfn

def phase_rassi_molcas_f(nel, csf, dico_ml_sign):
    '''
    molcas correction phase
    '''
    if nel in [1,2,3,4,5,6]:
        return dico_ml_sign
    if nel == 8:
        if csf in [2,4,6]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 9:
        if csf in [2,5,7,8,10,13,15,17,20]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] = -1*dico_ml_sign[conf]
    if nel == 10:
        if csf in [2,4,7,9,10,12,15,18,20,22,23,25,27,30,32,34]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 11:
        if csf in [2,4,6,9,11,13,14,16,18,21,24,26,27,29,32,34]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 12:
        if csf in [2,5,7,9,12,14,15,17,20]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 13:
        if csf in [2,4,6]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
 
    return dico_ml_sign

def phase_rassi_molcas_d(nel, csf, dico_ml_sign):
    '''
    molcas correction phase
    '''
    if nel in [1,2,3,4]:
        return dico_ml_sign
    if nel == 6:
        if csf in [1,3]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 7:
        if csf in [1,3,5,8]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 8:
        if csf in [1,4,6,8]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1
    if nel == 9:
        if csf in [1,3]:
            for i,conf in enumerate(dico_ml_sign):
                dico_ml_sign[conf] *= -1



    return dico_ml_sign

def error_molcas(mag_center):
    print('ERROR in the reading ion output')
    print(f'This output must be a minimal CAS -> CAS({mag_center.nel},{2*mag_center.l+1})')
    print('exemple OpenMolcas_v25.06 input:')
    if mag_center.l==3: OA='f'
    if mag_center.l==2: OA='d'
    n_cas = (2*mag_center.S+1)/(2*mag_center.l+2)
    n_cas *= math.comb(2*mag_center.l+2, int(mag_center.nel/2-mag_center.S))
    n_cas *= math.comb(2*mag_center.l+2, int(mag_center.nel/2+mag_center.S+1))
    print(f"""
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
  
""")
 
    print('EXIT Program')
    exit()
    return
 

def error_orca(mag_center):
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
end
 
""")
                                                                                                  
    print('EXIT Program')
    exit()
    return

