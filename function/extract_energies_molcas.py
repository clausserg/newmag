from function.helper_functions import *
import re
import numpy as np
import h5py

def get_energies_molcas(filename, level, root_sr_remove=[],soc_root_remove=[]):
    '''                                                                             
     reading the sr or so energies form molcas output                                  
                                                                                     
     INPUT                                                                           
     filename : main molcas output                                                   
     level : calculation level                                                       
     root_sr_remove : sr root not selected in the model                              
     soc_root_remove : sr root not selected in the model                           
                                                                                     
     OUTPUT                                                                          
     Energies : list of energies as a function of the level (in cm**-1) 
    '''                                                                             


    h5_file=filename.split(".")
    del h5_file[-1]
    h5_file=str('.'.join(h5_file))  
    if level=='casscf-sr':
        try: 
            marker="ROOT_ENERGIES"
            with h5py.File(h5_file+'.rasscf.h5', "r") as f:
                energies_hartree = list(f[marker][()])
            # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
            energies=[]
            for e in range(len(energies_hartree)):
                energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
            energies = np.delete(energies, root_sr_remove, axis=0)

        except FileNotFoundError:
            energies_hartree=[]
            marker_CASCI = False
            pattern = re.compile(r'::\s+RASSCF root number\s+(\d+)\s+Total energy:')
            with open(filename, "r") as f:
                for line in f:
                    if not marker_CASCI:
                        if "CASCI only, no orbital optimization will be done" in line:
                            marker_CASCI = True
                        continue
                    if marker_CASCI:
                        matches = pattern.search(line)
                        if matches: 
                            energies_hartree.append(float(line.split()[-1])) 
                    if " Molecular orbitals:" in line:
                        marker_CASCI = False
            # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
            energies=[]
            for e in range(len(energies_hartree)):
                energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
            energies = np.delete(energies, root_sr_remove, axis=0)

    elif level=='caspt2-sr':
        energies_hartree=[]
        pattern = re.compile(r'::\s+CASPT2 Root\s+(\d+)\s+Total energy:')
        with open(filename, "r") as f:
            for line in f:
                if 'Total CASPT2 energies:' in line:
                    continue
                matches = pattern.search(line)
                if matches:
                    energies_hartree.append(float(line.split()[-1]))    
        # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
        energies=[]
        for e in range(len(energies_hartree)):
            energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
        energies = np.delete(energies, root_sr_remove, axis=0)

    elif level in ['casscf-so', 'caspt2-so']:
        try :
            marker="SOS_ENERGIES"
            with h5py.File(h5_file+'.rassi.h5', "r") as f:
                energies_hartree = np.round(list(f[marker][()]),8)
        
            # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
            energies=[]
            for e in range(len(energies_hartree)):
                energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
            energies = np.delete(energies, soc_root_remove, axis=0)

        except FileNotFoundError:
            energies_hartree=[]
            pattern = re.compile(r'::\s+SO-RASSI State\s+(\d+)\s+Total energy:')
            with open(filename, "r") as f:
                for line in f:
                    matches = pattern.search(line)
                    if matches: 
                        energies_hartree.append(float(line.split()[-1])) 
            # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
            energies=[]
            for e in range(len(energies_hartree)):
                energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
            energies = np.delete(energies, soc_root_remove, axis=0)

    Energies=[]
    for e in range(len(energies)):
        Energies.append(energies[e]-min(energies))
    return Energies
