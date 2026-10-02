import re
import numpy as np  


def get_energies_orca(orca_out, level="casscf-sr", root_sr_remove=[], soc_root_remove=[]):
    '''                                                                             
     reading the sr or so energies form orca output                                  
                                                                                     
     INPUT                                                                           
     filename : main orca output                                                   
     level : calculation level                                                       
     root_sr_remove : sr root not selected in the model                              
     soc_root_remove : sr root not selected in the model                           
                                                                                     
     OUTPUT                                                                          
     Energies : list of energies as a function of the level (in cm**-1) 
    ''' 

    level = level.lower()
    if level not in {"casscf-sr", "nevpt2-sr", "casscf-so", "nevpt2-so"}:
        raise ValueError("level should be 'casscf-sr', 'nevpt2-sr', 'casscf-so' or 'nevpt2-so'")

    energies = []
    energies_hartree = []

    with open(orca_out, "r") as file:
        lines = file.readlines()

    if level == "casscf-sr":
        inside_block = False
        for line in lines:
            if "Spin-Determinant CI Printing" in line:
                inside_block = True
            if inside_block:
                #if "DENSITY MATRIX" or 'SA-CASSCF TRANSITION ENERGIES' in line:
                #    break
                if 'ROOT' and ":  E=" in line:
                    row = line.split()
                    energies_hartree.append(np.float64(row[3]))
            if "CAS-SCF STATES FOR BLOCK" in line:
                inside_block=False

        # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
        for e in range(len(energies_hartree)):
            energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
        energies = np.delete(energies, root_sr_remove, axis=0)

    elif level == "nevpt2-sr":
        inside_block = False
        for line in lines:
            if line.strip().startswith("NEVPT2 TOTAL ENERGIES"):
                inside_block = True
                continue
            if inside_block:
                if line.strip().startswith("NEVPT2 TRANSITION ENERGIES"):
                    break
                if "EDIAG" in line:
                    row = line.split()
                    energies_hartree.append(np.float64(row[-1]))
         # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
        for e in range(len(energies_hartree)):
            energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)
        energies = np.delete(energies, root_sr_remove, axis=0)

    elif level == "casscf-so":
        inside_block = False
        for line in lines:
            if "QDPT WITH CASSCF DIAGONAL ENERGIES" in line:
                inside_block = True
                continue
            if inside_block:
                if "COMPUTING QDPT PROPERTIE" in line or "SOC CORRECTED" in line:
                    break
                if "STATE" in line:
                    row = line.split()
                    energies.append(np.float64(row[-1]))
        energies = np.delete(energies, soc_root_remove, axis=0)

                       
    elif level == "nevpt2-so":
        inside_block = False
        for line in lines:
            if "QDPT WITH NEVPT2 DIAGONAL ENERGIES" in line:
                inside_block = True
                continue
            if inside_block:
                if "COMPUTING QDPT PROPERTIE" in line or "SOC CORRECTED" in line:
                    break
                if "STATE" in line:
                    row = line.split()
                    energies.append(np.float64(row[-1]))
        energies = np.delete(energies, soc_root_remove, axis=0)

    return energies
