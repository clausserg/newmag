from function.helper_functions import *  
from function.class_function import *
from itertools import permutations
from copy import deepcopy
import numpy as np
import pandas as pd

def H_ailft_orca(orca_output, mag, ailft_level):
    '''                                                                             
     reading the aiLFT matrix form orca out                                 
                                                                                     
     INPUT                                                                           
     orca_output : main orca output                                                   
     mag : class of functions possessing all the properties characteristic of the metallic center  
     ailft_level : calculation level of aiLFT block                                                       
                                                                                     
     OUTPUT                                        
     ailft_matrix : ailft matrix in the |l,ml> representation  
    ''' 
    ailft_matrix = []
    if ailft_level == 'casscf':
        with open(orca_output, "r") as f:
            Ailft_maker = False  
            for line in f:  
                if line.strip() == "AILFT MATRIX ELEMENTS (CASSCF)": 
                    Ailft_maker = True 
                elif line.strip() == "Slater-Condon Parameters (electronic repulsion) :":
                     if Ailft_maker:
                        break
                elif Ailft_maker: 
                    try: nums = [num for num in line.split()] 
                    except ValueError : continue
                    if len(nums) > 1 and  '-----' not in line and 'Orbital' not in line and 'Ligand field' not in line:  
                        ailft_matrix.append(nums[1:]) 

    elif ailft_level == 'nevpt2':
        with open(orca_output, "r") as f:
            Ailft_maker = False  
            for line in f:  
                if line.strip() == "AILFT MATRIX ELEMENTS (NEVPT2)": 
                    Ailft_maker = True 
                elif line.strip() == "Slater-Condon Parameters (electronic repulsion) :":
                    if Ailft_maker:
                        break
                elif Ailft_maker: 
                    try: nums = [num for num in line.split()] 
                    except ValueError : continue
                    if len(nums) > 1 and  '-----' not in line and 'Orbital' not in line and 'Ligand field' not in line:  
                        ailft_matrix.append(nums[1:]) 
 
    for i in range(len(ailft_matrix)):
        for j in range(len(ailft_matrix[i])):
            ailft_matrix[i][j] = float(ailft_matrix[i][j])*219474.63136314
    diag = np.diag(ailft_matrix)-np.mean(np.diag(ailft_matrix))
    for i in range(len(ailft_matrix)):
                ailft_matrix[i][i] = diag[i] 

    ailft_matrix = np.array(ailft_matrix)
    if len(ailft_matrix) == 0:
        return np.eye(2*mag.l+1)
    if len(ailft_matrix)==7:
        # inverse col
        col_temp = np.copy(ailft_matrix[:, 0])
        ailft_matrix[:, 0] = ailft_matrix[:, 6]
        ailft_matrix[:, 6] = col_temp

        col_temp = np.copy(ailft_matrix[:, 5])
        ailft_matrix[:, 5] = ailft_matrix[:, 6]
        ailft_matrix[:, 6] = col_temp

        col_temp = np.copy(ailft_matrix[:, 1])
        ailft_matrix[:, 1] = ailft_matrix[:, 4]
        ailft_matrix[:, 4] = col_temp
 
        col_temp = np.copy(ailft_matrix[:, 5])
        ailft_matrix[:, 5] = ailft_matrix[:, 3]
        ailft_matrix[:, 3] = col_temp
        
        # inverse line
        line_temp = np.copy(ailft_matrix[0, :])
        ailft_matrix[0, :] = ailft_matrix[6, :]
        ailft_matrix[6, :] = line_temp

        line_temp = np.copy(ailft_matrix[5, :])
        ailft_matrix[5, :] = ailft_matrix[6, :]
        ailft_matrix[6, :] = line_temp

        line_temp = np.copy(ailft_matrix[1, :])
        ailft_matrix[1, :] = ailft_matrix[4, :]
        ailft_matrix[4, :] = line_temp
 
        line_temp = np.copy(ailft_matrix[5, :])
        ailft_matrix[5, :] = ailft_matrix[3, :]
        ailft_matrix[3, :] = line_temp
        
        #correction phase orca
        ailft_matrix[:, 0] *=  -1
        ailft_matrix[:, 6] *=  -1
        ailft_matrix[0, :] *=  -1
        ailft_matrix[6, :] *=  -1
 
    if len(ailft_matrix)==5:
        # inverse col
        col_temp = np.copy(ailft_matrix[:, 0])
        ailft_matrix[:, 0] = ailft_matrix[:, 4]
        ailft_matrix[:, 4] = col_temp

        col_temp = np.copy(ailft_matrix[:, 1])
        ailft_matrix[:, 1] = ailft_matrix[:, 2]
        ailft_matrix[:, 2] = col_temp

        col_temp = np.copy(ailft_matrix[:, 2])
        ailft_matrix[:, 2] = ailft_matrix[:, 4]
        ailft_matrix[:, 4] = col_temp
 
        col_temp = np.copy(ailft_matrix[:, 3])
        ailft_matrix[:, 3] = ailft_matrix[:, 4]
        ailft_matrix[:, 4] = col_temp
        
        # inverse line
        line_temp = np.copy(ailft_matrix[0, :])
        ailft_matrix[0, :] = ailft_matrix[4, :]
        ailft_matrix[4, :] = line_temp

        line_temp = np.copy(ailft_matrix[1, :])
        ailft_matrix[1, :] = ailft_matrix[2, :]
        ailft_matrix[2, :] = line_temp

        line_temp = np.copy(ailft_matrix[2, :])
        ailft_matrix[2, :] = ailft_matrix[4, :]
        ailft_matrix[4, :] = line_temp
 
        line_temp = np.copy(ailft_matrix[3, :])
        ailft_matrix[3, :] = ailft_matrix[4, :]
        ailft_matrix[4, :] = line_temp
        
    if len(ailft_matrix)==3:
        # inverse col
        col_temp = np.copy(ailft_matrix[:, 0])
        ailft_matrix[:, 0] = ailft_matrix[:, 1]
        ailft_matrix[:, 1] = col_temp

        col_temp = np.copy(ailft_matrix[:, 0])
        ailft_matrix[:, 0] = ailft_matrix[:, 2]
        ailft_matrix[:, 2] = col_temp

        # inverse line
        line_temp = np.copy(ailft_matrix[0, :])
        ailft_matrix[0, :] = ailft_matrix[1, :]
        ailft_matrix[1, :] = line_temp

        line_temp = np.copy(ailft_matrix[0, :])
        ailft_matrix[0, :] = ailft_matrix[2, :]
        ailft_matrix[2, :] = line_temp

    U1 = harmonics_re_to_sph(MagneticCenter(mag.l, nel=1).basis_lmls)
    ailft_matrix = U1.T.conj() @ ailft_matrix @ U1
    ailft_matrix = ailft_matrix.conj()

    ######################
    # PRINT AILFT MATRIX #
    ######################

    #pretty_matrix(ailft_matrix, MagneticCenter(mag.l, nel=1).basis_lmls , MagneticCenter(3, nel=1), 2)
    #eigvals, eigvecs = np.linalg.eigh(ailft_matrix)
    #ailft_level = ailft_level+'-sr'
    #eigvals = eigvals - min(eigvals)
    #pretty_composition(eigvals, ailft_matrix, ailft_level, MagneticCenter(mag.l, nel=1), 2,False)
    
    return ailft_matrix
