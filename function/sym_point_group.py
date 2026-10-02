from pointgroup import PointGroup
import numpy as np
import pandas as pd
import re
from function.helper_functions import *

# difine two dict with all periodic table informations
list_all_atom = ['Q','H','He','Li',
    'Be', 'B'  , 'C'  , 'N', 'O',   
    'F' , 'Ne' , 'Na' , 'Mg' ,
    'Al', 'Si' , 'P'  , 'S'  , 
    'Cl', 'Ar' , 'K'  , 'Ca' ,
    'Sc', 'Ti' , 'V'  , 'Cr' ,
    'Mn', 'Fe' , 'Co' , 'Ni' ,
    'Cu', 'Zn' , 'Ga' , 'Ge' ,
    'As', 'Se' , 'Br' , 'Kr' ,
    'Rb', 'Sr' , 'Y'  , 'Zr' , 
    'Nb', 'Mo' , 'Tc' , 'Ru' ,
    'Rh', 'Pd' , 'Ag' , 'Cd' ,     
    'In', 'Sn' , 'Sb' , 'Te' ,
    'I' , 'Xe' , 'Cs' , 'Ba' ,     
    'La', 'Ce' , 'Pr' , 'Nd' ,    
    'Pm', 'Sm' , 'Eu' , 'Gd' ,
    'Tb', 'Dy' , 'Ho' , 'Er' ,     
    'Tm', 'Yb' , 'Lu' , 'Hf' ,    
    'Ta', 'W'  , 'Re' , 'Os' ,
    'Ir', 'Pt' , 'Au' , 'Hg' ,      
    'Tl', 'Pb' , 'Bi' , 'Po' ,
    'At', 'Rn' , 'Fr' , 'Ra' ,      
    'Ac', 'Th' , 'Pa' , 'U'  ,       
    'Np', 'Pu' , 'Am' , 'Cm']

Ln_An=['Ce' , 'Pr' , 'Nd' , 'Pm',
        'Sm' , 'Eu' , 'Gd' , 'Tb', 
        'Dy' , 'Ho' , 'Er' , 'Tm', 'Yb', 
        'Th' , 'Pa' , 'U'  ,       
        'Np', 'Pu' , 'Am' , 'Cm']
 
#covalent radii from Alvarez (2008)
#DOI: 10.1039/b801115j
covalent_radii = {
    'Q':0.0,'H': 0.31, 'He': 0.28, 'Li': 1.28,
    'Be': 0.96, 'B': 0.84, 'C': 0.76, 
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
    'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
    'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
    'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
    'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
    'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
    'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
    'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
    'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
    'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
    'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
    'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
    'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
    'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
    'Am': 1.80, 'Cm': 1.69                              
}   

def symmetry_part(filename, tolerance_eig=0.1 , tolerance_ang=5):
    """
    determination of symmetry point group of the molecule in the main outpute file
    the determination is done for the gobal symmetry and the local symmetry (metal + first coorination sphere)
    The determination is done via the pointgroup Python package:
    https://github.com/abelcarreras/pointgroup
    
    INPUTE
    - filename : main outpute file given by the user
    - tolerance_eig : tolerance of eigenvalues in the inertial matrice 
    - tolerance_ang : tolerance of the angle (in degree)
    
    OUTPUTE
    - nothing
    """
    print('The point group analysis is done by the PointGroup package')
    print(f'Inertia tensor precision = {tolerance_eig}')
    print(f'Angular tolerance in degrees = {tolerance_ang}')
    print('Warning: the point charge is not taken into account')
    print('https://github.com/abelcarreras/pointgroup',end='\n\n')
    if software(filename) == "orca":
        coords, symbols_atom = read_geom_xyz_orca(filename)
    elif software(filename) == "molcas":
        coords, symbols_atom = read_geom_xyz_molcas(filename)
    
    atom_center=[]
    for atom in symbols_atom:
            if atom in Ln_An:
                   atom_center.append(atom)
    if len(atom_center)>1:
        print('Your molecule contains: ',end='')
        print(",".join(["'" + item + "'" for item in atom_center]))
        print('NewMag can only consider one metal center...')
        print('Carefully review the results!')
        return
    elif len(atom_center)==0:
        print('Your molecule does not contain a metal center... ')
        print('Carefully review the results!')
        return
    elif len(atom_center)==1:
        atom_center=atom_center[0]
        coords_global = np.empty((len(symbols_atom),3), np.float64)
        df_xyz_tmp = pd.DataFrame(columns=['x', 'y', 'z'])
        for i in range(len(symbols_atom)):
                X=coords[i][0]-coords[symbols_atom.index(atom_center)][0]
                Y=coords[i][1]-coords[symbols_atom.index(atom_center)][1]
                Z=coords[i][2]-coords[symbols_atom.index(atom_center)][2]
                df_xyz_tmp.loc[i] = [X] + [Y] + [Z]
        #use the radii to selct the 1er coordination sphere of Ln or An
        coords_global=df_xyz_tmp.to_numpy()
        df_xyz_center = pd.DataFrame(columns=['x', 'y', 'z'])
        symbols_atom_local=[]
        i=0 
        for atom, xyz in zip(symbols_atom,coords_global):
            distance = np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
            if distance <= (covalent_radii[atom] + covalent_radii[atom_center]):
                df_xyz_center.loc[i] = [xyz[0]] + [xyz[1]] + [xyz[2]]
                symbols_atom_local.append(atom)
                i+=1
        coords_local=df_xyz_center.to_numpy()                                       

    try:
            pg_gobal = PointGroup(coords_global, symbols_atom, tolerance_eig, tolerance_ang)
            PG_gobal= pg_gobal.get_point_group()
            print(f'Global point group of your molecule: {PG_gobal}')
    
    except Exception: 
            print('Error determining the global point group')
    
    
    try:
        pg_local = PointGroup(coords_local, symbols_atom_local, tolerance_eig, tolerance_ang)
        PG_local = pg_local.get_point_group()
        print(f'Point group taking into account only the first coordination sphere of {atom_center}: {PG_local}')
        print('The covalent radii value come from Alvarez (2008) DOI: 10.1039/b801115j', end='\n\n')
        print('The following extraction of the crystal field is performed based on the')
        print('orientation of your molecule as shown in your output', end='\n\n')
        if PG_local not in ['C1', 'Ci', 'Cs']:
                print('Your magnetic core has local symmetry')
                print('We suggest recalculating using this orientation:')
                print(30*'-')
                coords_local = pg_local.get_standard_coordinates()
                try:
                    coords_global = pg_gobal.get_standard_coordinates()
                except UnboundLocalError:pass
                df_xyz_tmp = pd.DataFrame(columns=['x', 'y', 'z'])
                for i in range(len(symbols_atom)):
                        X=coords_global[i][0]-coords_global[symbols_atom.index(atom_center)][0]
                        Y=coords_global[i][1]-coords_global[symbols_atom.index(atom_center)][1]
                        Z=coords_global[i][2]-coords_global[symbols_atom.index(atom_center)][2]
                        df_xyz_tmp.loc[i] = [X] + [Y] + [Z]
        
                coords_global=df_xyz_tmp.to_numpy()
                pg_local = PointGroup(coords_local, symbols_atom_local, tolerance_eig, tolerance_ang)
                Rot=np.array(pg_local.get_principal_axis_of_inertia())
                new_coords_gobal = np.array(coords_global @ Rot.T)
                new_coords_local = np.array(coords_local @ Rot.T)
                I=pg_local.get_principal_moments_of_inertia()
        
                print(len(symbols_atom),end='\n\n')
                if PG_local in ['Cinfv', 'Dinfh']:
                        for atom, xyz in zip(symbols_atom,new_coords_gobal):
                                print("{} \t {:11.8f} \t {:11.8f} \t {:11.8f}".format(atom,xyz[2],xyz[1],xyz[0]))
                elif PG_local not in ['O','Oh','T','Td','Th','Ih']:
                        new_coords_local=np.array(new_coords_local)
                        dico_conf={
                                0:[0,1,2],
                                1:[0,2,1],
                                2:[1,0,2],
                                3:[1,2,0],
                                4:[2,0,1],
                                5:[2,1,0]}
                        n=int(list(PG_local)[1])
                        for conf in range(len(dico_conf)):
                                test_coords = np.empty((len(symbols_atom_local),3), np.float64)
                                test_coords[:,[0,1,2]] = new_coords_local[:,dico_conf[conf]]
                                Rz=check_Rz(n, test_coords, symbols_atom, tolerance_eig, tolerance_ang)
                                if Rz:
                                        for atom, xyz in zip(symbols_atom,new_coords_gobal):
                                                X=xyz[dico_conf[conf][0]]
                                                Y=xyz[dico_conf[conf][1]]
                                                Z=xyz[dico_conf[conf][2]]
                                                print("{} \t {:11.8f} \t {:11.8f} \t {:11.8f}".format(atom,X,Y,Z))
                                        break
                                elif conf==5:
                                    print("Internal Error. continue")
                else: 
                        for atom, xyz in zip(symbols_atom,new_coords_gobal):
                                print("{} \t {:11.8f} \t {:11.8f} \t {:11.8f}".format(atom,xyz[0],xyz[2],xyz[1]))

                print(30*'-')
        elif PG_local in ['C1', 'Ci', 'Cs']:
            print('Your magnetic core has no local symmetry')
            print('Be careful to use an orientation that makes physical sense')
    except Exception:
            print('Error determining the global point group')
            print('The following extraction of the crystal field is performed based on the')
            print('orientation of your molecule as shown in your output', end='\n\n')


    print('',end='\n\n')
    print('If you need to modify the xyz coordinates, we recommend using the xyzalign script')
    print('https://github.com/radi0sus/xyzalign')
    print('',end='\n\n')
    return


def check_Rz(n, coords, symbols_atom,  tolerance_eig, tolerance_ang):
        '''
        check if the R_z (C_n) operator exists 
        
        INPUT
        - n : integer, order of the rotation
        - coords : molecule coordinates
        - symbols_atom : list of atoms
        - tolerance_eig : tolerance of eigenvalues in the inertial matrice 
        - tolerance_ang : tolerance of the angle (in degree)
        
        OUTPUT
        - True or False
        '''

        Rz = np.array([
        [np.cos(2*np.pi/n), -np.sin(2*np.pi/n), 0],
        [np.sin(2*np.pi/n),  np.cos(2*np.pi/n), 0],
        [0,                                 0,  1]])
        error_abs_rad = tolerance_eig / np.clip(np.linalg.norm(coords, axis=1), tolerance_eig,None)

        Rz_coords = np.dot(coords, Rz)
        for idx, Rz_coord in enumerate(Rz_coords):

            norm_coor = np.linalg.norm(coords, axis=1)
            norm_Rz_coor = np.linalg.norm(Rz_coord)
            average_radii = np.clip((norm_coor + norm_Rz_coor) / 2, tolerance_eig, None)
            difference_rad = np.abs(norm_coor - norm_Rz_coor) / average_radii

            norm_coor = np.linalg.norm(coords, axis=1)
            norm_Rz_coor = np.linalg.norm(Rz_coord)
            angles = []
            for v, n in zip(np.dot(Rz_coord, coords.T), norm_coor*norm_Rz_coor):
                if n < tolerance_eig:
                    angles.append(0)
                else:
                    angles.append(np.arccos(np.clip(v/n, -1.0, 1.0)))
            difference_ang = np.array(angles)

            def check_diff(diff, diff2):
                for idx_2, (d1, d2) in enumerate(zip(diff, diff2)):
                    if symbols_atom[idx_2] != symbols_atom[idx]:
                        continue
                    # d_r = np.linalg.norm([d1, d2])
                    tolerance_total = tolerance_ang * tolerance_eig + error_abs_rad[idx_2]
                    if d1 < tolerance_total and d2 < tolerance_total:
                        return True
                return False

            if not check_diff(difference_ang, difference_rad):
                return False
        return True

def read_geom_xyz_molcas(molcas_out):
    '''
    reading the coordinate of the molecule for molcas outpute format
    
    INPUTE 
    - molcas_out : main outpute file given by the user
    
    OUTPUT
    - coords : array[x,y,z] with the coordinate
    - symbols_atom : list of the atom in the same order of the coords array
    '''

    symbols_atom=[]
    df_xyz = pd.DataFrame(columns=['x', 'y', 'z'])
    j=0
    with open(molcas_out, "r") as file:
        lines = file.readlines()
    inside_block_1 = False
    inside_block_2 = False

    for line in lines:
        if "Cartesian coordinates in angstrom:" in line :
            inside_block_1 = True
            continue
        if "No.  Label" in line :
            inside_block_2 = True
            continue
        if inside_block_1 and inside_block_2:
            if 'Nuclear repulsion energy =' in line:
                break
            try:
                
                if re.search(r"\b[A-Z]+\d*\b", line):
                    symbol_tmp=list(line.split()[1])
                    symbol_tmp_2=[]
                    for i in range(len(symbol_tmp)):
                        try :
                            int(symbol_tmp[i])
                            #symbol_tmp.pop(i)
                        except ValueError:
                            if i==0:
                                symbol_tmp_2.append(symbol_tmp[i].upper())
                            else:
                                symbol_tmp_2.append(symbol_tmp[i].lower())
                    symbols_atom.append("".join(symbol_tmp_2))
                    df_xyz.loc[j] = [float(line.split()[2])] + [float(line.split()[3])] + [float(line.split()[4])]
                    j+=1
            except IndexError: pass
    coords=df_xyz.to_numpy()
        
    return coords, symbols_atom

def read_geom_xyz_orca(orca_out):
    '''
    reading the coordinate of the molecule for orca outpute format
    
    INPUTE 
    - orca_out : main outpute file given by the user
    
    OUTPUT
    - coords : array[x,y,z] with the coordinate
    - symbols_atom : list of the atom in the same order of the coords array
    '''


    symbols_atom=[]
    df_xyz = pd.DataFrame(columns=['x', 'y', 'z'])
    i=0
    with open(orca_out, "r") as file:
        lines = file.readlines()
    inside_block = False
    for line in lines:
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line :
            inside_block = True
            continue
        if inside_block:
            if 'CARTESIAN COORDINATES (A.U.)' in line:
                break
            try:
                if str(line.split()[0]) in list(list_all_atom):
                    #print(line.strip()[0])
                    symbols_atom.append(str(line.split()[0]))
                    df_xyz.loc[i] = [float(line.split()[1])] + [float(line.split()[2])] + [float(line.split()[3])]
                    i+=1
            except IndexError: pass
    coords=df_xyz.to_numpy()

    return coords, symbols_atom
