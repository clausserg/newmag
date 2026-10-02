from function.class_function import *
from function.helper_functions import *

def input_check(out, ion_out, cfg, decimal, user_sr_root, user_so_root, user_ml_order, cas_mo):
    '''
    Check the input paramters given by the user with someinternal consistency test
    
    INPUTE
    - out : main outpute file given by the user
    - ion_out : ion outpute file given by the user
    - cfg : electronic congiguration of the metal -> f1, f4, d4,...
    - decimal : number of decimal printed in the NewMag out
    - user_sr_root : list of sr root selected by user
    - user_so_root : list of so root selected by user
    - user_ml_order : list of ml in the main output that make up the CSF
    - cas_mo : formating of the CSF
    
    OUTPUTE
    - mag : class of functions possessing all the properties characteristic of the metallic center 
    - cas_mo, user_sr_root, user_so_root : the same information of the INPUTE but another formating 
    
    '''
    # Print a beautiful NewMag logo ! :)
    print(r"""
__/\\\\\_____/\\\__________________________________          
 _\/\\\\\\___\/\\\__________________________________         
  _\/\\\/\\\__\/\\\__________________________________        
   _\/\\\//\\\_\/\\\_____/\\\\\\\\___/\\____/\\___/\\_       
    _\/\\\\//\\\\/\\\___/\\\/////\\\_\/\\\__/\\\\_/\\\_      
     _\/\\\_\//\\\/\\\__/\\\\\\\\\\\__\//\\\/\\\\\/\\\__     
      _\/\\\__\//\\\\\\_\//\\///////____\//\\\\\/\\\\\___    
       _\/\\\___\//\\\\\__\//\\\\\\\\\\___\//\\\\//\\\____   
        _\///_____\/////____\//////////_____\///__\///_____  
                                    __/\\\\____________/\\\\_____________________________        
                                     _\/\\\\\\________/\\\\\\_____________________________       
                                      _\/\\\//\\\____/\\\//\\\__________________/\\\\\\\\__      
                                       _\/\\\\///\\\/\\\/_\/\\\__/\\\\\\\\\_____/\\\////\\\_     
                                        _\/\\\__\///\\\/___\/\\\_\////////\\\___\//\\\\\\\\\_    
                                         _\/\\\____\///_____\/\\\___/\\\\\\\\\\___\///////\\\_   
      Version: 1.3                        _\/\\\_____________\/\\\__/\\\/////\\\___/\\_____\\\_  
      Date: 2026                           _\/\\\_____________\/\\\_\//\\\\\\\\/\\_\//\\\\\\\\__ 
                                            _\///______________\///___\////////\//___\////////___""")
    print(106*'_')
    print('\n')
    orb_dict = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
    if cfg[0] not in orb_dict.keys():
        print('Invalid configuration') 
        print('The orbital can only be: p d or f')
        print('example 1: -cfg f1')
        print('Here the configuration is the f type orbital with 1 electron')
        print('example 2: -cfg f8')
        print('Here the configuration is the f type orbital with 8 electron')
        print('example 3: -cfg d2')
        print('Here the configuration is the d type orbital with 2 electron')
        print('EXIT Program')
        exit()
    if cfg[0] in ['s']:
        print(f'The s orbitals are not implemented')
        print('EXIT Program')
        exit()
    if cfg[0] in ['p']:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('WARNING:') 
        print('p orbitals have been implemented but the code has not been tested for this. ')
        print('Results should be interpreted with caution at your own risk.')  
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')


    shell = orb_dict[cfg[0]]
    electrons = int(cfg[1:])
    
    if cfg=='f7':
        print('For the configuration f7, you should use a ZFS Hamiltonian instead')
        print('Unfortunately, this procedure has not yet been implemented in NewMag')
        print("Let's stop here... :)")
        print('EXIT Program')
        exit()
    elif shell==3:
        if electrons not in [1,2,3,4,5,6,8,9,10,11,12,13]:
            print('Invalid configuration')
            print('You have a f type orbital')
            print('Therefore the number of electrons has to be:')
            print('1 2 3 4 5 6 7 8 9 10 11 12 13')
            print('example 1: -cfg f1')               
            print('Here the configuration is the f type orbital with 1 electron')
            print('example 2: -cfg f8')               
            print('Here the configuration is the f type orbital with 8 electrons')
            print('EXIT Program') 
            exit()

    elif cfg=='d5':
        print('For the configuration d5, you shoud use a ZFS Hamiltonian instead')
        print('Unfortunately, this procedure has not yet been implemented in NewMag')
        print("Let's stop here... :)")
        print('EXIT Program')
        exit()
 
    elif shell==2:
        if electrons not in [1,2,3,4,5,6,7,8,9]:
            print('Invalid configuration')
            print('You have a d type orbital')
            print('Therefore the number of electrons has to be:')
            print('1 2 3 4 5 6 7 8 9')
            print('example 1: -cfg d1')               
            print('Here the configuration is the d type orbital with 1 electron')
            print('example 2: -cfg d8')               
            print('Here the configuration is the d type orbital with 8 electrons')
            print('EXIT Program') 
            exit()

    elif cfg=='p3':
        print('For the configuration p3, you shoud use a ZFS Hamiltonian instead')
        print('Unfortunately, this procedure has not yet been implemented in NewMag')
        print("Let's stop here... :)")
        print('EXIT Program')
        exit()
 
    elif shell==1:
        if electrons not in [1,2,4,5]:
            print('Invalid configuration')
            print('You have a p type orbital')
            print('Therefore the number of electron has to be:')
            print('1 2 4 5')
            print('example 1: -cfg p1')               
            print('Here the configuration is the p type orbital with 1 electron')
            print('example 2: -cfg p4')               
            print('Here the configuration is the p type orbital with 4 electrons')
            print('EXIT Program') 
            exit()
 
    mag = MagneticCenter(l=shell, nel=electrons)

    if decimal<0:
        print('Invalid input')
        print('The number of decimals must be a positive integer')
        print('example: -d 2')
        print('It is the defaut value of this variable')
        print('EXIT Program') 
        exit()

    if user_ml_order != None:
        if len(set(user_ml_order)) != 2*mag.l+1:
            print('Invalid ml_order input')
            print(f'It must have {2*mag.l+1} distinct elements')
            print('example for f type orbital: -ml -2 1 0 -3 3 -1 2') 
            print('example for d type orbital: -ml -2 1 0 -1 2') 
            print('EXIT Program') 
            exit()
        else:
            if shell==3:
                if sorted(set(user_ml_order)) != [-3, -2, -1, 0, 1, 2, 3]:
                    print('Invalid ml_order input')
                    print('The ml values must be in the list: [-3, -2, -1, 0, 1, 2, 3]')
                    print('example for f type orbital: -ml -2 1 0 -3 3 -1 2') 
                    print('EXIT Program')
                    exit()               
            elif shell==2:
                if sorted(set(user_ml_order)) != [-2, -1, 0, 1, 2]:
                    print('Invalid ml_order input')
                    print('The ml values must be in the list: [-2, -1, 0, 1, 2]')
                    print('example for d type orbital: -ml -2 1 0 -1 2') 
                    print('EXIT Program')
                    exit()               

    if user_sr_root != None:
        if len(set(user_sr_root)) != len(mag.basis_lmls):
            print('Invalid selected sr_root input')
            print(f'It must contain {len(mag.basis_lmls)} distinct elements')
            print('EXIT Program') 
            exit()      
        else:
            user_sr_root = sorted(set(user_sr_root))

    if user_so_root != None:
        if len(set(user_so_root)) != len(mag.basis_jmj):
            print('Invalid selected so_root input')
            print(f'It must contain {len(mag.basis_jmj)} distinct elements')
            print('EXIT Program') 
            exit()      
        else:
            user_so_root = sorted(set(user_so_root))
    
    if software(out) not in ["molcas", "orca"]:
        print('Invalid main output')
        print('NewMag works only with OpenMolcas and ORCA...')
        print('Ideally: OpenMolcas_v25.06 and ORCA_v6.1')
        print('EXIT Program')
        exit()
    if ion_out != None:
        if software(ion_out) not in ["molcas", "orca"]:
            print('Invalid ion output')
            print('NewMag works only with OpenMolcas and ORCA...')
            print('Ideally: OpenMolcas_v25.06 and ORCA_v6.1')
            print('EXIT Program')
            exit()
    
    if cas_mo != None:  
        mo_c=0
        for mo in cas_mo:
            nb_type_mo = list(mo)
            if nb_type_mo[-1] not in ['d','c','v','o']:
                print('Invalid cas_mo input')
                print('There are at most 4 categories of MOs in the active space and the sum of (d+c+o+v) must match the number of active orbitals:')
                print('d = doubly occupied -> Orbitals with a strict occupancy of 2 ("inactive")')
                print(f'c = centered orbitals -> Orbitals at the magnetic center used to build the model space (metal-centered p, d or f orbitals), the total number of c must be equal to {2*mag.l+1}')
                print('o = other occupied -> Orbitals belonging to the active space but not to the model one (other metal-centered orbitals or ligand-centered orbitals)')
                print('v = vacant orbitals -> Orbitals with a strict occupancy of 0 ("virtual")')
                print('EXIT Program')                  
                exit()                                 
            if len(nb_type_mo) == 1:
                nb_type_mo.insert(0,'1')

            if len(nb_type_mo) != 1:
                del nb_type_mo[-1]      
                tmp = ''.join(nb_type_mo)
                try:nb_mo=int(tmp)   
                except ValueError:
                    print('Invalid cas_mo input')
                    print('There are at most 4 categories of MOs in the active space and the sum of (d+c+o+v) must match the number of active orbitals:')
                    print('Before the d, c, o, v definition you should give an interger')
                    print(f'(The total sum of c type must be equal to {2*mag.l+1})')
                    print(f'example 1: -cas_mo 2d {2*mag.l+1}c')
                    print(f'Here there 2 first CFS have an occupancy of 2 and the {2*mag.l+1} next CFS are on the magnetic center') 
                    print(f'example 2: -cas_mo 4d {2*mag.l+1}c 3v')
                    print(f'Here there 4 first CFS have an occupancy of 2, the {2*mag.l+1} next CFS are on the magnetic center and the 3 last CFS have an occupancy of 0') 
                    print('EXIT Program')                  
                    exit()                                 
        for mo in cas_mo:
            nb_type_mo = list(mo) 
            if nb_type_mo[-1] == 'c':
                del nb_type_mo[-1]      
                tmp = ''.join(nb_type_mo)
                mo_c+=int(nb_type_mo[0])
        if mo_c != 2*mag.l+1:
            print('Invalid cas_mo input')
            print(f'You must give {2*mag.l+1}=2l+1 type c orbitals')
            print(f'example 1: -cas_mo {2*mag.l+1}c')
            print(f'Here it is the defaut value of the cas_mo variable') 
            print(f'example 2: -cas_mo 3d {2*mag.l+1}c 1v')
            print(f'Here the total c type take the value of 2l+1={2*mag.l+1}') 
            print('EXIT Program')
            exit()        
    else:
        cas_mo = [str(2*mag.l+1)+'c']

    return mag, cas_mo, user_sr_root, user_so_root
                 





   


   



