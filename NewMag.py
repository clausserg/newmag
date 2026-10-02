from function.build_hamiltonian import * 

from function.g_tensor import * 
from function.crystal_field import * 
from function.sym_point_group import *
from function.reading_outpute import *
from function.input_check import *
import argparse

parser = argparse.ArgumentParser(prog='NewMag program', description='Calculation of the crystal field parameters form ORCA or OpenMolcas output')
parser.add_argument("filename", type=str, help="Required: ORCA or OpenMolcas output file. Ideally: OpenMolcas_v25.06 and ORCA_v6.1. For OpenMolcas the name of this file must be: $Project.out (without any other dot in the name)\n")
parser.add_argument("-cfg", "--configuration", type=str, required=True, help="Required: Orbital configuration (e.g. f1, f3, d1, d4 ...)\n")
parser.add_argument("-ion", "--ion_filename", type=str, default=None, help="Optional: ORCA or OpenMolcas ion minimal CAS output file. Ideally: OpenMolcas_v25.06 and ORCA_v6.1. For OpenMolcas the name of this file must be: $Project.out\n")
parser.add_argument("-cas_mo", nargs='+', type=str, default=None, help="Optional: Order and type MO in the CFS (e.g. 2d 7c). Defaut value: {2l+1}c.\n")
parser.add_argument("-sr", "--selected_sr_root",nargs='+', type=int, default=None, help="Optional: List of selected scalar roots (the first root is 0, the second is 1, etc...)\n")
parser.add_argument("-so", "--selected_so_root",nargs='+', type=int, default=None, help="Optional: List of selected spin-orbit roots (the first root is 0, the second is 1, etc...)\n")
parser.add_argument("-ml", "--ml_order", nargs='+', type=int, default=None, help="Optional: order of ml that make up CFS. Defaut value for ORCA are the order of the AILFT. For OpenMolcas the order is reading in the filename\n")
parser.add_argument("-d", "--decimal", type=int, default=1, help="Optional: Specify the number of decimal digits to print. Default value is 1\n")
parser.add_argument("-soc_xyz", "--spin_orbit_anisotropy", action='store_true', help="Optional: calculate the spin-orbit anisotropy constants in terms of the x, y, z coordinates\n")
parser.add_argument("-p", "--print_large", action='store_true', help="Optional: Enable the large print option. Defaut value is False\n")
args = parser.parse_args()

out = args.filename
ion_out = args.ion_filename
cfg = args.configuration
decimal = args.decimal
print_large=args.print_large
user_sr_root = args.selected_sr_root
user_so_root = args.selected_so_root
user_ml_order = args.ml_order
cas_mo = args.cas_mo
soc_xyz = args.spin_orbit_anisotropy

# Checking input parameters for internal consistency 
mag, cas_mo, user_sr_root, user_so_root = input_check(out, ion_out, cfg, decimal, user_sr_root, user_so_root, user_ml_order, cas_mo)

#Symmetry and strcuture analyse
#The algorithm in this block is poorly designed; it can produce incorrect answers and simply crash.
#As a result, it is not a critical part of the program. 
print(40*'%')
print('Structural analysis')
print(40*'%')
try:symmetry_part(out)
except: pass

#Reading the outpute file 
#recovery the wavefunction/energie and the level of calculation
# if the file is an orca out and if there are an aiLFT block, the Heff contains the aiLFT hamiltonienne
print(40*'%')
print('Reading wavefunctions and energies')
print(40*'%')
print(mag, "\n")
wavefunction, energies, list_calc_level, Heff = reading_outpute(out, mag, ion_out, cas_mo, user_sr_root, user_so_root, user_ml_order)

#Building the effectif hamiltonienne form the ab initio outpute
#for each level of calculation
print(40*'%')
print('CRYSTAL FIELD EXTRACTION')
print(40*'%')
for level in list_calc_level:
    Heff[f'{level}-sr'], Heff[f'{level}-so'] = build_heff(wavefunction[f'{level}-sr'],wavefunction[f'{level}-so'],
                                                        energies[f'{level}-sr'],energies[f'{level}-so'],
                                                        mag, out, level, decimal, print_large)

#Extraction of CFPs (Bkq) using the ITO method
#the printing is in the function 
#the single_aniso reading is also in  the function
print(40*'%')
print('Extraction of the Crystal Field Parameters:')
print(40*'%')
print("(Please be patient, this could take a little time...)\n\n")
cfp = cfps_calc(Heff, mag, out, decimal)

# spin-orbit analyse to compute the SOC constante by the ITO method
# using the spherical approx. 
# Finally the eveluation of the model/ab-initio is done
print('\n\n')
print(40*'%')
print('Model spectra')
print(40*'%')
soc_analysies(Heff, cfp, mag, decimal, soc_xyz, print_large)        

print('\n\nNormal ending')
