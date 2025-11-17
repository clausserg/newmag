from aniso_parameters import *
from class_basis_states import *
from class_jmj import *
from class_lmlsms import *
from class_magnetic_center import *
from clebsch_gordan import *
from constants import *
from crystal_field_matrix import *
from dictionaries import *
from extract_caswfn_orca import *
from extract_energies_orca import *
from extract_sociwfn_orca import *
from g_tensor import *
from helper_functions import *
from inner_product import *
from ito_extraction import *
from lowdin_descloizeaux import *
from magnetic_moments import *
from real_to_spherical_harmonics import *
from scalar_hamiltonian import *
from soc_hamiltonian import *
from soc_matrix import *
from soc_operator import *
from state_composition import *
from stevens_operators import *
from susceptibility import *
from symbols import *
import argparse
from tabulate import tabulate
import pandas as pd


parser = argparse.ArgumentParser(
prog='NewMag program',
description='Calculation of crystal field parameters from ORCA outputs.\n !!!The program works currently only for the f1 configuration!!!')

answer = input("Do you have an f1 configuration? (yes/no): ").strip().lower()

if answer in ["yes", "y", "oui", "da"]:
    print("Great! Continuing with the program...")
else:
    print("!!!The program works currently only for the f1 configuration!!!\nFeature not implemented. Exiting...")
    exit()


parser.add_argument(
"filename",
type=str,
help="[ORCA] output file\n")
args = parser.parse_args()

out = args.filename

mag = MagneticCenter(l=3, nel=1)
print("\nThis is:", mag, "\n")

print("\nCASSCF numerical Hamiltonian:")
H_sr =  build_heff_sr(mag_center=mag, output=out, level="casscf-sr")
H_sr_print = pretty_matrix(H_sr, mag.basis_lmls)
print(tabulate(H_sr_print, headers='keys', tablefmt='rounded_grid'))

print("\nCASSCF compostion of states:")
compo_sr_print=pretty_composition(H_sr, level="casscf-sr")
print(tabulate(compo_sr_print, headers='keys', tablefmt='rounded_grid'))

print('\n\n')

print("\nSO-CASSCF numerical Hamiltonian:")
H_soc = build_heff_so(mag_center=mag, output=out, level="casscf-so")
H_soc_print = pretty_matrix(H_soc, mag.basis_jmj)
print(tabulate(H_soc_print, headers='keys', tablefmt='rounded_grid'))

print("\nSO-CASSCF compostion of states:")
compo_soc_print=pretty_composition(H_soc, level="casscf-so")
print(tabulate(compo_soc_print, headers='keys', tablefmt='rounded_grid'))

print("\n\nCystal Field Parameters:")
srp = extract_cf_prms(mag_center=mag, heff=H_sr, soc=False, Bkq=Bkq)
sop = extract_cf_prms(mag_center=mag, heff=H_soc, soc=True, Bkq=Bkq)
parm_print=print_cf_params(srp, sop, out,mag)
print(tabulate(parm_print, headers='keys', tablefmt='presto'))

print("\n\n g-values:")
g_principal={}
for root in range(7):
    g_principal[f'KD {root}']=g_tensor(mag, H_soc, root)
g_principal=pd.DataFrame(g_principal)
g_principal=g_principal.map(lambda x: f"{x.real:6.3f}")
print(tabulate(g_principal, headers='keys', tablefmt='rounded_grid'))

print("\n\n Magnetic susceptibility plot (chi*T vs. T):")
susceptibility_powder(mag, H_soc, 300, output="cas", exp_file=None)
