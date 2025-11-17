import numpy as np
import re
from dictionaries import *

def aniso_params(file_path: str, mag):

    b2 = [B22m, B21m, B20, B21, B22]
    b4 = [B44m, B43m, B42m, B41m, B40, B41, B42, B43, B44]     
    b6 = [B66m, B65m, B64m, B63m, B62m, B61m, B60, B61, B62, B63, B64, B65, B66]

    dict_single_aniso_L={}
#    dict_single_aniso_L['single_ansio_L']={}
    for group in (b2, b4, b6):      
        for par in group:           
            dict_single_aniso_L[par]=[]

    a, b, g = mag.factor_abg()
    a=np.float64(a)
    b=np.float64(b)
    g=np.float64(g)

    header_line = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L =  3"
    table_header = "k |  q  |    (K)^2    |         B(k,q)        |"
    stop_marker = "********************************************************************************"

    Bkq = {}
    in_section = False
    in_table = False

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()

            # Step 1: wait until the main header is found
            if not in_section:
                if header_line in line:
                    in_section = True
                continue

            # Step 2: wait for the right table header
            if in_section and not in_table:
                if table_header in line:
                    in_table = True
                continue

            # Step 3: parse the table lines until stop marker
            if in_table:
                # Stop if we reach the stop marker
                if stop_marker in line:
                    break

                # Skip separators
                if line.startswith("----") or line.startswith("------------------------------------------------|"):
                    continue

                # Match table entries
                match = re.match(
                    r"^\s*(\d+)\s*\|\s*(-?\d+)\s*\|\s*[\d.E+-]+\s*\|\s*([-\d.E+]+)\s*\|",
                    line,
                )
                if match:
                    k = int(match.group(1))
                    q = int(match.group(2))
                    B_value = np.float64(match.group(3))
                    if k==2:
                        if q==-2:
                            dict_single_aniso_L[B22m]=B_value/a
                        elif q==-1:
                            dict_single_aniso_L[B21m]=B_value/a
                        elif q==0:
                            dict_single_aniso_L[B20]=B_value/a
                        elif q==1:
                            dict_single_aniso_L[B21]=B_value/a
                        elif q==2:
                            dict_single_aniso_L[B22]=B_value/a
                    elif k==4:
                        if q==-4:
                            dict_single_aniso_L[B44m]=B_value/b
                        elif q==-3:
                            dict_single_aniso_L[B43m]=B_value/b
                        elif q==-2:
                            dict_single_aniso_L[B42m]=B_value/b
                        elif q==-1:
                            dict_single_aniso_L[B41m]=B_value/b
                        elif q==0:
                            dict_single_aniso_L[B40]=B_value/b
                        elif q==1:
                            dict_single_aniso_L[B41]=B_value/b
                        elif q==2:
                            dict_single_aniso_L[B42]=B_value/b
                        elif q==3:
                            dict_single_aniso_L[B43]=B_value/b
                        elif q==4:
                            dict_single_aniso_L[B44]=B_value/b
                    elif k==6:
                        if q==-6:
                            dict_single_aniso_L[B66m]=B_value/g
                        elif q==-5:
                            dict_single_aniso_L[B65m]=B_value/g
                        elif q==-4:
                            dict_single_aniso_L[B64m]=B_value/g
                        elif q==-3:
                            dict_single_aniso_L[B63m]=B_value/g
                        elif q==-2:
                            dict_single_aniso_L[B62m]=B_value/g
                        elif q==-1:
                            dict_single_aniso_L[B61m]=B_value/g
                        elif q==0:
                            dict_single_aniso_L[B60]=B_value/g
                        elif q==1:
                            dict_single_aniso_L[B61]=B_value/g
                        elif q==2:
                            dict_single_aniso_L[B62]=B_value/g
                        elif q==3:
                            dict_single_aniso_L[B63]=B_value/g
                        elif q==4:
                            dict_single_aniso_L[B64]=B_value/g
                        elif q==5:
                            dict_single_aniso_L[B65]=B_value/g
                        elif q==6:
                            dict_single_aniso_L[B66]=B_value/g

    return dict_single_aniso_L

#def aniso_params(file_path: str):
#    """
#    Extracts the B(k,q) coefficients from the Crystal-Field Hamiltonian table
#    that follows the header line about L = 3.
#
#    Args:
#        file_path (str): Path to the output file.
#
#    Returns:
#        dict[(int, int), float]: Dictionary mapping (k, q) -> B(k,q).
#    """
#    header_line = "CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L =  3"
#    table_header = "k |  q  |    (K)^2    |         B(k,q)        |"
#    stop_marker = "********************************************************************************"
#
#    Bkq = {}
#    in_section = False
#    in_table = False
#
#    with open(file_path, "r", encoding="utf-8") as f:
#        for line in f:
#            line = line.rstrip()
#
#            # Step 1: wait until the main header is found
#            if not in_section:
#                if header_line in line:
#                    in_section = True
#                continue
#
#            # Step 2: wait for the right table header
#            if in_section and not in_table:
#                if table_header in line:
#                    in_table = True
#                continue
#
#            # Step 3: parse the table lines until stop marker
#            if in_table:
#                # Stop if we reach the stop marker
#                if stop_marker in line:
#                    break
#
#                # Skip separators
#                if line.startswith("----") or line.startswith("------------------------------------------------|"):
#                    continue
#
#                # Match table entries
#                match = re.match(
#                    r"^\s*(\d+)\s*\|\s*(-?\d+)\s*\|\s*[\d.E+-]+\s*\|\s*([-\d.E+]+)\s*\|",
#                    line,
#                )
#                if match:
#                    k = int(match.group(1))
#                    q = int(match.group(2))
#                    B_value = float(match.group(3))
#                    Bkq[(k, q)] = B_value
#
#    return Bkq
