import re


def get_energies_orca(orca_out, level="casscf-sr"):
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
            if line.startswith("CAS-SCF STATES FOR"):
                inside_block = True
                continue
            if inside_block:
                if line.strip().startswith("Spin") or "DENSITY MATRIX" in line:
                    break
                if line.strip().startswith("ROOT"):
                    row = line.split()
                    energies_hartree.append(float(row[3]))
                    #row = line.split()
                    #if len(row) >= 8:
                    #    energies.append(float(row[7]))
                    #else:
                    #    energies.append(0.0)

        # convertion coming form : https://physics.nist.gov/cgi-bin/cuu/Value?hrminv|search_for=hartree
        for e in range(len(energies_hartree)):
            energies.append((energies_hartree[e]-energies_hartree[0])*219474.63136314)

#    elif level == "casscf-so":
#        inside_block = False
#        for line in lines:
#            if "QDPT WITH CASSCF DIAGONAL ENERGIES" in line:
#                inside_block = True
#                continue
#            if inside_block:
#                if "COMPUTING QDPT PROPERTIES" in line :
#                    break
#                if re.match(r"^\s+STATE\s+\d+:\s+\d+", line):
#                    row = line.split()
#                    energies.append(float(row[2]))
#
    elif level == "casscf-so":
        inside_block = False
        for line in lines:
            if "QDPT WITH CASSCF DIAGONAL ENERGIES" in line:
                inside_block = True
                continue
            if inside_block:
                if "The threshold for printing" in line or "Eigenvectors" in line:
                    break
                if re.match(r"^\s*\d+", line):
                    row = line.split()
                    if len(row) > 1:
                        energies.append(float(row[1]))

    if level == "nevpt2-sr":
        energies.append(0.0)
        inside_block = False
        for line in lines:
            if line.strip().startswith("NEVPT2 TRANSITION ENERGIES"):
                inside_block = True
                continue
            if inside_block:
                if line.strip().startswith("NEVPT2 CORRECTION TO THE TRANSITION ENERGY"):
                    break
                this_line = line.strip().split()
                if len(this_line) >= 3 and this_line[0][0] in tuple("0123456789"):
                    row = line.split()
                    if len(row) >= 6:
                        energies.append(float(row[5]))
                    else:
                        energies.append(0.0)
                        
    elif level == "nevpt2-so":
        inside_block = False
        for line in lines:
            if "QDPT WITH NEVPT2 DIAGONAL ENERGIES" in line:
                inside_block = True
                continue
            if inside_block:
                if "The threshold for printing" in line or "Eigenvectors" in line:
                    break
                if re.match(r"^\s*\d+", line):
                    row = line.split()
                    if len(row) > 1:
                        energies.append(float(row[1]))
    return energies
