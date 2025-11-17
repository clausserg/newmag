def casscf_wfn_orca(orca_output=None):
    ml_order = [0, 1, -1, 2, -2, 3, -3]
    all_wfns = {}  # {state_index: {config: coeff}}

    with open(orca_output, 'r') as file:
        content = file.readlines()

    # Locate the wavefunction printout section
    start, end = 0, 0
    for idx, line in enumerate(content):
        if "Spin-Determinant CI Printing" in line:
            start = idx
        elif any(x in line for x in [
            "CAS-SCF STATES FOR BLOCK", 
            "SA-CASSCF TRANSITION ENERGIES", 
            "DENSITY MATRIX"
        ]) and start != 0:
            end = idx
            break

    if start == 0 or end == 0:
        raise RuntimeError("Could not find Spin-Determinant CI section in the output.")

    wfn_content = content[start:end]
    current_state = None
    printed_wfn = []

    for line in wfn_content:
        line = line.strip()
        if line.startswith("ROOT") and ":" in line:
            if current_state is not None and printed_wfn:
                # save previous state
                state_dict = {}
                for cfg, coef in printed_wfn:
                    indexes = [key for key, val in cfg.items() if val == "u"]
                    if len(indexes) == 1:
                        state_dict[indexes[0]] = float(coef)
                    elif len(indexes) == 2:
                        state_dict[(indexes[0], indexes[1])] = float(coef)
                    else:
                        raise NotImplementedError("More than 2 orbitals not supported.")
                all_wfns[current_state] = state_dict
                printed_wfn = []

            current_state = int(line.split()[1].replace(":", ""))
            
        elif line.startswith("[") and "]" in line:
            row = line.split()
            if len(row) >= 2 and len(list(row[0].strip("[]"))) > 7 and row[0].strip("[]")[:2] == "22":
                key_str = row[0].strip("[]")[2:]
                coef = row[1]
                det = dict(zip(ml_order, key_str))
                printed_wfn.append((det, coef))   
            elif len(row) >= 2 and len(list(row[0].strip("[]"))) == 7:
                key_str = row[0].strip("[]")
                coef = row[1]
                det = dict(zip(ml_order, key_str))
                printed_wfn.append((det, coef))  

    # Save the last state
    if current_state is not None and printed_wfn:
        state_dict = {}
        for cfg, coef in printed_wfn:
            indexes = [key for key, val in cfg.items() if val == "u"]
            if len(indexes) == 1:
                state_dict[indexes[0]] = float(coef)
            elif len(indexes) == 2:
                state_dict[(indexes[0], indexes[1])] = float(coef)
            else:
                raise NotImplementedError("More than 2 orbitals not supported.")
        all_wfns[current_state] = state_dict
    return all_wfns
