from collections import defaultdict


def soci_wfn_orca(output, level="casscf-so"):

    level = level.lower()
    if level not in {"casscf-so", "nevpt2-so"}:
        raise ValueError("level should be 'casscf-so' or 'nevpt2-so'")
        
    with open(output, "r") as file:
        content = file.readlines()

    wavefunctions = {}
    state_index = None
    current_wfn = defaultdict(complex)

    # Map level to keyword in trigger line
    level_map = {
        "casscf-so": "CASSCF",
        "nevpt2-so": "NEVPT2",
    }

    keyword = level_map.get(level.lower(), "CASSCF")  # default fallback
    trigger = f"QDPT WITH {keyword} DIAGONAL ENERGIES"
    
    parsing = False
    
    for line in content:
        # Wait until trigger line appears
        if not parsing:
            if trigger in line:
                parsing = True
            continue  # skip lines until trigger is found
        
        if " STATE" in line and ":" in line:
            # When a new state starts, save the previous one
            if state_index is not None and current_wfn:
                wavefunctions[state_index] = dict(current_wfn)
                current_wfn.clear()

            # Extract state index
            try:
                state_index = int(line.split()[1].strip(":"))
            except (IndexError, ValueError):
                state_index = None

        elif state_index is not None:
            row = line.split()
            if len(row) == 8 and row[4] == "0":
                key = (row[5], row[6], row[7])
                coef = float(row[1]) + float(row[2]) * 1j
                current_wfn[key] += coef

        # Stop when passing known end blocks
        if "Center of nuclear charge           = (" in line or "COMPUTING QDPT PROPERTIES" in line:
            if state_index is not None and current_wfn:
                wavefunctions[state_index] = dict(current_wfn)
            break

    return wavefunctions
