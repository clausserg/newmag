import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from crystal_field_matrix import *
from clebsch_gordan import *

def extract_cf_prms(mag_center=None, heff=None, soc=False, Bkq=None):
    cfp_extracted = {}

    # Convert heff to numeric
    heff = np.array(heff).astype(np.complex128)

    if soc:
        basis = mag_center.basis_lmlsms
        bas_jmj = mag_center.basis_jmj
        ucg = cg_lml_jmj(basis, bas_jmj)
        ucg_conj_t = ucg.conj().T
        heff_eff = ucg_conj_t @ heff @ ucg  # rotated heff
    else:
        basis = mag_center.basis_lmls
        heff_eff = heff  # unrotated

    # Shared matrix_cf cache
    cf_cache = {}

    def process_key_val(key_val):
        key, val = key_val
        if key not in cf_cache:
            cf_cache[key] = np.array(matrix_cf(mag_center, basis, key)).astype(np.complex128)
        cf_mat = cf_cache[key]
        num = np.trace(heff_eff @ cf_mat)
        denom = np.trace(cf_mat @ cf_mat)
        cfp_value = np.real(num / denom)
        return val, cfp_value

    # Parallel execution
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_key_val, kv) for kv in Bkq.items()]
        for future in as_completed(futures):
            val, result = future.result()
            cfp_extracted[val] = result

    return cfp_extracted
