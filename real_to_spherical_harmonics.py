import numpy as np


def harmonics_re_to_sph(basis):
    # coded according to page 4 in https://doi.org/10.1016/S0166-1280(97)00185-1
    
    dim = len(basis)
    U = np.zeros((dim, dim), dtype=np.complex128)

    for i, m1 in enumerate(basis):
        for j, m2 in enumerate(basis):
            ml1, ml2 = m1.Ml, m2.Ml
            if abs(ml1) != abs(ml2):
                continue
            elif ml1 == 0 and ml2 == 0:
                U[i, j] = 1
            elif ml1 == ml2 and ml1 < 0:
                U[i, j] = 1j / np.sqrt(2)
            elif ml1 == ml2 and ml1 > 0:
                U[i, j] = ((-1) ** ml1) / np.sqrt(2)
            elif abs(ml1) == abs(ml2) and ml1 < 0 and ml2 > 0:
                U[i, j] = -1j * ((-1) ** ml2) / np.sqrt(2)
            elif abs(ml1) == abs(ml2) and ml1 > 0 and ml2 < 0:
                U[i, j] = 1 / np.sqrt(2)

    return U
