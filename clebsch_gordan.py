from sympy.physics.quantum.cg import CG
import numpy as np


def cg_lml_jmj(srbas, sobas):
    dim = len(sobas)
    Ucg = np.zeros((dim, len(srbas)), dtype=np.float64 )

    for i, jket in enumerate(sobas):
        j, mj = float(jket.J), float(jket.Mj)
        for jdx, lket in enumerate(srbas):
            l, ml, s, ms = float(lket.L), float(lket.Ml), float(lket.S), float(lket.Ms)

            if not np.isclose(ml + ms, mj):
                continue

            coeff = CG(l, ml, s, ms, j, mj).doit()
            if coeff != 0:
                Ucg[i, jdx] = float(coeff)

    return Ucg
