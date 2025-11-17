import numpy as np
from collections import defaultdict
import sympy as sp

def composition(H_cf, state_index, level):

    eigvals, eigvecs = np.linalg.eigh(H_cf)
    vec = eigvecs[:, state_index]     
    probs = np.abs(vec)**2
    Energy=eigvals[state_index]

    composition = defaultdict(float)
    if level=='casscf-sr':
         basis = [
         ("L=3", -3), ("L=3", -2), ("L=3", -1),
         ("L=3",  0), 
         ("L=3", +1), ("L=3", +2), ("L=3", +3)
         ]         
         for (L, M), p in zip(basis, probs):
            key = f"|{L}, Mj=±{abs(M)}>" if M != 0 else f"|{L}, Mj=0>"
            composition[key] += 100 * p

    elif level=='casscf-so':
        basis = [
        ("J=5/2", sp.Rational(-5,2)), ("J=5/2", sp.Rational(-3,2)), ("J=5/2", sp.Rational(-1,2)),
        ("J=5/2", sp.Rational(+1,2)), ("J=5/2", sp.Rational(+3,2)), ("J=5/2", sp.Rational(+5,2)),
        ("J=7/2", sp.Rational(-7,2)), ("J=7/2", sp.Rational(-5,2)), ("J=7/2", sp.Rational(-3,2)), ("J=7/2", sp.Rational(-1,2)),
        ("J=7/2", sp.Rational(+1,2)), ("J=7/2", sp.Rational(+3,2)), ("J=7/2", sp.Rational(+5,2)), ("J=7/2", sp.Rational(+7,2))
        ]   
        for (J, M), p in zip(basis, probs):
            key = f"|{J}, Mj=±{abs(M)}>" if M != 0 else f"|{J}, Mj=0>"
            composition[key] += 100 * p

    composition=dict(composition)
    return Energy,composition
