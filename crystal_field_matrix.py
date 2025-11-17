import sympy as sp
from dictionaries import *
from class_jmj import *
from class_lmlsms import *

def matrix_cf(magcent, basis, sto):
    Okq_mat = sp.zeros(len(basis), len(basis))
    for idx in range(len(basis)):
        for jdx in range(len(basis)):
            aket, bket = basis[idx], basis[jdx]
            expansion = Okq[sto](bket)
            for term in expansion:
                if sto[0] == 2 and isinstance(bket, JMJ):
                    term.coef *= magcent.factor_abg(bket)[0]
                if sto[0] == 2 and isinstance(bket, LMLSMS):
                    term.coef *= magcent.factor_abg()[0]
                if sto[0] == 4 and isinstance(bket, JMJ):
                    term.coef *= magcent.factor_abg(bket)[1]
                if sto[0] == 4 and isinstance(bket, LMLSMS):
                    term.coef *= magcent.factor_abg()[1]
                if sto[0] == 6 and isinstance(bket, JMJ):
                    term.coef *= magcent.factor_abg(bket)[2]
                if sto[0] == 6 and isinstance(bket, LMLSMS):
                    term.coef *= magcent.factor_abg()[2]
                if term.Ket == aket.Ket:
                    Okq_mat[idx, jdx] += term.coef
    return Okq_mat
