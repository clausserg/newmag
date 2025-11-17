from stevens_operators import *
from symbols import *

# Gather tuples for the rank and projection of the operator
Okq = {
    (2, 0): sto_O20, (2, 1): sto_O21, (2, 2): sto_O22, (2,-1): sto_O21m, (2,-2): sto_O22m,
    (4, 0): sto_O40, (4, 1): sto_O41, (4, 2): sto_O42, (4, 3): sto_O43, (4, 4): sto_O44, (4,-1): sto_O41m, (4,-2): sto_O42m, (4,-3): sto_O43m, (4,-4): sto_O44m,
    (6, 0): sto_O60, (6, 1): sto_O61, (6, 2): sto_O62, (6, 3): sto_O63, (6, 4): sto_O64, (6, 5): sto_O65, (6, 6): sto_O66, (6,-1): sto_O61m, (6,-2): sto_O62m, (6,-3): sto_O63m, (6,-4): sto_O64m, (6,-5): sto_O65m, (6,-6): sto_O66m
}

# Gather the CFPs in a dictionary
Bkq = {
    (2, 0): B20, (2, 1): B21, (2, 2): B22, (2,-1): B21m, (2,-2): B22m,
    (4, 0): B40, (4, 1): B41, (4, 2): B42, (4, 3): B43, (4, 4): B44, (4,-1): B41m, (4,-2): B42m, (4,-3): B43m, (4,-4): B44m,
    (6, 0): B60, (6, 1): B61, (6, 2): B62, (6, 3): B63, (6, 4): B64, (6, 5): B65, (6, 6): B66, (6,-1): B61m, (6,-2): B62m, (6,-3): B63m, (6,-4): B64m, (6,-5): B65m, (6,-6): B66m
}

# Make a dictionary relating the STOs with the proportionality constants used in the ITOs of SINGLE_ANISO
ito_cst = {
    (2, 0): 1, (2, 1): sp.sqrt(6), (2, 2): sp.sqrt(3)/sp.sqrt(2), (2,-1): sp.sqrt(6), (2,-2): sp.sqrt(3)/sp.sqrt(2),
    (4, 0): 1, (4, 1): 2*sp.sqrt(5), (4, 2): sp.sqrt(10), (4, 3): 2*sp.sqrt(35), (4, 4): sp.sqrt(35)/sp.sqrt(2), (4,-1): 2*sp.sqrt(5), (4,-2): sp.sqrt(10), (4,-3): 2*sp.sqrt(35), (4,-4): sp.sqrt(35)/sp.sqrt(2),
    (6, 0): 1, (6, 1): sp.sqrt(42), (6, 2): sp.sqrt(105)/2, (6, 3): sp.sqrt(105), (6, 4): 3*sp.sqrt(7)/sp.sqrt(2), (6, 5): 3*sp.sqrt(77), (6, 6): sp.sqrt(231)/2, (6,-1): sp.sqrt(42), (6,-2): sp.sqrt(105)/2, (6,-3): sp.sqrt(105), (6,-4): 3*sp.sqrt(7)/sp.sqrt(2), (6,-5): 3*sp.sqrt(77), (6,-6): sp.sqrt(231)/2
}
