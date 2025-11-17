import sympy as sp
from sympy.physics.quantum import Ket
from itertools import combinations


class BasisStates:
    def __init__(self, l=3, nel=0) -> None:
        self.l = l
        self.nel = nel
        self.nel_unpaired = (2*l+1) - abs(nel - (2*l+1))
        self.L = self._get_lvalues()
        self.S = ((2*self.l+1) - abs(self.nel-(2*self.l+1))) * sp.Rational(1,2)
        self.basis_lmlsms_uc = self._calculate_basis_lmlsms_uc()
        self.basis_lmlsms_c = self._calculate_basis_lmlsms_c()
        self.U = self._calculate_cg_matrix()
        self.casscf_coef = {}  # to store casscf coefficients from ourca outputs
    
    def _get_lvalues(self):
        if self.nel_unpaired in [1, 6]:
            return [3]
        elif self.nel_unpaired == 2:
            return [1, 3, 5]
        elif self.nel_unpaired in [3, 4]:
            return [0, 2, 3, 4, 6]
        elif self.nel_unpaired == 5:
            return [1, 3, 5]

    def _calculate_basis_lmlsms_uc(self):
        """Generate uncoupled basis states for up to 2 electrons."""
        ml_values = range(-self.l, self.l + 1)  # Magnetic quantum numbers
        ml_combs = combinations(ml_values, self.nel_unpaired)  # Electron distributions

        return [
            LMLSMS(L=self.l, Ml=ml_pair[0], S=1/2, Ms=1/2) if len(ml_pair) == 1
            else tuple(LMLSMS(L=self.l, Ml=ml, S=1/2, Ms=1/2) for ml in ml_pair)
            for ml_pair in ml_combs
        ]

    def _calculate_basis_lmlsms_c(self):
        """Generate coupled basis states."""
        states = []
        for l in self.L:  # Ensure l is defined before using it in range
            states.extend([
                LMLSMS(L=l, Ml=ml, S=self.S, Ms=self.S)
                for ml in range(-l, l + 1)
            ])
        return states
    

    def _calculate_cg_matrix(self):  # Max 2 electrons for now
        if self.nel_unpaired > 2:
            return "Sorry, no more than 2 electrons!"
        if self.nel_unpaired == 1:  # If 1 electron, return identity matrix
            return sp.eye(7)

        lml_bas = self.basis_lmlsms_uc  # Uncoupled basis
        LML_bas = self.basis_lmlsms_c  # Coupled basis
        dim = len(lml_bas)  # Basis dimension

        cg_mat = sp.zeros(dim, dim)  # Initialize zero matrix

        for idx, LML_ket in enumerate(LML_bas):
            lt, mlt = LML_ket.L, LML_ket.Ml  # Extract quantum numbers

            for jdx, lml_ket in enumerate(lml_bas):
                ml1, ml2 = lml_ket[0].Ml, lml_ket[1].Ml  # Extract uncoupled ML values

                # Compute Clebsch-Gordan coefficient
                cg_coeff = CG(self.l, ml1, self.l, ml2, lt, mlt).doit()
                if cg_coeff != 0:  # Only store non-zero values
                    cg_mat[idx, jdx] = cg_coeff

        return np.array(cg_mat).astype(np.float64)  # Convert to NumPy for efficiency
