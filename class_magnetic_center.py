import sympy as sp
import numpy as np
from sympy.physics.quantum.cg import CG
from class_lmlsms import LMLSMS
from class_jmj import JMJ
from stevens_operators import *


class MagneticCenter:
    def __init__(self, l=3, nel=0) -> None:
        self.nel = nel
        self.l = l
        self.S = ((2*self.l+1) - abs(self.nel-(2*self.l+1))) * sp.Rational(1,2)  
        self.L = self._calculate_lmax()
        self.basis_lmls, self.basis_lmlsms = self._construct_basis_lmlsms()
        self.basis_jmj, self.lande_factors = self._construct_basis_jmj()
        

    def __repr__(self):
        return f"Magnetic center with L={self.L} and S={self.S}"
    
    def _calculate_lmax(self):
        l_max = 0
        for idx, jdx in enumerate(range(self.l, -self.l, -1)):
            if idx+1 <= abs(self.nel - (2*self.l+1)):
                l_max += jdx
        return l_max
        
    
    def _construct_basis_lmlsms(self):
        """
        Returns a list of LMLSMS objects sorted by increasing Ml and Ms.
        """
        ms = [-(self.S - idx) for idx in range(int(2*self.S + 1))]
        ml = [-(self.L - idx) for idx in range(int(2*self.L + 1))]
        
        lmls   = [LMLSMS(L=self.L, Ml=jdx, S=self.S, Ms=self.S) for jdx in ml]
        lmlsms = [LMLSMS(L=self.L, Ml=jdx, S=self.S, Ms=idx) for idx in ms for jdx in ml]
        return lmls, lmlsms

    def _construct_basis_jmj(self):
        """
        Fills the basis_jmj attribute with JMJ objects and populates
        their lmlsms_expansion attributes using Clebsch-Gordan algebra.
        """
        # Possible J values: |L - S| to L + S
        j_values = [sp.Rational(j) / 2 for j in range(int(abs(self.L - self.S) * 2), int((self.L + self.S) * 2) + 1, 2)]

        basis_jmj = []
        lande = {}
        for idx in j_values:
            lande[idx] = self.lande_g(idx, self.S, self.L)

        for j in j_values:
            # Generate all possible Mj values for this J
            m_j_values = [-(j - idx) for idx in range(int(2 * j + 1))]
            for m_j in m_j_values:
                jmj_obj = JMJ(J=j, Mj=m_j)
                # Populate lmlsms_expansion with LMLSMS objects weighted by Clebsch-Gordan coefficients
                for lmlsms in self.basis_lmlsms:
                    if lmlsms.Ml + lmlsms.Ms == m_j:
                        # Calculate the Clebsch-Gordan coefficient
                        cg_coef = CG(lmlsms.L, lmlsms.Ml, lmlsms.S, lmlsms.Ms, j, m_j).doit()
                        if cg_coef != 0:  # Only include non-zero terms
                            jmj_obj.add_lmlsms(LMLSMS(L=lmlsms.L, Ml=lmlsms.Ml, S=lmlsms.S, Ms=lmlsms.Ms, coef=cg_coef))
                # attach the jmj object to the basis_jmj
                basis_jmj.append(jmj_obj)

        return basis_jmj, lande
    
    def lande_g(self, J, S, L):
        # Compute the Lande g-factor.
        numerator = J*(J+1) + S*(S+1) - L*(L+1)
        denominator = 2 * J * (J+1)
        return np.float64(1 + numerator / denominator)
    
    def factor_abg(self, myket=None):
        # test if myket is SF or SO basis function
        LaL, LbL, LgL, lgl = 0, 0, 0, sp.Rational(4,11*13*27)  # apply sign afterwards
        JaJ, JbJ, JgJ = 0, 0, 0

        LaL = (2*(2*self.l+1 - 4*self.S)) / ((2*self.l-1)*(2*self.l+3)*(2*self.L-1))
        LbL = LaL * ((3*(3*(self.l-1)*(self.l+2)-7*(self.l-2*self.S)*(self.l+1-2*self.S))) / (2*(2*self.l-3)*(2*self.l+5)*(self.L-1)*(2*self.L-3)))

        unpaired_el = list(range(0, self.S * 2, 1))  # total nr of unpaired eectrons
        ml_values = list(range(self.l, -self.l-1, -1))  # the ml value of the unpaired electrons
        el_ml = dict(zip(unpaired_el, ml_values))  # zip together the unpaired electrons with their ml value

        for idx, jdx in el_ml.items():
            LgL += sto_O60(JMJ(J=self.l, Mj=jdx, coef=1))[0].coef
        LgL *= lgl
        LgL /= sto_O60(JMJ(J=self.L, Mj=self.L, coef=1))[0].coef

        if self.nel < (2*self.l+1):  # 2l+1 gives 7 for f, 5 for d, etc. (change sign if smaller than half filled shell)
            LaL *= -1
            LbL *= -1
            LgL *= -1
        
        if myket:  # if yket is given, we deal with a JMJ manifold
            for term in myket.lmlsms_expansion:
                JaJ += term.coef**2 * sto_O20(JMJ(J=term.L, Mj=term.Ml, coef=1))[0].coef
                JbJ += term.coef**2 * sto_O40(JMJ(J=term.L, Mj=term.Ml, coef=1))[0].coef
                JgJ += term.coef**2 * sto_O60(JMJ(J=term.L, Mj=term.Ml, coef=1))[0].coef
            JaJ = (LaL * JaJ) / sto_O20(myket)[0].coef
            JbJ = (LbL * JbJ) / sto_O40(myket)[0].coef
            JgJ = (LgL * JgJ) / sto_O60(myket)[0].coef
            return tuple(0 if factor==sp.nan else factor for factor in (JaJ, JbJ, JgJ))
        return tuple(0 if factor==sp.nan else factor for factor in (LaL, LbL, LgL))
