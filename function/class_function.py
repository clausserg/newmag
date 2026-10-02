import sympy as sp
import numpy as np
from itertools import product
from sympy.physics.quantum.cg import CG
from itertools import combinations
from sympy.physics.quantum import Ket
from copy import deepcopy
#from function.crystal_field import sto_O20, sto_O40, sto_O60

def sto_O20(myket):            
    mj = myket.Mj              
    x = myket.J * (myket.J + 1)     
    newket = deepcopy(myket)   
    newket.coef = (3 * mj*mj - x)
    return [newket]            
def sto_O40(myket):                       
    mj = myket.Mj                         
    x = myket.J * (myket.J + 1)           
    new_coef = 35 * mj**4 - (30*x -25)*mj*mj + 3*x*x - 6*x
    newket = deepcopy(myket)              
    newket.coef = new_coef                
    return [newket]                       
def sto_O60(myket):                        
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    new_coef = 231*mj**6 - (315*x-735)*mj**4 + (105*x*x-525*x+294)*mj*mj - 5*x**3 + 40*x*x - 60*x
    newket = deepcopy(myket)               
    newket.coef = new_coef                 
    return [newket]                        


class MagneticCenter:
    def __init__(self, l=3, nel=0) -> None:
        self.nel = nel
        self.l = l
        self.S = ((2*self.l+1) - abs(self.nel-(2*self.l+1))) * sp.Rational(1,2)  
        self.L = self._calculate_lmax()
        self.basis_lmls, self.basis_lmlsms = self._construct_basis_lmlsms()
        self.basis_jmj, self.lande_factors = self._construct_basis_jmj()
        self.basis_ne_ml_uncoupled, self.basis_ne_lmms_uncoupled = self._construct_basis_lmlsms_uncoupled()

    def __repr__(self):
        return f"Magnetic center with L={self.L} and S={self.S}"
    
    def _calculate_lmax(self):
        l_max = 0
        for idx, jdx in enumerate(range(self.l, -self.l, -1)):
            if idx+1 <= abs(self.nel - (2*self.l+1)):
                l_max += jdx
        return l_max
        
    
    def _construct_basis_lmlsms_uncoupled(self):
        """
        Returns a tuple of tuple of lmssms uncoupled
        """
        ml_list=[-(self.l - idx) for idx in range(int(2*self.l + 1))]
        Ms_list = [str(-(self.S - idx)) for idx in range(int(2*self.S + 1))]
        lml_uncoupled=[]
        if self.nel < (2*self.l+1):
            for nel in range(self.nel):
                ne_ml_uncoupled=tuple(product(ml_list, repeat=self.nel))
        elif self.nel > (2*self.l+1):
            for nel in range(self.nel-2*(self.nel-int(2*self.l+1))):
                ne_ml_uncoupled=tuple(product(ml_list, repeat=self.nel-2*(self.nel-int(2*self.l+1))))
        ne_mlms_uncoupled=tuple(product(ne_ml_uncoupled,Ms_list,repeat=1))
        ne_mlms_uncoupled=sorted(ne_mlms_uncoupled, key = lambda x: str(x[1]))
        return ne_ml_uncoupled, ne_mlms_uncoupled

        
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
        if self.nel > (2*self.l+1):
            j_values = list(reversed(j_values))

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

    def Stevens_coeff(self):
        if self.l == 3:
            if self.nel == 1:
                return -17.50000, 157.50000, 0 
            elif self.nel == 2:
                return -47.59615, -1361.25000, 16395.05515
            elif self.nel == 3:
                return -155.57143, -3435.15441, -26324.13065
            elif self.nel == 4:
                return 129.64286, 2453.68172, 16452.58166
            elif self.nel == 5:
                return 24.23077, 399.80769, 0 
            elif self.nel == 6:
                return 0, 0, 0 
            elif self.nel == 8:
                return -99.00000, 8167.50000, -891891.00000
            elif self.nel == 9:
                return -157.50000, -16891.87500, 966215.25000
            elif self.nel == 10:
                return -450.00000, -30030.00000, -772972.20000
            elif self.nel == 11:
                return 393.75000, 22522.50000, 483107.62500
            elif self.nel == 12:
                return 99.00000, 6125.62500, -178378.20000
            elif self.nel == 13:
                return 31.50000, -577.50000, 6756.75000
        else: return 1,1,1

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

class LMLSMS:
    def __init__(self, L=0, Ml=0, S=sp.Rational(0), Ms=0, coef=1) -> None:
        self.L = sp.Rational(L)
        self.Ml = sp.Rational(Ml)
        # let's add J and Mj referencing L and Ml, to help with Stevens Ops.
        self.J = sp.Rational(L)
        self.Mj = sp.Rational(Ml)
        # done
        self.S   = sp.Rational(S)
        self.Ms = sp.Rational(Ms)
        self.coef = coef
        self.Ket = Ket(sp.Rational(L), sp.Rational(Ml), sp.Rational(S), sp.Rational(Ms))
        self.coef_casscf = None

    def __repr__(self):
        myself = "{}*|L={}, Ml={}, S={}, Ms={}>".format(self.coef, self.L, self.Ml, self.S, self.Ms)
        return "LMLSMS basis state:" + myself

    def __mul__(self, other):
        self.coef *= other
        return self

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if other != 0:
            self.coef /= other
        else:
            raise ValueError("Division by zero")
        return self

    # Spin ladder operators, S+ and S-
    def op_Sp(self, center=None):
        new_coef = self.coef * sp.sqrt(self.S * (self.S + 1) - self.Ms * (self.Ms + 1))
        newKet = LMLSMS(self.L, self.Ml, self.S, self.Ms+1, coef=new_coef)
        return newKet

    def op_Sm(self, center=None):
        new_coef = self.coef * sp.sqrt(self.S * (self.S + 1) - self.Ms * (self.Ms - 1))
        newKet = LMLSMS(self.L, self.Ml, self.S, self.Ms-1, coef=new_coef)
        return newKet

    # Operator Sz
    def op_Sz(self, center=None):
        new_coef = self.coef * self.Ms
        newKet = LMLSMS(self.L, self.Ml, self.S, self.Ms, coef=new_coef)
        return newKet

    # Operators Sx and Sy
    def op_Sx(self, center=None):
        splus = self.op_Sp()
        sminus = self.op_Sm()
        splus.coef *= sp.Rational(1,2)
        sminus.coef *= sp.Rational(1,2)
        return [idx for idx in [splus, sminus] if idx.coef !=0]

    def op_Sy(self, center=None):
        splus = self.op_Sp()
        sminus = self.op_Sm()
        splus.coef *= (+1)/(2*sp.I)
        sminus.coef *= (-1)/(2*sp.I)
        return [idx for idx in [splus, sminus] if idx.coef !=0]

    # Lz operator
    def op_Jz(self):
        new_coef = self.coef * self.Ml
        newKet = LMLSMS(self.L, self.Ml, self.S, self.Ms, coef=new_coef)
        return newKet

    # L+ operator
    def op_Jp(self):
        new_coef = self.coef * sp.sqrt(self.L * (self.L + 1) - self.Ml * (self.Ml + 1))
        newKet = LMLSMS(self.L, self.Ml+1, self.S, self.Ms, coef=new_coef)
        return newKet

    # L- operator
    def op_Jm(self):
        new_coef = self.coef * sp.sqrt(self.L * (self.L + 1) - self.Ml * (self.Ml - 1))
        newKet = LMLSMS(self.L, self.Ml-1, self.S, self.Ms, coef=new_coef)
        return newKet

    # Operators Lx and Ly
    def op_Jx(self, center=None):
        jplus = self.op_Jp()
        jminus = self.op_Jm()
        jplus.coef *= sp.Rational(1,2)
        jminus.coef *= sp.Rational(1,2)
        return [idx for idx in [jplus, jminus] if idx.coef !=0]

    def op_Jy(self, center=None):
        jplus = self.op_Jp()
        jminus = self.op_Jm()
        jplus.coef *= (+1)/(2*sp.I)
        jminus.coef *= (-1)/(2*sp.I)
        return [idx for idx in [jplus, jminus] if idx.coef !=0]

    # overlap with another ket of same class
    def ovl(self, other):
        tot_ovl = 0
        for aket in other:
            if aket.Ket == self.Ket:
                tot_ovl += (aket.coef * self.coef)
        return tot_ovl

    def times_cst(self, cst):
        new_coef = self.coef * cst
        newKet = LMLSMS(self.L, self.Ml, self.S, self.Ms, coef=new_coef)
        return newKet

class JMJ:
    def __init__(self, J, Mj, coef=1):
        self.J  = sp.Rational(J)
        self.Mj = sp.Rational(Mj)
        self.coef = coef
        self.Ket = Ket(sp.Rational(J), sp.Rational(Mj))
        self.lmlsms_expansion = []

    def __repr__(self):
        myself = "{}*|J={}, Mj={}>".format(self.coef, self.J, self.Mj)
        return "JMj basis state:" + myself

    # operator Jz
    def op_Jz(self):
        new_coef = self.coef * self.Mj
        newKet = JMJ(self.J, self.Mj, coef=new_coef)
        return newKet

    # ladder operators J+ and J-
    def op_Jp(self):
        new_coef = self.coef * sp.sqrt(self.J * (self.J + 1) - self.Mj * (self.Mj + 1))
        newKet = JMJ(self.J, self.Mj+1, coef=new_coef)
        return newKet

    def op_Jm(self):
        new_coef = self.coef * sp.sqrt(self.J * (self.J + 1) - self.Mj * (self.Mj - 1))
        newKet = JMJ(self.J, self.Mj-1, coef=new_coef)
        return newKet

    # Operators Jx and Jy
    def op_Jx(self, center=None):
        jplus = self.op_Jp()
        jminus = self.op_Jm()
        jplus.coef *= sp.Rational(1,2)
        jminus.coef *= sp.Rational(1,2)
        return [idx for idx in [jplus, jminus] if idx.coef !=0]

    def op_Jy(self, center=None):
        jplus = self.op_Jp()
        jminus = self.op_Jm()
        jplus.coef *= (+1)/(2*sp.I)
        jminus.coef *= (-1)/(2*sp.I)
        return [idx for idx in [jplus, jminus] if idx.coef !=0]

    def add_lmlsms(self, lmlsms):
        self.lmlsms_expansion.append(lmlsms)

    def void_lmlsms_expansion(self):
        self.lmlsms_expansion = []

    def times_cst(self, cst):
        new_coef = self.coef * cst
        newKet = JMJ(self.J, self.Mj, coef=new_coef)
        return newKet
