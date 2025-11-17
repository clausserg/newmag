import sympy as sp
from sympy.physics.quantum import Ket

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
