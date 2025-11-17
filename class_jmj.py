import sympy as sp
from sympy.physics.quantum import Ket

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
