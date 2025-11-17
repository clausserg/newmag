import sympy as sp
from copy import deepcopy


def sto_O20(myket):
    mj = myket.Mj
    x = myket.J * (myket.J + 1)
    newket = deepcopy(myket)
    newket.coef = (3 * mj*mj - x)
    return [newket]

def sto_O21(myket):
    # 1/4 [JzJp + JzJm + JpJz + JmJz]
    a = myket.op_Jp().op_Jz().times_cst(sp.Rational(1,4))
    b = myket.op_Jm().op_Jz().times_cst(sp.Rational(1,4))
    c = myket.op_Jz().op_Jp().times_cst(sp.Rational(1,4))
    d = myket.op_Jz().op_Jm().times_cst(sp.Rational(1,4))
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O21m(myket):
    # -i/4 [JzJp - JzJm + JpJz - JmJz]
    a = myket.op_Jp().op_Jz().times_cst(-sp.I/4)
    b = myket.op_Jm().op_Jz().times_cst(sp.I/4)
    c = myket.op_Jz().op_Jp().times_cst(-sp.I/4)
    d = myket.op_Jz().op_Jm().times_cst(sp.I/4)
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O22(myket):
    # 1/2 [JpJp + JmJm]
    a = myket.op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O22m(myket):
    # -i/2 [JpJp - JmJm]
    a = myket.op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O40(myket):
    mj = myket.Mj
    x = myket.J * (myket.J + 1)
    new_coef = 35 * mj**4 - (30*x -25)*mj*mj + 3*x*x - 6*x
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]

def sto_O41(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**3 - (3*x+1)*mj) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**3 - (3*x+1)*mj) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**3 - (3*x+1)*mj) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O41m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**3 - (3*x+1)*mj) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(-1*cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**3 - (3*x+1)*mj) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**3 - (3*x+1)*mj) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O42(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**2 - x -5) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**2 - x -5) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**2 - x -5) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O42m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (7*mj**2 - x -5) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (7*mj**2 - x -5) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (7*mj**2 - x -5) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O43(myket):
    # 1st term
    mj = myket.Mj
    cst = mj* sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * sp.Rational(1,4)
    mj = d.Mj
    d.coef *= mj * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O43m(myket):
    # 1st term
    mj = myket.Mj
    cst = mj * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * (-sp.I/4)
    mj = d.Mj
    d.coef *= mj * (sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O44(myket):
    # 1st term
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O44m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O60(myket):
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    new_coef = 231*mj**6 - (315*x-735)*mj**4 + (105*x*x-525*x+294)*mj*mj - 5*x**3 + 40*x*x - 60*x
    newket = deepcopy(myket)
    newket.coef = new_coef
    return [newket]

def sto_O61(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O61m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp()
    b = myket.times_cst(-1*cst).op_Jm()
    # 2nd term
    c = myket.op_Jp()
    d = myket.op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**5 - (30*x-15)*mj**3 + (5*x**2-10*x+12)*mj) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O62(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O62m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (33*mj**4 - (18*x+123)*mj**2 + x**2 + 10*x + 102) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O63(myket):
    # 1st term
    mj, x = myket.Mj, myket.J*(myket.J+1)
    cst = (11*mj**3 - (3*x+59)*mj) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J*(c.J+1)
    c.coef *= (11*mj**3 - (3*x+59)*mj) * sp.Rational(1,4)
    mj, x = d.Mj, d.J*(d.J+1)
    d.coef *= (11*mj**3 - (3*x+59)*mj) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O63m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J*(myket.J+1)
    cst = (11*mj**3 - (3*x+59)*mj) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J*(c.J+1)
    c.coef *= (11*mj**3 - (3*x+59)*mj) * (-sp.I/4)
    mj, x = d.Mj, d.J*(d.J+1)
    d.coef *= (11*mj**3 - (3*x+59)*mj) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O64(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (11*mj**2 -x - 38) * sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (11*mj**2 -x - 38) * sp.Rational(1,4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (11*mj**2 -x - 38) * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O64m(myket):
    # 1st term
    mj, x = myket.Mj, myket.J * (myket.J + 1)
    cst = (11*mj**2 -x - 38) * (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm()
    mj, x = c.Mj, c.J * (c.J + 1)
    c.coef *= (11*mj**2 -x - 38) * (-sp.I/4)
    mj, x = d.Mj, d.J * (d.J + 1)
    d.coef *= (11*mj**2 -x - 38) * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O65(myket):
    # 1st term
    mj = myket.Mj
    cst = mj* sp.Rational(1,4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(cst).op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * sp.Rational(1,4)
    mj = d.Mj
    d.coef *= mj * sp.Rational(1,4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O65m(myket):
    # 1st term
    mj = myket.Mj
    cst = mj* (-sp.I/4)
    a = myket.times_cst(cst).op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    b = myket.times_cst(-1*cst).op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    # 2nd term
    c = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp()
    d = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm()
    mj = c.Mj
    c.coef *= mj * (-sp.I/4)
    mj = d.Mj
    d.coef *= mj * (+sp.I/4)
    # sum equal terms and return
    a.coef += c.coef
    b.coef += d.coef
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O66(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(sp.Rational(1,2))
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(sp.Rational(1,2))
    return [idx for idx in [a, b] if idx.coef != 0]

def sto_O66m(myket):
    a = myket.op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().op_Jp().times_cst(-sp.I/2)
    b = myket.op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().op_Jm().times_cst(+sp.I/2)
    return [idx for idx in [a, b] if idx.coef != 0]
