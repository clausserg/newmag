import sympy as sp


"""
Sympy symbols used in the construction of the Stevens CF
and Zeeman matrices across the NewMag project
"""

# SO constant
Z = sp.symbols("Z",real=True)

# Magnetic couplings
J0, J1, J2, J3 = sp.symbols("J_0 J_1 J_2 J_3",real=True)

# Stevens CF parameters
B20, B21, B22 = sp.symbols("B_2^0 B_2^1 B_2^2",real=True)
B21m, B22m = sp.symbols("B_2^-1 B_2^-2",real=True)

B40, B41, B42, B43, B44 = sp.symbols("B_4^0 B_4^1 B_4^2 B_4^3 B_4^4", real=True)
B41m, B42m, B43m, B44m = sp.symbols("B_4^-1 B_4^-2, B_4^-3 B_4^-4", real=True)

B50, B51, B52, B53, B54, B55 = sp.symbols("B_5^0 B_5^1 B_5^2 B_5^3 B_5^4 B_5^5", real=True)
B51m, B52m, B53m, B54m, B55m = sp.symbols("B_5^-1 B_5^-2, B_5^-3 B_5^-4 B_5^-5", real=True)

B60, B61, B62, B63, B64, B65, B66 = sp.symbols("B_6^0 B_6^1 B_6^2 B_6^3 B_6^4 B_6^5 B_6^6", real=True)
B61m, B62m, B63m, B64m, B65m, B66m = sp.symbols("B_6^-1 B_6^-2, B_6^-3 B_6^-4 B_6^-5 B_6^-6", real=True)

# g-factors
gxx, gyy, gzz, gxy, gyx, gxz, gzx, gyz, gzy = sp.symbols("g_xx g_yy g_zz g_xy g_yx g_xz g_zx g_yz g_zy", real=True)

# Magnetic field in x, y, z directions
Bx, By, Bz = sp.symbols("B_x B_y B_z", real=True)

# Lande g-factor
gJ = sp.symbols("g_J", real=True)

# Bohr magneton
BM = sp.symbols("\\mu_B", real=True)
