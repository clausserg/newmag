def op_soc(jmj_ket):
    """
    Compute the spin-orbit coupling terms for a given JMJ ket.
    """
    def lz_sz(term):
        """Compute the Laz * Saz term."""
        return term.op_Sz().op_Jz()

    def lp_sm(term):
        """Compute the 1/2 (La+ * Sa-) term."""
        a = term.op_Sm().op_Jp()
        a.coef *= sp.Rational(1, 2)
        return a

    def lm_sp(term):
        """Compute the 1/2 (La- * Sa+) term."""
        a = term.op_Sp().op_Jm()
        a.coef *= sp.Rational(1, 2)
        return a

    # Generate all terms and filter out zero-coefficient results in a single step
    res = [
        term for op in (lz_sz, lp_sm, lm_sp) 
        for term in (op(t) for t in jmj_ket.lmlsms_expansion) 
        if term.coef != 0
    ]
    return res
