from class_jmj import JMJ
from class_lmlsms import LMLSMS
def inner_prod(bas_vec, myvec):
    """
    Compute the inner product between bas_vec and myvec.

    Parameters:
        bas_vec: A JMJ or LMLSMS ket (treated as 'bra').
        myvec: A list of JMJ or LMLSMS kets (treated as 'ket').

    Returns:
        result: The inner product value.
    """
    result = 0

    if isinstance(bas_vec, JMJ) and isinstance(myvec, JMJ):
        if bas_vec.Ket == myvec.Ket:
            result += bas_vec.coef.conjugate() * myvec.coef
        return result
    
    if isinstance(bas_vec, JMJ) and isinstance(myvec[0], LMLSMS):
        for idx in bas_vec.lmlsms_expansion:
            for jdx in myvec:
                if idx.Ket == jdx.Ket:
                    result += idx.coef.conjugate() * jdx.coef
        return result
                    
    if isinstance(bas_vec, LMLSMS) and isinstance(myvec[0], LMLSMS):
        for idx in myvec:
            if idx.Ket == bas_vec.Ket:
                result += bas_vec.coef.conjugate() * idx.coef
        return result
        
    if isinstance(bas_vec, JMJ) and isinstance(myvec[0], JMJ):
        for idx in myvec:
            if bas_vec.Ket == idx.Ket:
                result += bas_vec.coef.conjugate() * idx.coef
        return result
