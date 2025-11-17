def matrix_soc(basis, zeta):
    """
    Compute the Spin-Orbit (SO) matrix for the given basis set and SOC constant.

    Parameters:
        basis: A list of LMLSMS or JMJ kets.
        zeta: The spin-orbit coupling constant.

    Returns:
        SO_mat: The resulting Spin-Orbit matrix.
    """
    # If basis consists of LMLSMS kets, return a zero matrix immediately
    if isinstance(basis[0], LMLSMS):
        return sp.zeros(len(basis), len(basis))

    # Initialize the SO matrix
    SO_mat = sp.zeros(len(basis), len(basis))

    # Precompute op_soc for each bket in the basis
    soc_dict = {idx: op_soc(basis[idx]) for idx in range(len(basis))}

    # Loop over all pairs of basis elements
    for idx in range(len(basis)):
        for jdx in range(len(basis)):
            # Get the corresponding SOC expansion for the ket
            soc = soc_dict[jdx]  # Precomputed SOC terms for bket
            SO_mat[idx, jdx] = inner_prod(basis[idx], soc)  # Compute inner product

    # Return the matrix scaled by zeta
    return SO_mat * zeta
