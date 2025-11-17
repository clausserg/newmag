from magnetic_moments import *
import matplotlib.pyplot as plt
import os


def susceptibility_powder(mag_cntr, H_eff, T_max, output="cas", exp_file=None):
    """
    Compute powder susceptibility and chi*T vs T, optionally overlaying experimental data.

    Parameters
    ----------
    mag_cntr : object
        Magnetic center information.
    H_eff : ndarray
        Effective Hamiltonian.
    T_max : int
        Maximum temperature (K).
    output : str, optional
        Base name for output files.
    exp_file : str, optional
        Path to experimental data file with two columns: T(K) and chiT.
    """
    # let us get the magnetic moment operator matrices: -gJ*muB*J
    mu_ops_list = mu_operators_lmlsms(mag_cntr)  # gives a list of [mu_x, mu_y, mu_z]
    T_range = np.linspace(1, T_max, T_max)

    # Diagonalize the effective Hamiltonian
    evals, evecs = np.linalg.eigh(H_eff)
    N = len(evals)
    
    chi_powder = []
    for T in T_range:
        beta = 1.0 / (kB * T)

        shifted = evals - evals.min()  # work with relative energies for num stability
        weights = np.exp(-beta * shifted) 
        Z = np.sum(weights)  # statistical sum

        chi_dirs = []  # x y, z directions
        for mu in mu_ops_list:
            mu_mat = evecs.T.conj() @ mu @ evecs  # Transform mu into eigenbasis of Heff
            chi = 0.0

            for n in range(N):
                for m in range(N):
                    if n == m:  # Curie contribution
                        chi += (weights[n] / (kB * T * Z)) * (abs(mu_mat[n, n])**2)
                    else:  # Van Vleck contribution
                        denom = evals[m] - evals[n]
                        if abs(denom) > 1e-12:
                            chi += (weights[n] - weights[m]) * (abs(mu_mat[n, m])**2) / (Z * denom)

            chi_dirs.append(chi)

        # Powder average
        chi_avg = np.mean(chi_dirs)  # get the average
        chi_avg *= (muB_NA / muB)  # Convert to cm^3/mol
        chi_powder.append(chi_avg)

    # get chi*T
    chiT = np.array(chi_powder) * T_range

    # ----------------------------
    # Plot
    # ----------------------------
    plt.figure(figsize=(3.3, 2.5))  # one-column size in journals (~8.4 cm)
    plt.plot(T_range, chiT, color='black', label="NewMag", linewidth=1.0)

    # If experimental file provided, overlay it
    if exp_file is not None and os.path.exists(exp_file):
        exp_data = np.loadtxt(exp_file)
        if exp_data.shape[1] >= 2:  # need at least two columns
            T_exp, chiT_exp = exp_data[:, 0], exp_data[:, 1]
            plt.plot(T_exp, chiT_exp, color='blue', label="Orca")

    plt.xlabel("T (K)", fontsize=11)
    plt.ylabel(r"$\chi T$ (cm$^{3}$ K mol$^{-1}$)", fontsize=11)
    plt.legend(frameon=False, fontsize=9)

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='in', top=False, right=False)

    plt.tight_layout()
    fig_out = f"chiT_vs_T_{output}.png"
    plt.savefig(fig_out, dpi=600)  # high-resolution PNG
    # plt.show()

    # ----------------------------
    # Save to TXT (overwrite if exists)
    # ----------------------------
    txt_out = f"chiT_vs_T_{output}.txt"
    if os.path.exists(txt_out):
        os.remove(txt_out)

    with open(txt_out, "w") as f:
        f.write("T(K)\tchiT(emu K mol^-1)\n")
        for T, val in zip(T_range, chiT):
            f.write(f"{T:.6f}\t{val:.6f}\n")

    print(f"Data saved to {txt_out} and plot saved to {fig_out}")

