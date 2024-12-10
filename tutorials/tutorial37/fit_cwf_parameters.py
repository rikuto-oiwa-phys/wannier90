import sys
import os
import numpy as np
import lmfit

cwd = os.getcwd()


# ==================================================
def read_eig(seedname):
    """
    read seedname.eig file.

    Args:
        seedname (str, optional): seedname.

    Returns:
        ndarray: Kohn-Sham energies, E_{n}^{k}.
    """
    file_name = cwd + "/" + seedname + ".eig"
    if os.path.exists(file_name):
        with open(file_name) as fp:
            eig_data = fp.readlines()
    else:
        raise Exception("failed to read eig file: " + file_name)

    eig_data = [[v for v in lst.rstrip("\n").split(" ") if v != ""] for lst in eig_data]
    eig_data = [[float(v) if "." in v else int(v) for v in lst] for lst in eig_data]

    num_bands = np.max([v[0] for v in eig_data])
    num_k = np.max([v[1] for v in eig_data])
    Ek = np.array([[eig_data[k * num_bands + m][2] for m in range(num_bands)] for k in range(num_k)])

    return Ek


# ==================================================
def read_amn(seedname):
    """
    read seedname.amn file.

    Args:
        seedname (str, optional): seedname.

    Returns:
        ndarray: Overlap matrix elements, A_{mn}^{k}.
    """
    file_name = cwd + "/" + seedname + ".amn"
    if os.path.exists(file_name):
        with open(file_name) as fp:
            amn_data = fp.readlines()
    else:
        raise Exception("failed to read amn file: " + file_name)

    num_bands, num_k, num_wann = [int(x) for x in amn_data[1].split()]
    amn_data = np.genfromtxt(amn_data[2:]).reshape(num_k, num_wann, num_bands, 5)
    Ak = np.transpose(amn_data[:, :, :, 3] + 1j * amn_data[:, :, :, 4], axes=(0, 2, 1))

    return Ak


# ==================================================
def fermi_dist_func(x, kbT):
    if kbT == 0.0:
        return 1.0 * np.vectorize(float)(x < 0.0)

    return 0.5 * (1.0 - np.tanh(0.5 * x / kbT))


# ==================================================
def cwf_window_function(e, mu_min, mu_max, sigma_min, sigma_max):
    """Returns the window function according to the closest Wannier method"""
    delta = 10e-12
    return fermi_dist_func(mu_min - e, sigma_min) + fermi_dist_func(e - mu_max, sigma_max) - 1.0 + delta


# ==================================================
def optimize_params(seedname):
    """
    optimize the energy windows and smearing temperatures
    by using the projectability of each Kohn-Sham state in k-space.

    Args:
        seedname (str, optional): seedname.
    """
    Ek = read_eig(seedname)
    Ak = read_amn(seedname)

    num_k, num_bands, num_wann = Ak.shape

    # projectability of each Kohn-Sham state in k-space.
    Pk = np.real(np.diagonal(Ak @ Ak.transpose(0, 2, 1).conjugate(), axis1=1, axis2=2))

    ek = Ek.reshape(num_k * num_bands)
    pk = Pk.reshape(num_k * num_bands)

    # normalize
    pk = pk / np.max(pk)

    model = lmfit.Model(cwf_window_function)
    params = lmfit.Parameters()

    params_init_dict = {
        "mu_min": np.min(ek),
        "mu_max": np.max(ek),
        "sigma_min": 0.0,
        "sigma_max": 1.0,
    }

    params_bound_dict = {
        "mu_min": (-np.inf, np.inf),
        "mu_max": (-np.inf, np.inf),
        "sigma_min": (0, 100),
        "sigma_max": (0, 100),
    }

    for param in params_bound_dict.keys():
        init = params_init_dict[param]
        min, max = params_bound_dict[param]
        params.add(param, value=init, min=min, max=max)

    result = model.fit(pk, params, e=ek)

    mu_min_opt = result.params["mu_min"].value
    mu_max_fit = result.params["mu_max"].value
    sigma_min_opt = result.params["sigma_min"].value
    sigma_max_opt = result.params["sigma_max"].value

    mu_max_opt = mu_max_fit - 3.0 * sigma_max_opt

    msg = "# Optimized parameters: \n"
    msg += f"cwf_mu_max       = {mu_max_opt} = mu_max_fit - 3*sigma_max_opt \n"
    msg += f"cwf_mu_min       = {mu_min_opt}  \n"
    msg += f"cwf_sigma_max    = {sigma_max_opt}  \n"
    msg += f"cwf_sigma_min    = {sigma_min_opt} \n"
    
    f = open('cwf_parameters.txt', 'w')
    f.write(msg)
    f.close()
    
    print(msg)

    msg = "# energy [eV]   projectability   window function (fit)  window function (opt)\n"
    for ekn, pkn in zip(ek, pk):
        wkn_fit = cwf_window_function(ekn, mu_min_opt, mu_max_fit, sigma_min_opt, sigma_max_opt)
        wkn_opt = cwf_window_function(ekn, mu_min_opt, mu_max_opt, sigma_min_opt, sigma_max_opt)
        msg += f"{ekn}   {pkn}   {wkn_fit}   {wkn_opt} \n"

    f = open('p_w_vs_e.dat', 'w')
    f.write(msg)
    f.close()


# ==================================================
if __name__ == "__main__":
    seedname = sys.argv[1]
    optimize_params(seedname)
