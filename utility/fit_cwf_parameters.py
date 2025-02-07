#! /usr/bin/env python3
#
# This file is distributed as part of the Wannier90 code and
# under the terms of the GNU General Public License. See the
# file `LICENSE' in the root directory of the Wannier90
# distribution, or http://www.gnu.org/copyleft/gpl.txt
#
# The webpage of the Wannier90 code is www.wannier.org
#
# The Wannier90 code is hosted on GitHub:
#
# https://github.com/wannier-developers/wannier90
#
# Python3 script to fit the parameters for the closest Wannier method
# using the projectabilities calculated by Amn matrices.
# This method can be used if the pseudo-atomic orbitals (PAOs)
# are used for the initial guesses.
# For the detail of the theoretical background,
# see npj Comput. Mater. 6, 66 (2020).
#
# Written by Rikuto Oiwa (RIKEN)
# Last update December 11th, 2024
#
import os
import sys
import lmfit
import numpy as np

cwd = os.getcwd()

# ==================================================
def read_eig(seedname):
    """
    read seedname.eig file.

    Args:
        seedname (str): seedname.

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
        seedname (str): seedname.

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
def fermi_dist_func(e, kbT):
    """Returns the Fermi-dirac distribution function"""
    if kbT == 0.0:
        return 1.0 * np.vectorize(float)(e < 0.0)

    return 0.5 * (1.0 - np.tanh(0.5 * e / kbT))


# ==================================================
def cwf_window_function(e, mu_min, mu_max, sigma_min, sigma_max):
    """Returns the window function according to the closest Wannier method"""
    delta = 10e-12
    return fermi_dist_func(mu_min - e, sigma_min) + fermi_dist_func(e - mu_max, sigma_max) - 1.0 + delta


# ==================================================
def fit_cwf_parameters(seedname):
    """
    fit the parameters for the closest Wannier method
    using the projectability of each Kohn-Sham state in k-space.

    Args:
        seedname (str, optional): seedname.
    """
    Ek = read_eig(seedname)
    Ak = read_amn(seedname)
    # projectability of each Kohn-Sham state in k-space.
    Pk = np.real(np.diagonal(Ak @ Ak.transpose(0, 2, 1).conjugate(), axis1=1, axis2=2))

    num_k, num_bands, num_wann = Ak.shape

    ek = Ek.reshape(num_k * num_bands)
    pk = Pk.reshape(num_k * num_bands)

    # normalize
    # pk = pk / np.max(pk)

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

    mu_min_fit = result.params["mu_min"].value
    mu_max_fit = result.params["mu_max"].value
    sigma_min_fit = result.params["sigma_min"].value
    sigma_max_fit = result.params["sigma_max"].value

    # See npj Comput. Mater. 74, 195118 (2006)
    mu_min_opt = mu_min_fit
    mu_max_opt = mu_max_fit - 3.0 * sigma_max_fit
    sigma_min_opt = sigma_min_fit
    sigma_max_opt = sigma_max_fit

    msg = "# Optimized parameters: \n"
    msg += f"cwf_mu_max       = {mu_max_opt} (mu_max_fit - 3*sigma_max_fit) \n"
    msg += f"cwf_mu_min       = {mu_min_opt}  \n"
    msg += f"cwf_sigma_max    = {sigma_max_opt}  \n"
    msg += f"cwf_sigma_min    = {sigma_min_opt} \n"

    f = open('cwf_parameters.dat', 'w')
    f.write(msg)
    f.close()

    msg = "# energy [eV]           projectability        window func. (fit)    window func. (opt) \n"
    for ekn, pkn in zip(ek, pk):
        wkn_fit = cwf_window_function(ekn, mu_min_opt, mu_max_fit, sigma_min_opt, sigma_max_opt)
        wkn_opt = cwf_window_function(ekn, mu_min_opt, mu_max_opt, sigma_min_opt, sigma_max_opt)
        msg += " {0[0]:18.15f}    {0[1]:18.15f}    {0[2]:18.15f}    {0[3]:18.15f} \n".format([ekn, pkn, wkn_fit, wkn_opt])

    f = open('p_w_vs_e.dat', 'w')
    f.write(msg)
    f.close()

    msg = "optimal parameters are written in cwf_parameters.dat file. \n"
    msg += "projectabilities and the fitting window function are written in p_w_vs_e.dat file. \n"

    print(msg)


# ==================================================
if __name__ == "__main__":
    seedname = sys.argv[1]
    fit_cwf_parameters(seedname)

