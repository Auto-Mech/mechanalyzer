"""
Calculate rates with various fitting functions
"""

import numpy as np
from scipy.special import eval_chebyt


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K
RC2 = 0.0820573660809596  # Gas Constant in L.atm/mol.K


def single_arrhenius(a_par, n_par, ea_par,
                     t_ref, temp):
    """ calc value with single arrhenius function
    """
    kts = a_par * ((temp / t_ref)**n_par) * np.exp(-ea_par/(RC*temp))
    return kts


def double_arrhenius(a_par1, n_par1, ea_par1,
                     a_par2, n_par2, ea_par2,
                     t_ref, temp):
    """ calc value with single arrhenius function
    """
    kts = (
        a_par1 * ((temp / t_ref)**n_par1) * np.exp(-ea_par1/(RC*temp)) +
        a_par2 * ((temp / t_ref)**n_par2) * np.exp(-ea_par2/(RC*temp))
    )
    return kts


def arrhenius(params, t_ref, temp):
    """ simplified function whcih will call single or double arrhenis
        based on the number of params that are passed in
        only does single and double fits.
        params must be [a1, n1, ea1] or [a1, n2, ea1, a2, n2, ea2]
    """
    assert len(params) in (3, 6)

    if len(params) == 3:
        kts = single_arrhenius(
            params[0], params[1], params[2],
            t_ref, temp)
    else:
        kts = double_arrhenius(
            params[0], params[1], params[2],
            params[3], params[4], params[5],
            t_ref, temp)

    return kts


def lindemann(highp_ks, lowp_ks, pressures, temps):
    """ calculate pressure-dependence constants according to Lindemann
        model; no value for high
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = lindemann_rate_constants(
            highp_ks, lowp_ks, pressure, temps)

    return ktp_dct


def lindemann_rate_constants(highp_ks, lowp_ks, pressure, temps):
    """ calculate pressure-dependence constants according to Lindemann
        model
    """
    # Calculate the pr term
    pr_term = _pr_term(highp_ks, lowp_ks, pressure, temps)

    # Calculate Lindemann rate constants
    ktps = highp_ks * (pr_term / (1.0 + pr_term))

    return ktps


def troe(highp_ks, lowp_ks, pressures, temps,
         alpha, ts3, ts1, ts2=None):
    """ calculate pressure-dependence constants according to Troe
        model; no value for high
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = troe_rate_constants(
            highp_ks, lowp_ks, pressure, temps,
            alpha, ts3, ts1, ts2)

    return ktp_dct


def troe_rate_constants(highp_ks, lowp_ks, pressure, temp,
                        alpha, ts3, ts1, ts2=None):
    """ calculate pressure-dependence constants according to Troe
        model
    """
    # Calculate the pr term and broadening factor
    pr_term = _pr_term(highp_ks, lowp_ks, pressure, temp)
    f_term = _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temp)
    # Calculate Troe rate constants
    ktps = highp_ks * (pr_term / (1.0 + pr_term)) * f_term

    return ktps


def plog(plog_dct, t_ref, pressures, temps):
    """ calculate the rate constant using a dictionary of plog params
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = plog_rate_constants(
            plog_dct, t_ref, pressure, temps)

    return ktp_dct


def plog_rate_constants(plog_dct, t_ref, pressure, temps):
    """ calculate the rate constant using a dictionary of plog params
    """
    plog_pressures = [pressure for pressure in plog_dct]

    # Check if pressure is in plog dct; use plog pressure for numerical stab
    pressure_defined = False
    for plog_pressure in plog_pressures:
        if np.isclose(pressure, plog_pressure, atol=1.0e-3):
            pressure_defined = True
            plog_params = plog_dct[plog_pressure]

    # print('plog test')
    # print(plog_dct)
    # print(pressure)
    # print(plog_pressures)

    # If pressure equals value use, arrhenius expression
    if pressure_defined:
        ktps = arrhenius(plog_params, t_ref, temps)
    # Find which two PLOG pressures our pressure of interest sits between
    else:
        for i, _ in enumerate(plog_pressures):
            if i != len(plog_pressures)-1:
                if plog_pressures[i] < pressure < plog_pressures[i+1]:
                    plow = plog_pressures[i]
                    phigh = plog_pressures[i+1]
                    plow_params = plog_dct[plow]
                    phigh_params = plog_dct[phigh]
                    break

        kt_low = arrhenius(plow_params, t_ref, temps)
        kt_high = arrhenius(phigh_params, t_ref, temps)
        pres_term = (
            (np.log10(pressure) - np.log10(plow)) /
            (np.log10(phigh) - np.log10(plow))
        )
        logkt = (
            np.log10(kt_low) +
            ((np.log10(kt_high) - np.log10(kt_low)) * pres_term)
        )

        ktps = 10**(logkt)

    return ktps


def chebyshev(alpha, tmin, tmax, pmin, pmax, pressures, temps):
    """ computes the rate constants using the chebyshev polynomials
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = chebyshev_rate_constants(
            temps, pressure, alpha, tmin, tmax, pmin, pmax)

    return ktp_dct


def chebyshev_rate_constants(temps, pressure, alpha, tmin, tmax, pmin, pmax):
    """ computes the rate constants using the chebyshev polynomials
    """
    alpha_nrows, alpha_ncols = alpha.shape

    ktps = np.zeros(len(temps))
    for i, temp in enumerate(temps):
        ctemp = (
            (2.0 * temp**(-1) - tmin**(-1) - tmax**(-1)) /
            (tmax**(-1) - tmin**(-1))
        )
        cpress = (
            (2.0 * np.log10(pressure) - np.log10(pmin) - np.log10(pmax)) /
            (np.log10(pmax) - np.log10(pmin))
        )

        logktp = 0.0
        for j in range(alpha_nrows):
            for k in range(alpha_ncols):
                logktp += (
                    alpha[j][k] *
                    eval_chebyt(j, ctemp) *
                    eval_chebyt(k, cpress)
                )

        ktps[i] = 10**(logktp)

    return ktps


def _pr_term(highp_rateks, lowp_rateks, pressure, temp):
    """ calculate the corrective pr term used for Lindemann and Troe
        pressure-dependent forms
    """
    pr_term = (lowp_rateks / highp_rateks) * (pressure / (RC2 * temp))
    return pr_term


def _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temp):
    """ calculate the F broadening factor used for Troe
        pressure-dependent forms
    """

    # Calculate Fcent term
    f_cent = ((1.0 - alpha) * np.exp(-temp / ts3) +
              alpha * np.exp(-temp / ts1))
    if ts2 is not None:
        f_cent += np.exp(-ts2 / temp)

    # Calculate the Log F term
    c_val = -0.4 - 0.67 * np.log10(f_cent)
    n_val = 0.75 - 1.27 * np.log10(f_cent)
    d_val = 0.14
    val = ((np.log10(pr_term) + c_val) /
           (n_val - d_val * (np.log10(pr_term) + c_val)))**2
    logf = (1.0 + val)**(-1) * np.log10(f_cent)

    # Calculate F broadening term
    f_term = 10**(logf)

    return f_term
