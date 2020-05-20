"""
Calculate rates with various fitting functions
"""

import numpy as np
from scipy.special import eval_chebyt


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K
RC2 = 0.0820573660809596  # Gas Constant in L.atm/mol.K


def single_arrhenius(a_par, n_par, ea_par,
                     t_ref, temps, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using
        a single Arrhenius functional expression:
            k(T) = A*[(T/T_ref)^n]*exp(-Ea/RT)
        :param float a_par: pre-exponential A parmater
        :param float n_par: temperature exponent n parmater
        :param float ea_par: activation energy Ea parmater (kcal.mol-1)
        :param float t_ref: Reference temperature (K)
        :param numpy.ndarray temps: List of Temperatures (K)
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """
    kts = a_par * ((temps / t_ref)**n_par) * np.exp(-ea_par/(rval*temps))
    return kts


def double_arrhenius(a_par1, n_par1, ea_par1,
                     a_par2, n_par2, ea_par2,
                     t_ref, temps, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using
        a double Arrhenius functional expression:
            k(T) = { A1*[(T/T_ref)^n1]*exp(-Ea1/RT) +
                     A2*[(T/T_ref)^n2]*exp(-Ea2/RT)}
        :param float a_par1: 1st pre-exponential A parmater
        :param float n_par1: 1st temperature exponent n parmater
        :param float ea_par1: 1st activation energy Ea parmater (kcal.mol-1)
        :param float a_par2: 2nd pre-exponential A parmater
        :param float n_par2: 2nd temperature exponent n parmater
        :param float ea_par2: 2nd activation energy Ea parmater (kcal.mol-1)
        :param float t_ref: Reference temperature (K)
        :param numpy.ndarray temps: List of Temperatures (K)
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """
    kts = (
        a_par1 * ((temps / t_ref)**n_par1) * np.exp(-ea_par1/(rval*temps)) +
        a_par2 * ((temps / t_ref)**n_par2) * np.exp(-ea_par2/(rval*temps))
    )
    return kts


def arrhenius(params, t_ref, temp):
    """ Calculates T-dependent rate constants [k(T)]s using
         a either a single or double Arrhenius functional expression,
         depending on the number of input parameters.
         params must be [a1, n1, ea1] or [a1, n2, ea1, a2, n2, ea2]
         :param list params: pre-exponential A parmater
         :param float t_ref: Reference temperature (K)
         :param numpy.ndarray temp: List of Temperatures (K)
         :return kts: T-dependent rate constants
         :rtype: numpy.ndarray
    """

    assert len(params) in (1, 2)
    for param_set in params:
        assert len(param_set) == 3

    if len(params) == 1:
        kts = single_arrhenius(
            params[0][0], params[0][1], params[0][2],
            t_ref, temp)
    else:
        kts = double_arrhenius(
            params[0][0], params[0][1], params[0][2],
            params[1][0], params[1][1], params[1][2],
            t_ref, temp)

    return kts


def lindemann(highp_ks, lowp_ks, pressures, temps):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Lindemann functional expression:
          k(T,P) = kinf(T) * [ Pr / 1 + Pr ]
          where Pr is the reduced pressure.
        :param list highp_ks: k(T)s determined at high-pressure
        :param list lowp_ks: k(T)s determined at low-pressure
        :param list pressures: Pressures used to calculate k(T,P)s
        :param numpy.ndarray temps: Temps used to calculate high- and low-k(T)s
        :return ktp_dct: T-dependent rate constants
        :rtype: dct {Pressure1: TempsArray1, Pressure2: TempsArray2, ...}
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
         alpha, ts3, ts1, ts2=None, collid_factor=1.0):
    """ calculate pressure-dependence constants according to Troe
        model; no value for high
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = troe_rate_constants(
            highp_ks, lowp_ks, pressure, temps,
            alpha, ts3, ts1, ts2=ts2, collid_factor=collid_factor)

    return ktp_dct


def troe_rate_constants(highp_ks, lowp_ks, pressure, temp,
                        alpha, ts3, ts1, ts2=None, collid_factor=1.0):
    """ calculate pressure-dependence constants according to Troe
        model
    """
    # Calculate the pr term and broadening factor
    pr_term = _pr_term(highp_ks, lowp_ks, pressure, temp)
    f_term = _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temp)
    # Calculate Troe rate constants
    ktps = highp_ks * (pr_term / (1.0 + pr_term)) * f_term * collid_factor

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
    plog_pressures = list(plog_dct.keys())

    # Check if pressure is in plog dct; use plog pressure for numerical stab
    pressure_defined = False
    for plog_pressure in plog_pressures:
        if np.isclose(pressure, plog_pressure, atol=1.0e-3):
            pressure_defined = True
            plog_params = plog_dct[plog_pressure]

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


def _pr_term(highp_rateks, lowp_rateks, pressure, temps, rval=RC2):
    """ Calculates the reduced pressure term for a single pressure:
           [k0(T) / kinf(T)] * [P / RT]
        used for Lindemann and Troe P-dependent functional expressions
        :param list highp_ks: k(T)s determined at high-pressure
        :param list lowp_ks: k(T)s determined at low-pressure
        :param list pressure: Pressures used to calculate k(T,P)
        :param numpy.ndarray temps: Temps used to calculate high- and low-k(T)s
        :return ktp_dct: T-dependent rate constants
        :rtype: dct {Pressure1: TempsArray1, Pressure2: TempsArray2, ...}
    """
    pr_term = (lowp_rateks / highp_rateks) * (pressure / (rval * temps))
    return pr_term


def _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temp):
    """ Calculates the F broadening factor term
           log F = { 1 + [ (logPr + c)/(n-d(logPr)+c)) ]^2 }^(-1)
           c = -0.4 - 0.67*logFcent
           n = 0.75 - 1.27*logFcent
           d = 0.14
           Fcent = (1-alpha)*exp(-T/ts3)+
                   alpha*exp(-T/ts1)+
                   exp(-ts2/T)
        :param float pr_term: reduced pressure
        :param float alpha: Troe alpha parameter
        :param float ts3: Troe T3 parameter
        :param float ts1: Troe T1 parameter
        :param float ts2: Troe T2 parameter
        :param float temp: temperature
        :return f_cent: F broadening term
        :rtype: float
        used for Troe P-dependent functional expressions
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
