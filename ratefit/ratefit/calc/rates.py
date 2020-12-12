"""
Calculate rates with various fitting functions
"""

import numpy as np
from scipy.special import eval_chebyt
from lib.phydat import phycon


RC = phycon.RC_cal  # gas constant in cal/(mol.K) 
RC2 = phycon.RC_atm  # gas constant in cm^3.atm/(mol.K) 


def single_arrhenius(a_par, n_par, ea_par,
                     t_ref, temps, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using
        a single Arrhenius functional expression.

        :param a_par: pre-exponential A parmater
        :type a_par: float
        :param n_par: temperature exponent n parmater
        :type n_par: float
        :param ea_par: activation energy Ea parmater (kcal.mol-1)
        :type ea_par: float
        :param t_ref: Reference temperature (K)
         :type t_ref: float
        :param temps: List of Temperatures (K)
        :type temps: numpy.ndarray
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """
    kts = a_par * ((temps / t_ref)**n_par) * np.exp(-ea_par/(rval*temps))
    return kts


def double_arrhenius(a_par1, n_par1, ea_par1,
                     a_par2, n_par2, ea_par2,
                     t_ref, temps, rval=RC):
    """ Calculates T-dependent rate constants [k(T)]s using
        a double Arrhenius functional expression.

        :param a_par1: 1st pre-exponential A parmater
        :type a_par1: float
        :param n_par1: 1st temperature exponent n parmater
        :type n_par1: float
        :param ea_par1: 1st activation energy Ea parmater (kcal.mol-1)
        :type ea_par1: float
        :param a_par2: 2nd pre-exponential A parmater
        :type a_par2: float
        :param n_par2: 2nd temperature exponent n parmater
        :type n_par2: float
        :param ea_par2: 2nd activation energy Ea parmater (kcal.mol-1)
        :type ea_par2: float
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: List of Temperatures (K)
        :type temps: numpy.ndarray
        :return kts: T-dependent rate constants
        :rtype: numpy.ndarray
    """
    kts = (
        a_par1 * ((temps / t_ref)**n_par1) * np.exp(-ea_par1/(rval*temps)) +
        a_par2 * ((temps / t_ref)**n_par2) * np.exp(-ea_par2/(rval*temps))
    )
    return kts


def arrhenius(params, t_ref, temps):
    """ Calculates T-dependent rate constants [k(T)]s using
         a either a single or double Arrhenius functional expression,
         depending on the number of input fitting parameters.

         :param params: fitting parameters
         :type params: list(float)
         :param t_ref: Reference tempserature (K)
         :type t_ref: float
         :param temps: List of Temperatures (K)
         :type temps: numpy.ndarray
         :return kts: T-dependent rate constants
         :rtype: numpy.ndarray
    """

    assert len(params) in (1, 2)
    for param_set in params:
        assert len(param_set) == 3

    if len(params) == 1:
        kts = single_arrhenius(
            params[0][0], params[0][1], params[0][2],
            t_ref, temps)
    else:
        kts = double_arrhenius(
            params[0][0], params[0][1], params[0][2],
            params[1][0], params[1][1], params[1][2],
            t_ref, temps)

    return kts


def lowp_limit(highp_rateks, temps, pressures, collid_factor=1.0, rval=RC2):
    """ Calculates T,P-dependent rate constants [k(T,P)]s assuming
        the reaction occurs in the low-pressure regime where the
        rates are linear with pressure.

        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = lowp_limit_one_pressure(
            highp_rateks, temps, pressure,
            collid_factor=collid_factor, rval=rval)

    return ktp_dct


def lowp_limit_one_pressure(highp_rateks, temps, pressure,
                            collid_factor=1.0, rval=RC2):
    """ Calculates the reduced pressure term for a single pressure
        used for Lindemann and Troe P-dependent functional expressions.

        :param list highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :temps: numpy.ndarray
        :param pressure: Pressure used to calculate reduced pressure
        :type pressure: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :rtype: numpy.ndarray
    """
    return highp_rateks * p_to_m(pressure, temps, rval=rval) * collid_factor


def lindemann(highp_ks, lowp_ks, temps, pressures, collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Lindemann functional expression.

        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = lindemann_one_pressure(
            highp_ks, lowp_ks, temps, pressure, collid_factor=collid_factor)

    return ktp_dct


def lindemann_one_pressure(highp_ks, lowp_ks, temps, pressure,
                           collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Lindemann functional expression, at a given pressure,
        across several temperatures.

        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressure: Pressure used to calculate k(T,P)s
        :type pressure: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
    """
    # Calculate the pr term
    pr_terms = _pr_term(highp_ks, lowp_ks, temps, pressure,
                        collid_factor=collid_factor)

    # Calculate Lindemann rate constants
    ktps = highp_ks * (pr_terms / (1.0 + pr_terms))

    return ktps


def troe(highp_ks, lowp_ks, temps, pressures,
         alpha, ts3, ts1, ts2=None, collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Troe functional expression.

        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressure: list(float)
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = troe_one_pressure(
            highp_ks, lowp_ks, temps, pressure,
            alpha, ts3, ts1, ts2=ts2, collid_factor=collid_factor)

    return ktp_dct


def troe_one_pressure(highp_ks, lowp_ks, temps, pressure,
                      alpha, ts3, ts1, ts2=None, collid_factor=1.0):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Troe functional expression, at a given pressure,
        across several temperatures.

        :param highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressure used to calculate k(T,P)s
        :type pressure: float
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param collid_factor: Buffer enhancement collision factor
        :type collid_factor: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
        :return ktps: Set of k(T,P)s at given pressure
        :rtype: dict[pressure: temps]
    """

    # Calculate the pr term and broadening factor
    pr_term = _pr_term(highp_ks, lowp_ks, temps, pressure, collid_factor)
    f_term = _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temps)

    # Calculate Troe rate constants (collision factor could be wrong)
    ktps = highp_ks * (pr_term / (1.0 + pr_term)) * f_term

    return ktps


def plog(plog_dct, t_ref, temps, pressures):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a PLOG functional expression.

        :param plog_dct: Arrhenius fitting parameters at several pressures
        :type plog_dct: dict[pressure: [fit_params]]
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """

    # Set the plog pressures to see if pressure is in range
    plog_pressures = list(plog_dct.keys())

    ktp_dct = {}
    for pressure in pressures:
        if min(plog_pressures) <= pressure <= max(plog_pressures):
            ktp_dct[pressure] = plog_one_pressure(
                plog_dct, t_ref, temps, pressure)

    return ktp_dct


def plog_one_pressure(plog_dct, t_ref, temps, pressure):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a PLOG functional expression, at a given pressure,
        across several temperatures.

        :param plog_dct: Arrhenius fitting parameters at several pressures
        :type plog_dct: dict[pressure: [fit_params]]
        :param t_ref: Reference temperature (K)
        :type t_ref: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressure: Pressure used to calculate k(T,P)s
        :type pressure: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
    """

    plog_pressures = list(plog_dct.keys())

    # Check if pressure is in plog dct; use plog pressure for numerical stab
    pressure_defined = False
    for plog_pressure in plog_pressures:
        if np.isclose(pressure, plog_pressure, atol=1.0e-3):
            pressure_defined = True
            plog_params = plog_dct[plog_pressure]

    # If pressure equals value use, arrhenius expression
    # pressure_defined = False
    # print('def', pressure_defined)
    if pressure_defined:
        ktps = arrhenius(plog_params, t_ref, temps)
    else:
        # Calculate pressure term for PLOG expression
        # Use two PLOG pressures our pressure of interest sits between
        for i, _ in enumerate(plog_pressures):
            if i != len(plog_pressures)-1:
                if plog_pressures[i] < pressure < plog_pressures[i+1]:
                    plow = plog_pressures[i]
                    phigh = plog_pressures[i+1]
                    plow_params = plog_dct[plow]
                    phigh_params = plog_dct[phigh]
                    break
        pres_term = (
            (np.log10(pressure) - np.log10(plow)) /
            (np.log10(phigh) - np.log10(plow))
        )

        # Calculate k(T)s at high-P and low-P with Arrhenius expressions
        kt_low = arrhenius(plow_params, t_ref, temps)
        kt_high = arrhenius(phigh_params, t_ref, temps)

        # Calculate K(T,P)s with PLOG expression
        logkt = (
            np.log10(kt_low) +
            ((np.log10(kt_high) - np.log10(kt_low)) * pres_term)
        )
        ktps = 10**(logkt)

    return ktps


def chebyshev(alpha, tmin, tmax, pmin, pmax, temps, pressures):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Chebyshev functional expression.

        :param alpha: Chebyshev coefficient matrix
        :type alpha: numpy.ndarray
        :param tmin: minimum temperature Chebyshev model is defined
        :type tmin: float
        :param tmax: maximum temperature Chebyshev model is defined
        :type tmax: float
        :param pmin: minimum pressure Chebyshev model is defined
        :type pmin: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressures: Pressures used to calculate k(T,P)s
        :type pressures: list(float)
        :return ktp_dct: k(T,Ps) at all temps and pressures
        :rtype: dict[pressure: temps]
    """
    ktp_dct = {}
    for pressure in pressures:
        ktp_dct[pressure] = chebyshev_one_pressure(
            alpha, tmin, tmax, pmin, pmax, temps, pressure)

    return ktp_dct


def chebyshev_one_pressure(alpha, tmin, tmax, pmin, pmax, temps, pressure):
    """ Calculates T,P-dependent rate constants [k(T,P)]s using
        a Chebyshev functional expression, at a given pressure,
        across several temperatures.

        :param alpha: Chebyshev coefficient matrix
        :type alpha: numpy.ndarray
        :param tmin: minimum temperature Chebyshev model is defined
        :type tmin: float
        :param tmax: maximum temperature Chebyshev model is defined
        :type tmax: float
        :param pmin: minimum pressure Chebyshev model is defined
        :type pmin: float
        :param pmax: maximum pressure Chebyshev model is defined
        :type pmax: float
        :param temps: Temps used to calculate high- and low-k(T)s
        :type temps: numpy.ndarray
        :param pressure: Pressure used to calculate k(T,P)s
        :type pressure: float
        :return ktps: Set of k(T,P)s at given pressure
        :rtype numpy.ndarray
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


def _pr_term(highp_rateks, lowp_rateks, temps, pressure,
             collid_factor=1.0, rval=RC2):
    """ Calculates the reduced pressure term for a single pressure
        used for Lindemann and Troe P-dependent functional expressions.

        :param list highp_ks: k(T)s determined at high-pressure
        :type highp_ks: numpy.ndarray
        :param list lowp_ks: k(T)s determined at low-pressure
        :type lowp_ks: numpy.ndarray
        :param temps: Temps used to calculate high- and low-k(T)s
        :temps: numpy.ndarray
        :param pressure: Pressure used to calculate reduced pressure
        :type pressure: float
        :rtype: numpy.ndarray
    """

    pr_term = (
        (lowp_rateks / highp_rateks) *
        p_to_m(pressure, temps, rval=rval) *
        collid_factor
    )

    return pr_term


def _f_broadening_term(pr_term, alpha, ts3, ts1, ts2, temp):
    """ Calculates the F broadening factor term used for Troe expressions.

        :param pr_term: reduced pressure
        :type pr_term: float
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param temp: temperature
        :type temp: float
        :return f_cent: F broadening term
        :rtype: float
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


def p_to_m(pressure, temps, rval=RC2):
    """ Convert the pressure to the concentration of a gas [M]
        assuming an ideal gas form where [M] ~ P/RT.

        :param pressure: Pressure of gas (in atm)
        :type pressure: float
        :return mconc: Conncentration of Gas (mol/cm^3)
        :rtype: float
    """
    # print(pressure, type(pressure))
    # print(temps, type(temps))
    # print(rval, type(rval))
    return pressure / (rval * temps)
