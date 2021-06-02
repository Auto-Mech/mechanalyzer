"""
Calculate rates with various fitting functions
"""

import copy
import numpy as np
from scipy.special import eval_chebyt
from phydat import phycon
import ratefit.calc


RC = phycon.RC_CAL  # gas constant in cal/(mol.K)
RC2 = phycon.RC_ATM  # gas constant in cm^3.atm/(mol.K)


def eval_rxn_param_dct(rxn_param_dct, pressures, temps):
    """ Loop through all rxns in a rxn_param_dct and get a ktp_dct for
        each one. Return a rxn_ktp_dct.

        :param rxn_param_dct:
        :type rxn_param_dct:
        :param pressures:
        :type pressures:
        :param temps:
        :type temps:
    """
    def add_ktp_dcts(ktp_dct1, ktp_dct2):
        """ Add the rates in two ktp_dcts. The input dcts should have identical P and T values.
            However, this may not always be true if one is from PLOG and one is from an expression
            with a pressure-independent rate.
        """
        if ktp_dct1 == {}:  # if the starting dct is empty, just copy the addition
            added_dct = copy.deepcopy(ktp_dct2)
        else:
            added_dct = {}
            for pressure, (temps, kts1) in ktp_dct1.items():
                if pressure in ktp_dct2.keys():
                    (_, kts2) = ktp_dct2[pressure]  # unpack ktp_dct2 values
                else:  # if the pressure is not in ktp_dct2
                    kts2 = np.zeros(np.shape(kts1))  # use an array of zeros
                added_kts = kts1 + kts2
                added_dct[pressure] = (temps, added_kts)  # store added values

        return added_dct

    check_p_t(pressures, temps)  # enforce formatting rules
    rxn_ktp_dct = {}
    for rxn, param_tups in rxn_param_dct.items():
        ktp_dct = {}
        for param_tup in param_tups:
            new_ktp_dct = eval_param_tup(param_tup, pressures, temps)
            ktp_dct = add_ktp_dcts(ktp_dct, new_ktp_dct)
        rxn_ktp_dct[rxn] = ktp_dct

    return rxn_ktp_dct


def eval_param_tup(param_tup, pressures, temps):
    """ Look through a param_tup and evaluate k(T,P) based on the contents.
        Return a ktp_dct.

    :param param_tup:
    :type param_tup:
    :param pressures:
    :type pressures:
    :param temps:
    :type temps:
    :return ktp_dct: rate constant as a function of temp and pressure
    :rtype: dct {pressure1: {temps1, kts1}, pressures2: ...}
    """
    if param_tup[3] is not None:  # Chebyshev
        alpha = param_tup[3]['alpha_elm']
        t_limits = param_tup[3]['t_limits']
        p_limits = param_tup[3]['p_limits']
        ktp_dct = ratefit.calc.chebyshev(
            alpha, t_limits[0], t_limits[1], p_limits[0], p_limits[1],
            temps, pressures)

    elif param_tup[4] is not None:  # PLOG
        plog_dct = param_tup[4]
        ktp_dct = ratefit.calc.plog(plog_dct, temps, pressures, t_ref=1.0)

    elif param_tup[2] is not None:  # Troe
        assert param_tup[0] is not None, (
            'Troe parameters are included,',
            'but the high-P parameters are absent'
            )
        assert param_tup[1] is not None, (
            'Troe and high-P parameters are included,',
            'but the low-P parameters are absent'
            )
        highp_params = param_tup[0]
        lowp_params = param_tup[1]
        troe_params = param_tup[2]
        alpha = troe_params[0]
        ts3 = troe_params[1]  # T***
        ts1 = troe_params[2]  # T*
        if len(troe_params) == 4:
            ts2 = troe_params[3]  # T**; this one is commonly omitted
        else:
            ts2 = None

        ktp_dct = ratefit.calc.troe(
            highp_params, lowp_params, temps, pressures,
            alpha, ts3, ts1, ts2, collid_factor=1.0)

    elif param_tup[1] is not None:  # Lindemann
        assert param_tup[0] is not None, (
            'Low-P parameters are included,',
            'but the high-P parameters are absent'
            )
        highp_params = param_tup[0]
        lowp_params = param_tup[1]
        ktp_dct = ratefit.calc.lindemann(
            highp_params, lowp_params,
            temps, pressures, collid_factor=1.0)

    else:  # Arrhenius
        assert param_tup[0] is not None, (
            'The param_tup does not seem to contain any useful information.'
            )
        # Case is unique as kTP dict contains only one key-value pair
        highp_params = param_tup[0]
        kts = ratefit.calc.arrhenius(highp_params, temps, t_ref=1.0)
        ktp_dct = {}
        ktp_dct['high'] = (temps, kts)

    return ktp_dct


def ktp(kp_dct, temps, highp_params=None, t_ref=1.0):
    """ Creates a kTP dictionary from a kP dictionary and an
        array of temperatures

        :param kp_dct: dct of the form {P:[k@T1, k@T2]}
        :type kp_dct: {float: np.array}
        :param temps: array of temperatures, either 1- or 2-D
        :type temps: np.ndarray
        :return ktp_dct:
        :rtype:

    """
    ktp_dct = {}
    for index, pressure in enumerate(kp_dct):
        kts = kp_dct[pressure]
        if np.ndim(temps) == 1:  # if the dimensionality of temps is 1
            ktp_dct[pressure] = (temps,kts)
        else:  # if the dimensionality of temps is 2
            ktp_dct[pressure] = (temps[index],kts)

    # Add the high-P k(T)s to the kTP dictionary if needed
    if highp_params:
        if np.ndim(temps) == 1:
            highp_kts = arrhenius(highp_params, temps, t_ref)
            ktp_dct['high'] = (temps,highp_kts)
        else:
            highp_kts = arrhenius(highp_params, temps[-1], t_ref)  # use the last value in temps
            ktp_dct['high'] = (temps[-1],highp_kts)

    return ktp_dct


def check_p_t(pressures, temps):
    """ Enforces certain rules regarding the pressure and temperature arrays.

        :param pressures: array of pressures
        :type pressures: numpy.ndarray
        :param temps: array of temps
        :type temps: numpy.ndarray
    """
    # Check that the dimensionality of the temps array is either 1 or 2
    temp_dim = np.ndim(temps)
    assert temp_dim in (1,2), (
        f'The dimensionality of temps is {temp_dim}; it should be either 1 or 2'
        )

    # If temps is 2-D, enforce that the # of values in each temp array matches the # of pressures
    if temp_dim == 2:
        len_temps = np.shape(temps)[1]
        len_pressures = len(pressures)
        assert len_pressures ==  len_temps,(
            f'# of pressures is {len_pressures}, while # of temps in each array is {len_temps}'
            )
