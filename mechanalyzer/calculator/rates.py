"""
Calculate rates with various fitting functions
"""

import copy
import numpy
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
        """ Add the rates in two ktp_dcts.

            The input dcts should have identical P and T values. May not always
            be true if one is from PLOG and one is from an expression
            with a pressure-independent rate.
        """

        if not ktp_dct1:  # if starting dct is empty, copy the addition
            added_dct = copy.deepcopy(ktp_dct2)
        else:
            added_dct = {}
            for pressure, (temps, kts1) in ktp_dct1.items():
                if pressure in ktp_dct2.keys():
                    (_, kts2) = ktp_dct2[pressure]  # unpack ktp_dct2 values
                else:  # if the pressure is not in ktp_dct2
                    kts2 = numpy.zeros(numpy.shape(kts1))
                added_kts = kts1 + kts2
                added_dct[pressure] = (temps, added_kts)  # store added values
        return added_dct

    ratefit.ktpdct.check_p_t(pressures, temps)  # enforce formatting rules
    rxn_ktp_dct = {}
    for rxn, param_tups in rxn_param_dct.items():
        ktp_dct = {}
        for param_tup in param_tups:
            new_ktp_dct = eval_param_tup(param_tup, pressures, temps)
            ktp_dct = add_ktp_dcts(ktp_dct, new_ktp_dct)
        rxn_ktp_dct[rxn] = ktp_dct

    return rxn_ktp_dct


def eval_param_tup(param_tup, pressures, temps, t_ref=1.0):
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
        ktp_dct = ratefit.calc.plog(plog_dct, t_ref, temps, pressures)

    elif param_tup[2] is not None:  # Troe
        assert param_tup[0] is not None, (
            'Troe parameters are included,',
            'but the high-P parameters are absent'
            )
        assert param_tup[1] is not None, (
            'Troe and high-P parameters are included,',
            'but the low-P parameters are absent'
            )
        troe_params = param_tup[2]
        alpha = troe_params[0]
        ts3 = troe_params[1]  # T***
        ts1 = troe_params[2]  # T*
        if len(troe_params) == 4:
            ts2 = troe_params[3]  # T**; this one is commonly omitted
        else:
            ts2 = None

        highp_kts = ratefit.calc.arrhenius((param_tup[0],), t_ref, temps)
        lowp_kts = ratefit.calc.arrhenius((param_tup[1],), t_ref, temps)
        ktp_dct = ratefit.calc.troe(
            highp_kts, lowp_kts, temps, pressures,
            alpha, ts3, ts1, ts2, collid_factor=1.0)

    elif param_tup[1] is not None:  # Lindemann
        assert param_tup[0] is not None, (
            'Low-P parameters are included,',
            'but the high-P parameters are absent'
            )
        highp_kts = ratefit.calc.arrhenius((param_tup[0],), t_ref, temps)
        lowp_kts = ratefit.calc.arrhenius((param_tup[1],), t_ref, temps)
        ktp_dct = ratefit.calc.lindemann(
            highp_kts, lowp_kts,
            temps, pressures, collid_factor=1.0)

    else:  # Arrhenius
        assert param_tup[0] is not None, (
            'The param_tup does not seem to contain any useful information.'
            )
        # Case is unique as kTP dict contains only one key-value pair
        highp_params = param_tup[0]

        kts = ratefit.calc.arrhenius((highp_params,), t_ref, temps)
        ktp_dct = {}
        ktp_dct['high'] = (temps, kts)

    return ktp_dct
