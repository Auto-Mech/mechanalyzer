"""
Calculate rates with various fitting functions
"""

import copy
import numpy
from phydat import phycon
import ratefit.calc
from autoreact.params import RxnParams  # do I actually need this?

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


    ratefit.ktpdct.check_p_t(pressures, temps)  # enforce formatting rules
    rxn_ktp_dct = {}
    for rxn, param_tups in rxn_param_dct.items():
        ktp_dct = {}
        for param_tup in param_tups:
            new_ktp_dct = eval_param_tup(param_tup, pressures, temps)
            ktp_dct = add_ktp_dcts(ktp_dct, new_ktp_dct)
        rxn_ktp_dct[rxn] = ktp_dct

    return rxn_ktp_dct


def eval_params(params, pressures, temps, t_ref=1.0):
    """ Look through a params and evaluate k(T,P) based on the contents.
        Return a ktp_dct.

        :param params:
        :type params: autochem/autoreact RxnParams object
        :param pressures:
        :type pressures:
        :param temps:
        :type temps:
        :return ktp_dct: rate constant as a function of temp and pressure
        :rtype: dct {pressure1: {temps1, kts1}, pressures2: ...}
    """
    # Get the list of existing forms
    rxn_forms = params.existing_parameters()
    assert rxn_forms != (), 'The params object is empty'

    # Loop over each rxn type (usually only one, but trying to be general)
    ktp_dct = {}
    for rxn_form in rxn_forms:
        new_ktp_dct = {}
        if rxn_form == 'arr':
            arr_tuples = params.arr
            kts = ratefit.calc.arrhenius(arr_tuples, temps, t_ref)
            new_ktp_dct['high'] = (temps, kts)

        if rxn_form == 'plog':
            plog_dct = params.plog
            new_ktp_dct = ratefit.calc.plog(plog_dct, t_ref, temps, pressures)

        if rxn_form == 'cheb':
            cheb_dct = params.cheb
            alpha = cheb_dct['alpha']
            tlim = cheb_dct['tlim']
            plim = cheb_dct['plim']
            new_ktp_dct = ratefit.calc.cheb(alpha, tlim, plim, temps, pressures)

        # Allow for summations of ktp_dcts in case of more than one rxn_form
        ktp_dct = add_ktp_dcts(new_ktp_dct, ktp_dct)                

# DON'T DELETE THIS
#    elif params[2] is not None:  # Troe
#        assert params[0] is not None, (
#            'Troe parameters are included,',
#            'but the high-P parameters are absent'
#            )
#        assert params[1] is not None, (
#            'Troe and high-P parameters are included,',
#            'but the low-P parameters are absent'
#            )
#        troe_params = params[2]
#        alpha = troe_params[0]
#        ts3 = troe_params[1]  # T***
#        ts1 = troe_params[2]  # T*
#        if len(troe_params) == 4:
#            ts2 = troe_params[3]  # T**; this one is commonly omitted
#        else:
#            ts2 = None
#
#        highp_kts = ratefit.calc.arrhenius(params[0], t_ref, temps)
#        lowp_kts = ratefit.calc.arrhenius(params[1], t_ref, temps)
#        ktp_dct = ratefit.calc.troe(
#            highp_kts, lowp_kts, temps, pressures,
#            alpha, ts3, ts1, ts2, collid_factor=1.0)
#
#    elif params[1] is not None:  # Lindemann
#        assert params[0] is not None, (
#            'Low-P parameters are included,',
#            'but the high-P parameters are absent'
#            )
#        highp_kts = ratefit.calc.arrhenius(params[0], t_ref, temps)
#        lowp_kts = ratefit.calc.arrhenius(params[1], t_ref, temps)
#        ktp_dct = ratefit.calc.lindemann(
#            highp_kts, lowp_kts,
#            temps, pressures, collid_factor=1.0)
#
#    else:  # Arrhenius
#        assert params[0] is not None, (
#            'The params does not seem to contain any useful information.'
#            )
#        # Case is unique as kTP dict contains only one key-value pair
#        kts = ratefit.calc.arrhenius(params[0], t_ref, temps)
#        ktp_dct = {}
#        ktp_dct['high'] = (temps, kts)

    return ktp_dct


def add_ktp_dcts(ktp_dct1, ktp_dct2):
    """ Add the rates in two ktp_dcts.

        The input dcts should have identical P and T values. May not always
        be true if one is from PLOG and one is from an expression
        with a pressure-independent rate.
    """

    # If either starting dct is empty, simply copy the other
    if ktp_dct1 == {}: 
        added_dct = copy.deepcopy(ktp_dct2)
    elif ktp_dct2 == {}: 
        added_dct = copy.deepcopy(ktp_dct1)
    # Otherwise, add the dcts
    else:
        added_dct = {}
        # Loop over the pressures in the first dct
        for pressure, (temps, kts1) in ktp_dct1.items():
            if pressure in ktp_dct2.keys():
                (_, kts2) = ktp_dct2[pressure]  # unpack ktp_dct2 values
            else:  # if the pressure is not in ktp_dct2
                kts2 = numpy.zeros(numpy.shape(kts1))
            added_kts = kts1 + kts2
            added_dct[pressure] = (temps, added_kts)  # store added values
        # Loop over the pressures in the second dct and see if any are new
        for pressure, (temps, kts2) in ktp_dct2.items():
            # If the pressure is not in added_dct, it is unique to ktp_dct
            if pressure not in added_dct.keys():
                (_, kts2) = ktp_dct2[pressure]  # unpack ktp_dct2 values
                added_dct[pressure] = (temps, kts2)  # store new values

    return added_dct
