""" Get all of the required reactions and rate constants
"""

import copy
import itertools
import mess_io
from ratefit.fit._util import filter_ktp_dct
from ratefit.fit._pdep import pressure_dependent_ktp_dct


def gen_reaction_pairs(mess_out_str, label_dct):
    """ Generate pairs of reactions
    """

    sorted_rxn_pairs = mess_io.reader.rates.reactions(
        mess_out_str,
        read_rev=True,
        read_fake=False,
        read_self=False)

    if label_dct is None:
        labels = set(itertools.chain(*sorted_rxn_pairs))
        label_dct = dict(zip(labels, labels))

    return sorted_rxn_pairs, label_dct


def read_rates(mess_out_str, pdep_dct, rct_lab, prd_lab,
               fit_temps=None, fit_pressures=None,
               fit_tunit=None, fit_punit=None):
    """ Read the rate constants from the MESS output and
        (1) filter out the invalid rates that are negative or undefined
        and obtain the pressure dependent values
    """

    # Initialize vars
    ktp_dct = {}
    bimol = bool('W' not in rct_lab)

    # Read temperatures, pressures and rateks from MESS output
    mess_temps, tunit = mess_io.reader.rates.temperatures(mess_out_str)
    mess_press, punit = mess_io.reader.rates.pressures(mess_out_str)

    fit_temps = fit_temps if fit_temps is not None else mess_temps
    fit_pressures = fit_pressures if fit_pressures is not None else mess_press
    fit_tunit = fit_tunit if fit_tunit is not None else tunit
    fit_punit = fit_punit if fit_punit is not None else punit

    fit_temps = list(set(list(fit_temps)))
    fit_temps.sort()
    assert set(fit_temps) <= set(mess_temps)
    assert set(fit_pressures) <= set(mess_press)

    # Read all k(T,P) values from MESS output; filter negative/undefined values
    calc_ktp_dct = mess_io.reader.rates.ktp_dct(
        mess_out_str, rct_lab, prd_lab)
    print(
        '\nRemoving invalid k(T,P)s from MESS output that are either:\n',
        '  (1) negative, (2) undefined [***], or (3) below 10**(-21) if',
        'reaction is bimolecular')
    filt_ktp_dct = filter_ktp_dct(calc_ktp_dct, bimol)

    # Filter the ktp dictionary by assessing the presure dependence
    if filt_ktp_dct:
        if list(filt_ktp_dct.keys()) == ['high']:
            print('\nValid k(T)s only found at High Pressure...')
            ktp_dct['high'] = filt_ktp_dct['high']
        else:
            if pdep_dct:
                print(
                    '\nUser requested to assess pressure dependence',
                    'of reaction.')
                ktp_dct = pressure_dependent_ktp_dct(filt_ktp_dct, **pdep_dct)
            else:
                ktp_dct = copy.deepcopy(filt_ktp_dct)

    return ktp_dct, fit_temps
