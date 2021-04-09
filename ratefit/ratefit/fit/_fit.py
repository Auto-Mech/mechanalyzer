"""
  Fit the rate constants read from the MESS output to
  Arrhenius, Plog, Troe, and Chebyshev expressions
"""

import copy
import mess_io
import ioformat
from ratefit.fit import arrhenius as arrfit
from ratefit.fit import chebyshev as chebfit
# from ratefit.fit import troe as troefit
from ratefit.fit import filter_ktp_dct
from ratefit.fit import pressure_dependent_ktp_dct


DEFAULT_PDEP_DCT = {
    'pdep_temps': (500, 100),
    'pdep_tol': 20.0,
    'no_pdep_pval': 1.0,
    'pdep_low': None,
    'pdep_high': None
}
DEFAULT_ARRFIT_DCT = {
    'dblarr_tolerance': 15.0,
    'dblarr_check': 'max'
}
DEFAULT_TROE_DCT = {
    'param_list': ('ts1', 'ts2', 'ts3', 'alpha')
}
DEFAULT_CHEB_DCT = {
    'tdeg': 6,
    'pdeg': 4,
    'fit_tolerance': 20.0
}


def fit_ktp_dct(mess_path, pes_formula, fit_method,
                pdep_dct=DEFAULT_PDEP_DCT,
                arrfit_dct=DEFAULT_ARRFIT_DCT,
                chebfit_dct=DEFAULT_CHEB_DCT,
                troefit_dct=DEFAULT_TROE_DCT,
                label_dct=None,
                fit_temps=None, fit_pressures=None,
                fit_tunit='K', fit_punit='atm'):
    """ Parse the MESS output and fit the rates to
        Arrhenius expressions written as CHEMKIN strings
    """

    # Read the mess input and output strings using the path
    mess_out_str = ioformat.pathtools.read_file(mess_path, 'rate.out')

    # Read the label dct from the MESS file if unprovided
    if label_dct is None:
        labels = mess_io.reader.rates.labels(mess_out_str, read_fake=False)
        label_dct = dict(zip(labels, labels))
    rxn_pairs = gen_reaction_pairs(label_dct)

    # Loop through reactions, fit rates, and write ckin strings
    chemkin_str_dct = {}
    for (name_i, lab_i), (name_j, lab_j) in rxn_pairs:

        # Set the name and A conversion factor
        reaction = name_i + '=' + name_j
        print('Reading and Fitting Rates for {}'.format(reaction))

        # Read the rate constants out of the mess outputs
        print('\nReading k(T,P)s from MESS output...')
        ktp_dct = read_rates(
            mess_out_str, pdep_dct, lab_i, lab_j,
            fit_temps=fit_temps, fit_pressures=fit_pressures,
            fit_tunit=fit_tunit, fit_punit=fit_punit)

        # Check the ktp dct and fit_method to see how to fit rates
        fit_method = _assess_fit_method(ktp_dct, fit_method)

        # Get the desired fits in the form of CHEMKIN strs
        if fit_method is None:
            continue
        if fit_method == 'arrhenius':
            chemkin_str = arrfit.pes(
                ktp_dct, reaction, mess_path, **arrfit_dct)
        elif fit_method == 'chebyshev':
            chemkin_str = chebfit.pes(
                ktp_dct, reaction, mess_path, **chebfit_dct)
            # ktp_dct, inp_temps, reaction, mess_path)
            if not chemkin_str:
                chemkin_str = arrfit.pes(
                    ktp_dct, reaction, mess_path, **arrfit_dct)
        # elif fit_method == 'troe':
        #     chemkin_str += troefit.pes(
        #         ktp_dct, reaction, mess_path, **troefit_dct)

        # Update the chemkin string dct
        print('Final Fitting Parameters in CHEMKIN Format:', chemkin_str)
        ridx = pes_formula + '_' + reaction.replace('=', '_')
        chemkin_str_dct.update({ridx: chemkin_str})

    return chemkin_str_dct


def gen_reaction_pairs(label_dct):
    """ Generate pairs of reactions
    """

    rxn_pairs = ()
    for name_i, lab_i in label_dct.items():
        if 'F' not in lab_i and 'B' not in lab_i:
            for name_j, lab_j in label_dct.items():
                if 'F' not in lab_j and 'B' not in lab_j and lab_i != lab_j:
                    rxn_pairs += (((name_i, lab_i), (name_j, lab_j)),)

    # Only grab the forward reactions, remove the reverse reactions
    sorted_rxn_pairs = ()
    for pair in rxn_pairs:
        rct, prd = pair
        if (rct, prd) in sorted_rxn_pairs or (prd, rct) in sorted_rxn_pairs:
            continue
        sorted_rxn_pairs += ((rct, prd),)

    return sorted_rxn_pairs


# Readers
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
    assert fit_temps <= mess_temps
    assert fit_pressures <= mess_press

    # Read all k(T,P) values from MESS output; filter negative/undefined values
    calc_ktp_dct = mess_io.reader.rates.ktp_dct(
        mess_out_str, rct_lab, prd_lab)

    print(
        'Removing invalid k(T,P)s from MESS output that are either:\n',
        '  (1) negative, (2) undefined [***], or (3) below 10**(-21) if',
        'reaction is bimolecular', newline=1)
    filt_ktp_dct = filter_ktp_dct(calc_ktp_dct, bimol)

    # Filter the ktp dictionary by assessing the presure dependence
    if filt_ktp_dct:
        if list(filt_ktp_dct.keys()) == ['high']:
            print('\nValid k(T)s only found at High Pressure...')
            ktp_dct['high'] = filt_ktp_dct['high']
        else:
            if pdep_dct:
                print(
                    'User requested to assess pressure dependence',
                    'of reaction.', newline=1)
                ktp_dct = pressure_dependent_ktp_dct(
                    filt_ktp_dct,
                    tolderance=pdep_dct['pdep_tol'],
                    pdep_temps=pdep_dct['pdep_temps'],
                    plow=pdep_dct['plow'],
                    phigh=pdep_dct['phigh'],
                    no_pdep_pval=DEFAULT_PDEP_DCT['no_pdep_pval'])
            else:
                ktp_dct = copy.deepcopy(filt_ktp_dct)

    return ktp_dct


def _assess_fit_method(ktp_dct, inp_fit_method):
    """ Assess if there are any rates to fit and if so, check if
        the input fit method should be used, or just simple Arrhenius
        fits will suffice because there is only one pressure for which
        rates exist to be fit.

        # If only one pressure (outside HighP limit), just run Arrhenius
    """

    if ktp_dct:
        pressures = list(ktp_dct.keys())
        npressures = len(pressures)
        if npressures == 1 or (npressures == 2 and 'high' in pressures):
            fit_method = 'arrhenius'
    else:
        fit_method = None

    # Print message to say what fitting will be done
    if fit_method == 'arrhenius':
        if inp_fit_method != 'arrhenius':
            print(
                'Rates at not enough pressures for Troe/Chebyshev.', newline=1)
        print(
            'Fitting k(T,P)s to PLOG/Arrhenius Form....', newline=1)
    elif fit_method == 'chebyshev':
        print(
            'Fitting k(T,P)s to Chebyshev Form...', newline=1)
    elif fit_method == 'troe':
        print(
            'Fitting k(T,P)s to Tree Form...', newline=1)
    elif fit_method is None:
        print(
            'No valid k(T,Ps)s from MESS output to fit.',
            'Skipping to next reaction...', newline=1)

    return fit_method
