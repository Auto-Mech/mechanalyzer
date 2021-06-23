"""
  Fit the rate constants read from the MESS output to
  Arrhenius, Plog, Troe, and Chebyshev expressions
"""

import ioformat
from ratefit.fit import arrhenius as arrfit
from ratefit.fit import chebyshev as chebfit
from ratefit.fit import troe as troefit
from ratefit.fit._read import gen_reaction_pairs
from ratefit.fit._read import read_rates


DEFAULT_PDEP_DCT = {
    'temps': (500.0, 1000.0),
    'tol': 20.0,
    'pval': 1.0,
    'plow': None,
    'phigh': None
}
DEFAULT_ARRFIT_DCT = {
    'dbltol': 15.0,
    'dblcheck': 'max'
}
DEFAULT_TROE_DCT = {
    'params': ('ts1', 'ts2', 'ts3', 'alpha'),
    'tol': 20.0
}
DEFAULT_CHEB_DCT = {
    'tdeg': 6,
    'pdeg': 4,
    'tol': 20.0
}


def fit_ktp_dct(mess_path, inp_fit_method,
                pdep_dct=None,
                arrfit_dct=None,
                chebfit_dct=None,
                troefit_dct=None,
                label_dct=None,
                fit_temps=None, fit_pressures=None,
                fit_tunit='K', fit_punit='atm'):
    """ Parse the MESS output and fit the rates to
        Arrhenius expressions written as CHEMKIN strings
    """

    # Read the mess input and output strings using the path
    mess_out_str = ioformat.pathtools.read_file(mess_path, 'rate.out')

    # Set dictionaries if they are unprovided
    pdep_dct = pdep_dct or DEFAULT_PDEP_DCT
    arrfit_dct = arrfit_dct or DEFAULT_ARRFIT_DCT
    chebfit_dct = chebfit_dct or DEFAULT_CHEB_DCT
    troefit_dct = troefit_dct or DEFAULT_TROE_DCT

    rxn_pairs, label_dct = gen_reaction_pairs(mess_out_str, label_dct)

    # Loop through reactions, fit rates, and write ckin strings
    chemkin_str_dct = {}
    for (lab_i, lab_j) in rxn_pairs:

        # Set the name and A conversion factor
        reaction = label_dct[lab_i] + '=' + label_dct[lab_j]
        print('------------------------------------------------\n')
        print('Reading and Fitting Rates for {}'.format(reaction))

        # Read the rate constants out of the mess outputs
        print('\nReading k(T,P)s from MESS output...')
        ktp_dct, cheb_fit_temps = read_rates(
            mess_out_str, pdep_dct, lab_i, lab_j,
            fit_temps=fit_temps, fit_pressures=fit_pressures,
            fit_tunit=fit_tunit, fit_punit=fit_punit)

        # Check the ktp dct and fit_method to see how to fit rates
        fit_method = _assess_fit_method(ktp_dct, inp_fit_method)

        # Get the desired fits in the form of CHEMKIN strs
        if fit_method is None:
            continue
        if fit_method == 'arrhenius':
            chemkin_str = arrfit.pes(
                ktp_dct, reaction, mess_path, **arrfit_dct)
        elif fit_method == 'chebyshev':
            chemkin_str = chebfit.pes(
                ktp_dct, reaction, mess_path,
                fit_temps=cheb_fit_temps, **chebfit_dct)
            if not chemkin_str:
                chemkin_str = arrfit.pes(
                    ktp_dct, reaction, mess_path, **arrfit_dct)
        elif fit_method == 'troe':
            chemkin_str += troefit.pes(
                ktp_dct, reaction, mess_path, **troefit_dct)

        # Update the chemkin string dct {PES FORMULA: [PES CKIN STRS]}
        print('\nFinal Fitting Parameters in CHEMKIN Format:', chemkin_str)
        ridx = reaction.replace('=', '_')
        chemkin_str_dct.update({ridx: chemkin_str})

    return chemkin_str_dct


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
            fit_method = inp_fit_method
    else:
        fit_method = None

    # Print message to say what fitting will be done
    if fit_method == 'arrhenius':
        if inp_fit_method != 'arrhenius':
            print(
                '\nRates at not enough pressures for Troe/Chebyshev.')
        print(
            '\nFitting k(T,P)s to PLOG/Arrhenius Form....')
    elif fit_method == 'chebyshev':
        print(
            '\nFitting k(T,P)s to Chebyshev Form...')
    elif fit_method == 'troe':
        print(
            '\nFitting k(T,P)s to Tree Form...')
    elif fit_method is None:
        print(
            '\nNo valid k(T,Ps)s from MESS output to fit.',
            'Skipping to next reaction...')

    return fit_method
