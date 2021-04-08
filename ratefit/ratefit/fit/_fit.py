"""
  Fit the rate constants read from the MESS output to
  Arrhenius, Plog, Troe, and Chebyshev expressions
"""

import copy
import numpy
import ratefit
import mess_io
from phydat import phycon
from mechlib.amech_io import writer
from mechlib.amech_io import printer as print
from mechroutines.pf.ktp.fit import _arr as arr
from mechroutines.pf.ktp.fit import _cheb as cheb


INI_ARRFIT_DCT = {}
INI_TROE_DCT = {}


def fit_ktp_dct(mess_path, pes_formula,
                pdep_dct={},
                arrfit_dct={}, troe_dct={},
                label_dct=None,
                fit_temps=None, fit_pressure=None
                fit_tunit='K', fit_punit='atm'):
    """ Parse the MESS output and fit the rates to
        Arrhenius expressions written as CHEMKIN strings
    """

    # Read the mess input and output strings using the path
    mess_inp_str, mess_out_str = _read_mess(mess_path)

    # Read the label dct from the MESS file if unprovided
    if label_dct is None:
        label_dct = _mess_labels(mess_path)

    # Loop through reactions, fit rates, and write ckin strings
    rxn_pairs = gen_reaction_pairs(label_dct)
    for (name_i, lab_i), (name_j, lab_j) in rxn_pairs:

        # Set the name and A conversion factor
        reaction = name_i + '=' + name_j
        a_conv_factor = phycon.NAVO if 'W' not in lab_i else 1.00
        bimol = False if 'W' not in lab_i else True

        print('Reading and Fitting Rates for {}'.format(reaction))

        # Read the rate constants out of the mess outputs
        print('\nReading k(T,P)s from MESS output...')
        ktp_dct = read_rates(
            inp_temps, inp_pressures, inp_tunit, inp_punit,
            lab_i, lab_j, mess_path, pdep_fit,
            bimol=bimol)

        # Check the ktp dct and fit_method to see how to fit rates
        fit_method = _assess_fit_method(ktp_dct, inp_fit_method)

        # Get the desired fits in the form of CHEMKIN strs
        if fit_method is None:
            continue
        if fit_method == 'arrhenius':
            chemkin_str = arr.perform_fits(
                ktp_dct, reaction, mess_path,
                a_conv_factor, arrfit_thresh)
        elif fit_method == 'chebyshev':
            chemkin_str = cheb.perform_fits(
                ktp_dct, inp_temps, reaction, mess_path,
                a_conv_factor)
            if not chemkin_str:
                chemkin_str = arr.perform_fits(
                    ktp_dct, reaction, mess_path,
                    a_conv_factor, arrfit_thresh)
        # elif fit_method == 'troe':
        #     # chemkin_str += troe.perform_fits(
        #     #     ktp_dct, reaction, mess_path,
        #     #     troe_param_fit_lst,
        #     #     a_conv_factor, err_thresh)

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
def read_rates(mess_out_str, pdep_fit, rct_lab, prd_lab,
               bimol=False,
               fit_temps=None, fit_pressures=None
               fit_tunit=None, fit_punit=None):
    """ Read the rate constants from the MESS output and
        (1) filter out the invalid rates that are negative or undefined
        and obtain the pressure dependent values
    """

    # Read temperatures, pressures and rateks from MESS output
    mess_temps, tunit = mess_io.reader.rates.temperatures(output_str)
    mess_pressures, punit = mess_io.reader.rates.pressures(output_str)
    
    fit_temps = fit_temps if fit_temps is not None else mess_temps
    fit_pressures = fit_pressures if fit_pressures is not None else mess_pressures
    fit_tunit = fit_tunit if fit_tunit is not None else mess_tunit
    fit_punit = fit_punit if fit_punit is not None else mess_punit

    fit_temps = list(set(list(fit_temps)))
    fit_temps.sort()
    assert inp_temps <= mess_temps
    assert inp_pressures <= mess_pressures

    # Loop over the pressures obtained from the MESS output
    calc_ktp_dct = mess_io.reader.rates.ktp_dct(
        output_string, rct_lab, prd_lab)

    # Remove k(T) vals at each P where where k is negative or undefined
    # If ANY valid k(T,P) vals at given pressure, store in dct
    print(
        'Removing invalid k(T,P)s from MESS output that are either:\n',
        '  (1) negative, (2) undefined [***], or (3) below 10**(-21) if',
        'reaction is bimolecular', newline=1)
    filt_ktp_dct = ratefit.fit.filter_ktp_dct(calc_ktp_dct, bimol)

    # Filter the ktp dictionary by assessing the presure dependence
    if filt_ktp_dct:
        if list(filt_ktp_dct.keys()) == ['high']:
            print(
                'Valid k(T)s only found at High Pressure...', newline=1)
            ktp_dct['high'] = filt_ktp_dct['high']
        else:
            # if pdep_dct:
            if pdep_fit:
                print(
                    'User requested to assess pressure dependence',
                    'of reaction.', newline=1)
                rxn_is_pdependent = ratefit.fit.assess_pressure_dependence(
                    filt_ktp_dct,
                    tolderance=pdep_dct['pdep_tol'],
                    pdep_temps=pdep_dct['pdep_temps'],
                    plow=pdep_dct['plow'], plow=pdep_dct['plow'])
                if rxn_is_pdependent:
                    print(
                        'Reaction found to be pressure dependent.',
                        'Fitting all k(T)s from all pressures',
                        'found in MESS.', indent=1/2.)
                    ktp_dct = copy.deepcopy(filt_ktp_dct)
                else:
                    no_pdep_pval = pdep_fit['no_pdep_pval']
                    print(
                        'No pressure dependence detected.',
                        'Grabbing k(T)s at {} {}'.format(
                            no_pdep_pval, punit), newline=1)
                    if no_pdep_pval in filt_ktp_dct:
                        ktp_dct['high'] = filt_ktp_dct[no_pdep_pval]
            else:
                ktp_dct = copy.deepcopy(filt_ktp_dct)

        # Patchy way to get high-pressure rates in dct if needed
        if 'high' not in ktp_dct and 'high' in filt_ktp_dct.keys():
            ktp_dct['high'] = filt_ktp_dct['high']
    else:
        ktp_dct = None

    return ktp_dct


def _read_mess_path(mess_path):
    """ Read the input and output strings from the MESS path
    """

    

def read_labels(mess_path):
    """ get the labels from the mess path
    """
    labels = mess_io.reader.labels(mess_inp_str, read_fake=False)
    label_dct = dict(zip(label_dct, label_dct))

    return label_dct


def _assess_fit_method(ktp_dct, fit_method):
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
