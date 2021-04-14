"""
Functions to handle various aspects of pressure dependence
"""

import copy
import numpy as np


INI_PDEP_DCT = {
    'pdep_temps': (500, 100),
    'pdep_tol': 20.0,
    'no_pdep_pval': 1.0,
    'plow': None,
    'phigh': None
}


def pressure_dependent_ktp_dct(inp_ktp_dct,
                               tolerance=INI_PDEP_DCT['pdep_tol'],
                               pdep_temps=INI_PDEP_DCT['pdep_temps'],
                               plow=INI_PDEP_DCT['plow'],
                               phigh=INI_PDEP_DCT['phigh'],
                               no_pdep_pval=INI_PDEP_DCT['no_pdep_pval']):
    """ Takes a full ktp dictionary, assesses if there is significant
        pressure dependnce in the rate constants.

        If so, return the pressure dependent rates.
        If not, return rates at a single pressure equal to some valeu

        :param ktp_dct:
        :type ktp_dct:
    """

    # Assess the pressure dependence of the rate constants
    rxn_is_pdependent = assess_pressure_dependence(
        inp_ktp_dct,
        tolerance=tolerance,
        assess_pdep_temps=pdep_temps,
        plow=plow,
        phigh=phigh)

    # Build the rate constants

    # no pdep_dct is amde if no pdependence found and no_pdep_pval not in filtered ktp dct
    # no high may be in there because rates may be undefined
    # pdependecne check could fail if assess temp ranges not big enough
    # pdep check lowest/highest temp

    if rxn_is_pdependent:
        print('Reaction found to be pressure dependent.',
              'Fitting all k(T)s from all pressures',
              'found in MESS.')
        pdep_ktp_dct = copy.deepcopy(inp_ktp_dct)
    else:
        print('No pressure dependence detected.',
              'Grabbing k(T)s at {} atm'.format(no_pdep_pval))
        # print('pval', no_pdep_pval)
        # print('ktpdct\n', inp_ktp_dct)
        if no_pdep_pval in inp_ktp_dct:
            pdep_ktp_dct = {'high': inp_ktp_dct[no_pdep_pval]}
        else:
            pdep_ktp_dct = None

    return pdep_ktp_dct


def assess_pressure_dependence(tk_dct, assess_pdep_temps,
                               tolerance=20.0, plow=None, phigh=None):
    """ Assess if there are significant differences between k(T,P) values
        at low-pressure and high-pressure, signaling that a given reaction is
        pressure dependent. Function assesses these changes across all
        temperatures at both pressures.

        :param tk_dct: T,k pairs
        :type tk_dict: dict[pressure] = [temps, k(T, P)s]
        :param assess_pdep_temps: Temperatures to assess P-dependence
        :type assess_pdep_temps: list(float)
        :param tolerance: %-difference cutoff for changes in k(T,P)
        :type tolerance: float
        :param plow: Minimum pressure used to assess P-dependence
        :type plow: float
        :param phigh: Maximum pressure used to assess P-dependence
        :type phigh: float
        :return is_pressure_dependent: variable signaling P-dependence
        :rtype: bool
    """
    # Get list of the sorted pressures, ignoring the high-pressure limit rates
    pressures = [pressure for pressure in tk_dct
                 if pressure != 'high']
    pressures.sort()

    # Set the low- and high-pressure if not specified by user
    if plow is None:
        plow = min(pressures)
    if phigh is None:
        phigh = max(pressures)

    # Check % difference for k(T, P) vals
    is_pressure_dependent = False
    if plow in tk_dct and phigh in tk_dct:

        # Loop over temps to examine for large % dif in k(T) at low- and high-P
        for temp_compare in assess_pdep_temps:
            # For low- and high-P, find the idx for the temp in temp_compare
            temps_low = tk_dct[plow][0]
            temps_high = tk_dct[phigh][0]
            temp_low_match = np.where(np.isclose(temps_low, temp_compare))[0]
            temp_high_match = np.where(np.isclose(temps_high, temp_compare))[0]
            if temp_low_match.size > 0 and temp_high_match.size > 0:
                temp_low_idx = temp_low_match[0]
                temp_high_idx = temp_high_match[0]
                # Grab the k(T, P) vale for the approprite temp and pressure
                ktp_low = tk_dct[plow][1][temp_low_idx]
                ktp_high = tk_dct[phigh][1][temp_high_idx]
                # Calculate the % difference and see if above threshold
                ktp_dif = (abs(ktp_low - ktp_high) / ktp_low) * 100.0
                if ktp_dif > tolerance:
                    is_pressure_dependent = True

    return is_pressure_dependent
