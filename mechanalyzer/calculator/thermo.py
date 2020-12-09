""" Calculates thermodynamic quantities using thermo data strings
"""


import numpy as np
from chemkin_io.parser import thermo as thm_parser
from lib.phydat import phycon


RC = phycon.RC_cal  # gas constant in cal/(mol.K) 


def create_spc_thermo_dct(spc_nasa7_dct, temps, rval=RC):
    """ Create a spc_thermo_dct. If left with the default input rval=phycon.RC_kcal, the 
        thermo quantities will have units of kcal/mol for h(T) and g(T) and units of kcal/mol-K
        for cp(T) and s(T). 
        
        :param spc_nasa7_dct: dictionary describing the NASA-7 polynomials
        :type spc_nasa7_dct: dct {spc:[list of NASA-7 params]) 
        :param temps: list of temperatures at which thermo is to be evaluated
        :type temps: list
        :param rval: universal gas constant (units decided by the user)
        :type rval: float
        :return spc_thermo_dct: dictionary indexed by species containing arrays of T, h(T), cp(T), s(T), and g(T)
        :rtype: dct of lists {spc:[temps, h_t, cp_t, s_t, g_t]} 
    """
    spc_thermo_dct = {}
    for spc, nasa7_params in spc_nasa7_dct.items():
        h_t, cp_t, s_t, g_t, = [], [], [], [] 
        for temp in temps:
            h_t.append(enthalpy(nasa7_params, temp, rval=rval))
            cp_t.append(heat_capacity(nasa7_params, temp, rval=rval)) 
            s_t.append(entropy(nasa7_params, temp, rval=rval)) 
            g_t.append(gibbs(nasa7_params, temp, rval=rval)) 

            if h_t[-1] is None or cp_t[-1] is None or s_t[-1] is None or g_t[-1] is None:
                print(f'Failure to calculate thermo at {temp} K for the species {spc}, likely due to an invalid temp.')
            
        spc_thermo_dct[spc] = [temps, h_t, cp_t, s_t, g_t]

    return spc_thermo_dct

    
def enthalpy(nasa7_params, temp, rval=RC):
    """ Calculate the enthalpy of a species using the
        coefficients of its NASA-7 polynomial.

        :param nasa7_params: values describing a NASA-7 polynomial
        :type nasa7_params: list
        :param temp: temperature to calculate the enthalpy (K)
        :type temp: float
        :param rval: universal gas constant (units decided by the user)
        :type rval: float
        :return h_t: value for the enthalpy (units same as the input rval)
        :rtype: float
    """
    cfts = coeffs_for_specific_temp(nasa7_params, temp)
    if cfts:
        h_t = (
            cfts[0] +
            ((cfts[1] * temp) / 2.0) +
            ((cfts[2] * temp**2) / 3.0) +
            ((cfts[3] * temp**3) / 4.0) +
            ((cfts[4] * temp**4) / 5.0) +
            (cfts[5] / temp)
        )
        h_t *= (rval * temp)
    else:
        h_t = None

    return h_t


def heat_capacity(nasa7_params, temp, rval=RC):
    """ Calculate the heat capacity of a species using the
        coefficients of its NASA-7 polynomial.

        :param nasa7_params: values describing a NASA-7 polynomial
        :type nasa7_params: list
        :param temp: temperature to calculate heat capacity (K)
        :type temp: float
        :param rval: universal gas constant (units decided by the user)
        :type rval: float
        :return cp_t: value for the heat capacity (units same as the input rval) 
        :rtype: float
    """
    cfts = coeffs_for_specific_temp(nasa7_params, temp)
    if cfts:
        cp_t = (
            cfts[0] +
            (cfts[1] * temp) +
            (cfts[2] * temp**2) +
            (cfts[3] * temp**3) +
            (cfts[4] * temp**4)
        )
        cp_t *= rval
    else:
        cp_t = None

    return cp_t


def entropy(nasa7_params, temp, rval=RC):
    """ Calculate the entropy of a species using the
        coefficients of its NASA-7 polynomial.

        :param nasa7_params: values describing a NASA-7 polynomial
        :type nasa7_params: list
        :param float temp: temperature to calculate the entropy (K)
        :type temp: float
        :param rval: universal gas constant (units decided by the user)
        :type rval: float
        :return s_t: value for the entropy (units same as the input rval)
        :rtype: float
    """
    cfts = coeffs_for_specific_temp(nasa7_params, temp)
    if cfts is not None:
        s_t = (
            (cfts[0] * np.log(temp)) +
            (cfts[1] * temp) +
            ((cfts[2] * temp**2) / 2.0) +
            ((cfts[3] * temp**3) / 3.0) +
            ((cfts[4] * temp**4) / 4.0) +
            (cfts[6])
        )
        s_t *= rval
    else:
        s_t = None

    return s_t


def gibbs(nasa7_params, temp, rval=RC):
    """ Calculate the Gibbs free energy of a species using the
        coefficients of its NASA-7 polynomial.

        :param nasa7_params: values describing a NASA-7 polynomial
        :type nasa7_params: list
        :param temp: temperature to calculate the Gibbs free energy
        :type temp: float
        :param rval: universal gas constant (units decided by the user)
        :type rval: float
        :return g_t: value for the Gibbs free energy (units same as the input rval)
        :rtype: float
    """
    h_t = enthalpy(nasa7_params, temp, rval=rval)
    s_t = entropy(nasa7_params, temp, rval=rval)
    if h_t is not None and s_t is not None:
        g_t = h_t - (s_t * temp)
    else:
        g_t = None

    return g_t


def coeffs_for_specific_temp(nasa7_params, temp):
    """ Gets either the low-T or high-T NASA-7 polynomial coefficients 
        depending on the temperature

    """
    cutoff_temps = nasa7_params[3]
    low_temp, high_temp, mid_temp = cutoff_temps  # the order seems odd, but it's the NASA format

    if low_temp <= temp <= mid_temp:  
        cfts = nasa7_params[4][1]
    elif mid_temp < temp <= high_temp:
        cfts = nasa7_params[4][0]
    else:
        cfts = None
        print(f'The desired temp, {temp}, is outside the indicated temperature limits')

    return cfts

