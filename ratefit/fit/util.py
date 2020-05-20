"""
utility functions used for the fitting functions
"""

import numpy as np


def get_valid_tk(temps, rate_constants, bimol,
                 tmin=None, tmax=None):
    """ Takes in lists of temperature-rate constant pairs [T,k(T)] 
        and removes invalid pairs for which
           (1) k(T) < 0
           (2) k(T) is undefined from Master Equation (i.e. k(T) == '***')
           (3) k(T) < 1.0e-21 for a bimolecular reaction, or
           (4) T is outside the cutoff
        :param numpy.ndarray temps: temps
        :param numpy.ndarray rate_constants: rate constants
        :param bool bimol: Parameter indicating bimolecular reation
        :param float tmin: minimum temperature cutoff for valid T,k(T) pairs
        :param float tmax: maximum temperature cutoff for valid T,k(T) pairs
        :return valid_t: List of vaild temperatures
        :rtype numpy.ndarray
        :return valid_k: List of vaild rate constants
        :rtype numpy.ndarray
        """

    # Convert temps and rate constants to floats
    temps = [float(temp) for temp in temps]
    rate_constants = [float(rate_constant)
                      if rate_constant != '***' else rate_constant
                      for rate_constant in rate_constants]

    # Set tmin and tmax
    if tmin is None:
        tmin = min(temps)
    if tmax is None:
        tmax = max(temps)
    assert tmin in temps and tmax in temps

    # Grab the temperature, rate constant pairs which correspond to
    # temp > 0, temp within tmin and tmax, rate constant defined (not ***)
    valid_t, valid_k = [], []
    for temp, rate_constant in zip(temps, rate_constants):
        if rate_constant == '***':
            continue
        kthresh = 0.0 if not bimol else 1.0e-21
        if float(rate_constant) > kthresh and tmin <= temp <= tmax:
            valid_t.append(temp)
            valid_k.append(rate_constant)

    # Convert the lists to numpy arrays
    valid_t = np.array(valid_t, dtype=np.float64)
    valid_k = np.array(valid_k, dtype=np.float64)

    return valid_t, valid_k


def flip_ktp_dct(ktp_dct):
    """ Invert the keys and values of a k(T,P) dictionary
        such that the dct changes it index from pressures to temperatures:
            ktp[temp] = [[p1, k1], ... , [pn, kn]]
        :param dct ktp_dct: ktp_dct[press] = [[t1, k1], ... , [tn, kn]]
        :return inv_ktp_dct: ktp_dct[temp] = [[p1, k1], ... , [pn, kn]]
        :rtype: dct
    """

    inv_ktp_dct = {}
    for pressure, tk_arr in ktp_dct.items():

        # Set the temperatures and rate constants
        temps = tk_arr[0]
        rate_constants = tk_arr[1]

        for temp, rate in zip(temps, rate_constants):
            if temp not in inv_ktp_dct:
                # Set new temperature lst in dct
                inv_ktp_dct[temp] = [[pressure], [rate]]
            else:
                # Add pressure & rate k to lsts for temp in dct
                [p_arr, k_arr] = inv_ktp_dct[temp]
                p_arr.append(pressure)
                k_arr.append(rate)
                inv_ktp_dct[temp] = [p_arr, k_arr]

    return inv_ktp_dct
