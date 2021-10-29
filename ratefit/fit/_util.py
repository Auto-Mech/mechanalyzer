"""
utility functions used for the fitting functions
"""


import numpy as np
from phydat import phycon


def filter_ktp_dct(inp_ktp_dct, bimol, tmin=None, tmax=None):
    """ Filters out negative and undefined rates from a ktp dictionary
    """

    filt_ktp_dct = {}
    for pressure, kt_lst in inp_ktp_dct.items():
        temps, kts = kt_lst[0], kt_lst[1]
        filt_temps, filt_kts = get_valid_tk(
            temps, kts, bimol, tmin=tmin, tmax=tmax)
        if filt_kts.size > 0:
            filt_ktp_dct[pressure] = (filt_temps, filt_kts)

    return filt_ktp_dct


def get_valid_tk(temps, rate_kts, bimol,
                 tmin=None, tmax=None, bimolthresh=1.0e-24):
    """ Takes in lists of temperature-rate constant pairs [T,k(T)]
        and removes invalid pairs for which
        (1) k(T) < 0
        (2) k(T) is undefined from Master Equation (i.e. k(T) is None)
        (3) k(T) < 1.0e-24 for a bimolecular reaction, or
        (4) T is outside the cutoff

        :param temps: Temperatures (K)
        :type temps: list(float)
        :param rate_kts: rate constants (s-1 or cm^3.s-1)
        :type rate_kts: list(str, float)
        :param numpy.ndarray temps: temps
        :param numpy.ndarray rate_constants: rate constants
        :param bool bimol: Parameter indicating bimolecular reation
        :param float tmin: minimum temperature cutoff for valid T,k(T) pairs
        :param float tmax: maximum temperature cutoff for valid T,k(T) pairs
        :return valid_t: List of vaild temperatures
        :return valid_k: List of vaild rate constants
        :rtype: (numpy.ndarray, numpy.ndarray)
        """

    # Check inputs
    # assert tmin in temps and tmax in temps

    # Set max temperature to user input, if none use max of input temperatures
    if tmax is None:
        tmax = max(temps)
    else:
        assert tmax in temps, ('{} not in temps: {}'.format(tmax, temps))

    # Set min temperature to user input, if none use either
    # min of input temperatures or
    # if negative kts are found, set min temp to be just above highest neg.
    if tmin is None:
        max_neg_idx = None
        for kt_idx, rate_kt in enumerate(rate_kts):
            # find idx for max temperature for which kt is negative, if any
            frate_kt = float(rate_kt) if rate_kt is not None else 1.0
            if frate_kt < 0.0:
                max_neg_idx = kt_idx
        # Set tmin to highest T where k(T) is non-negative
        # Set tmin to None (return no valid T,k(T)) if highest T has neg. k(T)
        if max_neg_idx is not None:
            if max_neg_idx+1 < len(temps):
                tmin = temps[max_neg_idx+1]
        else:
            tmin = min(temps)
    else:
        assert tmin in temps, ('{} not in temps: {}'.format(tmin, temps))

    # Grab the temperature, rate constant pairs which correspond to
    # temp > 0, temp within tmin and tmax, rate constant defined (not ***)
    valid_t, valid_k = [], []
    if tmin is not None and tmax is not None:
        for temp, rate_kt in zip(temps, rate_kts):
            if rate_kt is None:
                continue
            kthresh = 0.0 if not bimol else bimolthresh
            if float(rate_kt) > kthresh and tmin <= temp <= tmax:
                valid_t.append(temp)
                valid_k.append(rate_kt)

    # Convert the lists to numpy arrays
    valid_t = np.array(valid_t, dtype=np.float64)
    valid_k = np.array(valid_k, dtype=np.float64)

    return valid_t, valid_k


def invert_ktp_dct(ktp_dct):
    """ Invert the keys and values of a k(T,P) dictionary
        such that the dct changes it index from pressures to temperatures:

        ktp_dct[press] = [[t1, k1], ... , [tn, kn]]
        inv_ktp_dct[temp] = [[p1, k1], ... , [pn, kn]]

        :param ktp_dct: rates dictionary indexed by pressure
        :type ktp_dct: dict[float: tuple(numpy.ndarray, numpy.ndarray)]
        :return inv_ktp_dct: rates dictionary indexed by temperature
        :rtype ktp_dct: dict[float: tuple(numpy.ndarray, numpy.ndarray)]
    """

    inv_ktp_dct = {}
    for pressure, tk_arr in ktp_dct.items():

        # Set the temperatures and rate constants
        temps = tk_arr[0]
        rate_constants = tk_arr[1]

        for temp, rate in zip(temps, rate_constants):
            if temp not in inv_ktp_dct:
                # Set new temperature lst in dct
                inv_ktp_dct[temp] = ((pressure,), (rate,))
            else:
                # Add pressure & rate k to lsts for temp in dct
                p_arr, k_arr = inv_ktp_dct[temp]
                p_arr += (pressure,)
                k_arr += (rate,)
                inv_ktp_dct[temp] = (p_arr, k_arr)

    return inv_ktp_dct


def fitting_temperatures_dct(ktp_dct):
    """ Build a fit temp dct
    """

    _dct = {}
    for pressure, tk_arr in ktp_dct.items():
        temps = tk_arr[0]
        _dct[pressure] = (min(temps), max(temps))

    return _dct


def pull_highp_from_dct(param_dct):
    """ seperate the high pressure rates from the param
    """

    # Build Pressure dependent ktp dct
    pdep_dct = {}
    pressures = [pressure for pressure in param_dct.keys()
                 if pressure != 'high']
    for pressure in pressures:
        pdep_dct[pressure] = param_dct[pressure]

    # Get the high pressure parameters
    if 'high' in param_dct:
        highp_params = param_dct['high']
    else:
        highp_params = tuple()

    return highp_params, pdep_dct, pressures


def set_a_conversion_factor(reaction):
    """ Use reaction to set A conversion factor
    """
    [lab_i, _] = reaction.split('=')
    return phycon.NAVO if 'W' not in lab_i else 1.00


def get_temps_and_pressures(ktp_dct):
    """ Reads a ktp_dct and gets the list of pressure and corresponding list of
        temperature arrays

        :param ktp_dct: k(T,P) values
        :type ktp_dct: dict {pressure: (temps, kts)}
        :return temps_lst: list of temperature arrays at each pressure (K)
        :rtype: list [numpy.ndarray1, numpy.ndarray2, ...]
        :return pressures: list of pressures (atm)
        :rtype: list
    """

    temps_lst = []
    pressures = []
    for pressure, (temps, kts) in ktp_dct.items():
        temps_lst.append(temps)
        pressures.append(pressure)
        
    return temps_lst, pressures
