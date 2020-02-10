"""
utility functions used for the fitting functions
"""

import numpy as np


def get_valid_tk(temps, rate_constants, bimol,
                 tmin=None, tmax=None):
    """ this subroutine takes in a array of rate constants and
        returns the subset of this array that is positive,
        along with the corresponding Temperature array """

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
        else:
            #print('bimol:', bimol, 'rate_constant', float(rate_constant))
            if not bimol and float(rate_constant) > 0.0 and tmin <= temp <= tmax:
                valid_t.append(temp)
                valid_k.append(rate_constant)
            if bimol and float(rate_constant) > 1.e-21 and tmin <= temp <= tmax:
                valid_t.append(temp)
                valid_k.append(rate_constant)

    # Convert the lists to numpy arrays
    valid_t = np.array([valid_t], dtype=np.float64)
    valid_k = np.array([valid_k], dtype=np.float64)

    return valid_t, valid_k
