"""
Calculate rates with various fitting functions
"""

import numpy as np


def calc_sse_and_mae(calc_ks, fit_ks):
    """ (1) get the sum of square error (SSE) useful when determining
            which double plog routine will be used to initialize
            the nonlinear solver
        (2) also get the mean absolute error (MAE), which is written
            to the plog file
        Only need to assess error if there are 2 or more rate constants 
    """

    abs_err = []
    if len(calc_ks) > 2:
        for calc_k, fit_k in zip(calc_ks, fit_ks):
            abs_err.append(np.abs((calc_k - fit_k) / calc_k))
        abs_err = np.array(abs_err, dtype=np.float64)
        mean_abs_err = np.mean(abs_err) * 100.0
        max_abs_err = np.max(abs_err) * 100.0
    else:
        mean_abs_err = 0.0
        max_abs_err = 0.0

    return mean_abs_err, max_abs_err


def assess_pressure_dependence(tk_dct, assess_pdep_temps,
                               tolerance=20.0, plow=None, phigh=None):
    """ Assess how much the rate constants change from
        a low-pressure to high-pressure regime

        tk_dct[pressure] = [temps, k(T, P)s]
        we assume the temps and pressures give all positive, defined rates
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
            
        for temp_compare in assess_pdep_temps:
            # For the low- and high-P, find the idx for the temp in temp_compare
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


if __name__ == '__main__':
    PAIRS = [
        [800, 27103.4],
        [900, 193348],
        [1000, 955781],
        [1100, 3.60E+06],
        [1200, 1.10E+07],
        [1300, 2.85E+07],
        [1400, 6.51E+07],
        [1500, 1.34E+08],
        [1600, 2.52E+08],
        [1700, 4.43E+08],
        [1800, 7.34E+08],
        [1900, 1.16E+09],
        [2000, 1.74E+09],
        [2100, 2.53E+09],
        [2200, 3.56E+09],
        [2300, 4.86E+09],
        [2400, 6.48E+09],
        [2500, 8.45E+09],
        [2600, 1.08E+10],
        [2700, 1.36E+10],
        [2800, 1.68E+10],
        [2900, 2.06E+10],
        [3000, 2.48E+10]
    ]
    
    KTP_DCT = {
    #     1:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
    #     10:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
    #      10:  [[pair[0] for pair in PAIRS2], [pair[1] for pair in PAIRS2]],
          10:  np.array([[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]]),
    #     100:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
    #     'high':  numpy.array([[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]]),
    }
    
    assess_pressure_dependence(KTP_DCT, [1000.0],
                               tolerance=20.0, plow=None, phigh=None)
