"""
Calculate the errors introduced by rates using fitting them to
some functional formwith various fitting functions
"""

import numpy as np


def fitting_errors(calc_ks, fit_ks):
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
