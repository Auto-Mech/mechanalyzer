"""
Calculate the errors introduced by rates using fitting them to
some functional formwith various fitting functions
"""

import numpy as np


def fitting_errors(calc_ks, fit_ks):
    """ Calculates the error associated with fitting a set of
        rate constants [k(T,P)s] to a functional form.

        :param calc_ks: Original k(T,P)s fit to functional form
        :type calc_ks: numpy.ndarray
        :param fit_ks: k(T,P)s calc'd from a functional form
        :type fit_ks: numpy.ndarray
        :return mean_abs_err: Mean absoulte error from fit
        :rtype: float
        :return max_abs_err: Maximum absoulte error from fit
        :rtype: float
    """

    assert len(calc_ks) == len(fit_ks)

    abs_err = []
    # if len(calc_ks) > 2:  # error SHOULD be zero for one ktp, but will run
    for calc_k, fit_k in zip(calc_ks, fit_ks):
        abs_err.append(np.abs((calc_k - fit_k) / calc_k))
    abs_err = np.array(abs_err, dtype=np.float64)
    mean_abs_err = np.mean(abs_err) * 100.0
    max_abs_err = np.max(abs_err) * 100.0
    # else:
    #     mean_abs_err = 0.0
    #     max_abs_err = 0.0

    return mean_abs_err, max_abs_err
