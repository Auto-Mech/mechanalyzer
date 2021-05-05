"""
Calculate the errors introduced by rates using fitting them to
some functional formwith various fitting functions
"""

from statistics import mean


def fitting_error_dct(calc_ktp_dct, fit_ktp_dct, err_set='all'):
    """ Calculate the fitting errors for two k(T,P) dictionaries
        consisting of calculated and fitted values at each pressure.

        :param calc_ktp_dct: calculate rate k(T,P)s
        :type calc_ktp_dct: dict[float: (tuple(float), tuple(float))
        :param fit_ktp_dct: fitted rate k(T,P)s
        :type fit_ktp_dct: dict[float: (tuple(float), tuple(float))
        :rtype: dict[float: (float, float)]
    """

    assert set(calc_ktp_dct.keys()) == set(fit_ktp_dct.keys()), (
        'pressures of two k(T,P) dcts are not the same')
    assert all(len(calc_ktp_dct[p][0]) == len(fit_ktp_dct[p][0])
               for p in calc_ktp_dct.keys()), ('temps not the same')

    _fit_err_dct = {}
    for pressure in calc_ktp_dct:

        # Get k(T) values at pressure
        calc_ks = calc_ktp_dct[pressure][1]
        fit_ks = fit_ktp_dct[pressure][1]

        # Assess the errors using some subset of the rate constants
        test_calc_ks, test_fit_ks = _gen_err_set(
            calc_ks, fit_ks, err_set=err_set)
        mean_avg_err, max_avg_err = fitting_errors(
            test_calc_ks, test_fit_ks)

        # Store in a dictionary
        _fit_err_dct[pressure] = (mean_avg_err, max_avg_err)

    return _fit_err_dct


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
    for calc_k, fit_k in zip(calc_ks, fit_ks):
        abs_err.append(abs((calc_k - fit_k) / calc_k))

    mean_abs_err = mean(abs_err) * 100.0
    max_abs_err = max(abs_err) * 100.0

    return mean_abs_err, max_abs_err


def _gen_err_set(calc_ks, fit_ks, err_set='all'):
    """ Set the ranges of the k(T) values that are
        used to assess the fitting errors.

        :param calc_ks: calculated k(T) values
        :type calc_ks: numpy.ndarray
        :param calc_ks: fitted k(T) values
        :type calc_ks: numpy.ndarray
        :rtype: (numpy.ndarray, numpy.ndarray)
    """

    assert len(calc_ks) == len(fit_ks)

    if err_set == 'all':
        test_calc_ks = calc_ks
        test_fit_ks = fit_ks
    else:
        test_calc_ks = calc_ks[1:-2]
        test_fit_ks = fit_ks[1:-2]

    return test_calc_ks, test_fit_ks
