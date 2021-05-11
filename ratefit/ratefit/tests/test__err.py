""" test error
"""

import numpy
from ratefit.fit import fitting_error_dct


CALC_KTP_DCT = {
    0.1: ((250.0, 500.0, 750.0, 1000.0, 1250.0), (1.0, 2.0, 3.0, 4.0, 5.0)),
    1.0: ((250.0, 500.0, 750.0, 1000.0, 1250.0), (2.0, 3.0, 4.0, 5.0, 6.0)),
    10.: ((250.0, 500.0, 750.0, 1000.0, 1250.0), (3.0, 4.0, 5.0, 6.0, 7.0))
}
FIT_KTP_DCT = {
    0.1: ((250.0, 500.0, 750.0, 1000.0, 1250.0), (1.2, 2.2, 3.2, 4.2, 5.2)),
    1.0: ((250.0, 500.0, 750.0, 1000.0, 1250.0), (2.5, 3.5, 4.5, 5.5, 6.5)),
    10.: ((250.0, 500.0, 750.0, 1000.0, 1250.0), (3.8, 4.8, 5.8, 6.8, 7.8))
}


def test__ktp_err():
    """ test ratefit.fit.err
    """

    ref_fit_err_dct1 = {
        0.1: (9.133333333333336, 19.999999999999996),
        1.0: (14.499999999999998, 25.0),
        10.0: (17.48571428571428, 26.66666666666666)
    }
    ref_fit_err_dct2 = {
        0.1: (8.33333333333334, 10.000000000000009),
        1.0: (14.583333333333332, 16.666666666666664),
        10.0: (17.999999999999996, 19.999999999999996)
    }

    fit_err_dct1 = fitting_error_dct(
        CALC_KTP_DCT, FIT_KTP_DCT)
    fit_err_dct2 = fitting_error_dct(
        CALC_KTP_DCT, FIT_KTP_DCT, err_set='skip')

    for key1, key2 in zip(fit_err_dct1.keys(), ref_fit_err_dct1.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(fit_err_dct1[key1], ref_fit_err_dct1[key2])

    for key1, key2 in zip(fit_err_dct2.keys(), ref_fit_err_dct2.keys()):
        assert numpy.isclose(key1, key2)
        assert numpy.allclose(fit_err_dct2[key1], ref_fit_err_dct2[key2])
