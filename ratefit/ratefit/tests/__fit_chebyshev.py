"""
Test an Arrhenius fit of T, k(T,P) rates to a single Arrhenius function
where the fits are performed using Python and dsarrfit
"""

import os
import tempfile
import numpy
import pandas
import ratefit


# Set NA and the T0 value in the (T/T0)^n term in the Arrhenius expr.
NA = 6.0221409e23
NA_INV = (1.0 / NA)
TDEG = 6
PDEG = 4

TMIN, TMAX = (300.000, 2200.000)
PMIN, PMAX = (0.010, 98.702)
TEMPS = numpy.arange(TMIN, TMAX, 100.0)
PRESSURES = (0.1, 1.0, 10.0, 20.0)



def test__chebyshev_fit():
    """ test ratefit.fit.troe.std
    """

    # Run a Chebyshev fit
    alpha, trange, prange = ratefit.fit.chebyshev.kfit(
        TEMPS, RATEKS, tdeg=TDEG, pdeg=PDEG)
    print('alpha\n', alpha)
    # assert numpy.allclose(fit_params, ref_fit_params)

    # # Calculate fitted rate constants using the fitted parameters
    fit_ktps = ratefit.calc.chebyshev(
        alpha, TMIN, TMAX, PMIN, PMAX, TEMPS, PRESSURES)
    # print('fit ks', fit_ks)
    # assert numpy.allclose(fit_ks, ref_fit_ks)

    # # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_errs, max_avg_errs = [], []
    for pressure in PRESSURES:
        rate_kts = RATEKS[pressure][1]
        fit_kts = numpy.array(fit_ktps[pressure])
        mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
            rate_kts, fit_kts)
        mean_avg_errs.append(mean_avg_err)
        max_avg_errs.append(max_avg_err)
        print('rate kts')
        print(rate_kts)
        print('fit kts')
        print(fit_kts)
        print('errs')
        print(mean_avg_err, max_avg_err)

    # ref_mean_avg_err = 2.504321595549729
    # ref_max_avg_err = 7.827700474963309

    # assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    # assert numpy.allclose(max_avg_err, ref_max_avg_err)

    assert True == False


if __name__ == '__main__':
    test__chebyshev_fit()
