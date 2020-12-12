"""
Test an Arrhenius fit of T, k(T,P) rates to a single Arrhenius function
where the fits are performed using Python and dsarrfit
"""

import os
import tempfile
import numpy
import pandas
import ratefit


def _read_csv(filename):
    """ read csv values from arrhenius.csv
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()

    return data


# Set path to data files
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
RATEK_FILE_NAME = 'rateks.csv'

# Read csv file for data
RATEK_DATA = _read_csv(
    os.path.join(DATA_PATH, RATEK_FILE_NAME))

# Set the rate constant data
TEMPS = RATEK_DATA.temps.values
RATEKS = RATEK_DATA.rateks.values

# Set NA and the T0 value in the (T/T0)^n term in the Arrhenius expr.
NA = 6.0221409e23
NA_INV = (1.0 / NA)
T_REF = 1.0
TMIN = None
TMAX = None

# Set paths to run the dsarrfit code
DSARRFIT_PATH = tempfile.mkdtemp()


def test__single_arrhenius_fit_python():
    """ test ratefit.fit.arrhenius.single with Python fitter
    """

    # Run a single Arrhenius fit
    fit_params = ratefit.fit.arrhenius.single(
        TEMPS, RATEKS, T_REF, 'python')

    ref_fit_params = [30559230.626540944,
                      1.3816123272407292,
                      25.773337886526104]

    assert numpy.allclose(fit_params, ref_fit_params)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.calc.single_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        T_REF, TEMPS)

    ref_fit_ks = [
        8.62436E+01, 2.34060E+03, 2.85285E+04, 2.03360E+05,
        9.93878E+05, 3.68621E+06, 1.11047E+07, 2.84839E+07,
        6.43507E+07, 1.31272E+08, 2.46371E+08, 4.31568E+08,
        7.13541E+08, 1.12346E+09, 1.69654E+09, 2.47143E+09,
        3.48964E+09, 4.79480E+09, 6.43203E+09, 8.44734E+09,
        1.08870E+10, 1.37972E+10, 1.72234E+10, 2.12101E+10,
        2.58005E+10
    ]

    assert numpy.allclose(fit_ks, ref_fit_ks)

    # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
        RATEKS, fit_ks)

    ref_mean_avg_err = 2.504321595549729
    ref_max_avg_err = 7.827700474963309

    assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    assert numpy.allclose(max_avg_err, ref_max_avg_err)


def test__single_arrhenius_fit_dsarrfit():
    """ test ratefit.fit.arrhenius.single using dsarrfit code
    """

    fit_params = ratefit.fit.arrhenius.single(
        TEMPS, RATEKS, T_REF, 'dsarrfit',
        dsarrfit_path=DSARRFIT_PATH,
        a_conv_factor=1.0)

    ref_fit_params = [30543100.0,
                      1.38168,
                      25.773244352868108]

    assert numpy.allclose(fit_params, ref_fit_params)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.calc.single_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        T_REF, TEMPS)

    ref_fit_ks = [
        8.62422E+01, 2.34056E+03, 2.85281E+04, 2.03356E+05,
        9.93865E+05, 3.68617E+06, 1.11046E+07, 2.84838E+07,
        6.43504E+07, 1.31272E+08, 2.46372E+08, 4.31569E+08,
        7.13545E+08, 1.12347E+09, 1.69655E+09, 2.47146E+09,
        3.48969E+09, 4.79488E+09, 6.43215E+09, 8.44751E+09,
        1.08873E+10, 1.37976E+10, 1.72239E+10, 2.12107E+10,
        2.58012E+10
    ]

    assert numpy.allclose(fit_ks, ref_fit_ks)

    # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
        RATEKS, fit_ks)

    ref_mean_avg_err = 2.504192817017522
    ref_max_avg_err = 7.82923446743916

    assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    assert numpy.allclose(max_avg_err, ref_max_avg_err)


if __name__ == '__main__':
    test__single_arrhenius_fit_python()
    test__single_arrhenius_fit_dsarrfit()
