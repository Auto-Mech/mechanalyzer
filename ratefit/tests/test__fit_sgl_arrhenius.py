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
TEMPS = RATEK_DATA.temps
RATEKS = RATEK_DATA.rateks

# Set NA and the T0 value in the (T/T0)^n term in the Arrhenius expr.
NA = 6.0221409e23
NA_INV = (1.0 / NA)
T_REF = 1.0
TMIN = None
TMAX = None

# Set paths to run the dsarrfit code
DSARRFIT_PATH = tempfile.mkdtemp()


def test__no_rates_to_fit():
    """ test ratefit.fit.util.get_valid_tk for no rates
    """

    # Filter the temperatures and rate constants to get valid values
    # k > 0 and k != *** and tmin <= T <= tmax
    temps, calc_ks = ratefit.fit.util.get_valid_tk(
        TEMPS[:3], RATEKS[:3], tmin=TMIN, tmax=TMAX, bimol=False)

    # Print header for python fitting
    assert calc_ks.size == 0 and temps.size == 0


def test__single_arrhenius_fit_python():
    """ test ratefit.fit.arrhenius.single with Python fitter
    """

    # Filter the temperatures and rate constants to get valid values
    temps, calc_ks = ratefit.fit.util.get_valid_tk(
        TEMPS, RATEKS, tmin=TMIN, tmax=TMAX, bimol=False)

    ref_temps = [
        600., 700., 800., 900., 1000., 1100., 1200.,
        1300., 1400., 1500., 1600., 1700., 1800., 1900.,
        2000., 2100., 2200., 2300., 2400., 2500., 2600.,
        2700., 2800., 2900., 3000.
    ]
    ref_calc_ks = [
        9.35678e+01, 2.28702e+03, 2.71034e+04,
        1.93348e+05, 9.55781e+05, 3.60000e+06,
        1.10000e+07, 2.85000e+07, 6.51000e+07,
        1.34000e+08, 2.52000e+08, 4.43000e+08,
        7.34000e+08, 1.16000e+09, 1.74000e+09,
        2.53000e+09, 3.56000e+09, 4.86000e+09,
        6.48000e+09, 8.45000e+09, 1.08000e+10,
        1.36000e+10, 1.68000e+10, 2.06000e+10,
        2.48000e+10
    ]

    assert numpy.allclose(temps, ref_temps)
    assert numpy.allclose(calc_ks, ref_calc_ks)

    # Run a single Arrhenius fit
    fit_params = ratefit.fit.arrhenius.single(
        temps, calc_ks, T_REF, 'python')

    ref_fit_params = [30559230.626540944,
                      1.3816123272407292,
                      25.773337886526104]

    assert numpy.allclose(fit_params, ref_fit_params)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.calc.single_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        T_REF, temps)

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
    mean_avg_err, max_avg_err = ratefit.calc.fitting_errors(
        calc_ks, fit_ks)

    ref_mean_avg_err = 2.504321595549729
    ref_max_avg_err = 7.827700474963309

    assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    assert numpy.allclose(max_avg_err, ref_max_avg_err)


def test__single_arrhenius_fit_dsarrfit():
    """ test ratefit.fit.arrhenius.single using dsarrfit code
    """

    # Filter the temperatures and rate constants to get valid values
    temps, calc_ks = ratefit.fit.util.get_valid_tk(
        TEMPS, RATEKS, tmin=TMIN, tmax=TMAX, bimol=False)

    fit_params = ratefit.fit.arrhenius.single(
        temps, calc_ks, T_REF, 'dsarrfit',
        dsarrfit_path=DSARRFIT_PATH,
        a_conv_factor=1.0)

    ref_fit_params = [30543100.0,
                      1.38168,
                      25.773244352868108]

    assert numpy.allclose(fit_params, ref_fit_params)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.calc.single_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        T_REF, temps)

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
    mean_avg_err, max_avg_err = ratefit.calc.fitting_errors(
        calc_ks, fit_ks)

    ref_mean_avg_err = 2.504192817017522
    ref_max_avg_err = 7.82923446743916

    assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    assert numpy.allclose(max_avg_err, ref_max_avg_err)


if __name__ == '__main__':
    test__no_rates_to_fit()
    test__single_arrhenius_fit_python()
    test__single_arrhenius_fit_dsarrfit()
