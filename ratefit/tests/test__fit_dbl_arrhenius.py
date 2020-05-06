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
DATA_PATH = os.path.join(PATH, '../', 'data')
RATEK_FILE_NAME = 'rateks2.csv'

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
TMIN = 500.0
TMAX = None

# Set paths to run the dsarrfit code
DSARRFIT_PATH = tempfile.mkdtemp()
print('dsarffit', DSARRFIT_PATH)


def test__double_arrhenius_fit_python():
    """ test ratefit.fit.arrhenius.double
    """

    # Filter the temperatures and rate constants to get valid values
    temps, calc_ks = ratefit.fit.util.get_valid_tk(
        TEMPS, RATEKS, tmin=TMIN, tmax=TMAX, bimol=False)

    # Run a single Arrhenius fit to build a guess for the double fit
    sgl_fit = ratefit.fit.arrhenius.single(
        temps, calc_ks, T_REF, 'python')

    test_sgl_fit = [1.18698937208059e-15,
                    1.6600455223689599,
                    9.28586359726521]

    assert numpy.allclose(sgl_fit, test_sgl_fit)

    # Run a double Arrhenius fit
    fit_params = ratefit.fit.arrhenius.double(
        temps, calc_ks, T_REF, 'python',
        a_guess=sgl_fit[0],
        n_guess=sgl_fit[1],
        ea_guess=sgl_fit[2])

    test_fit_params = [
        8.036185028297988e-10, 0.0761979228320102, 16.66413018792614,
        1.8274992480059228e-13, 0.8886973965169201, 9.488903662805047
    ]

    assert numpy.allclose(fit_params, test_fit_params)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.calc.double_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        fit_params[3], fit_params[4], fit_params[5],
        T_REF, temps)

    test_fit_ks = [
        3.32484E-15, 2.87796E-14, 1.31244E-13, 4.13627E-13,
        1.02511E-12, 3.13576E-12, 7.26684E-12, 1.39700E-11,
        2.35473E-11, 3.60671E-11, 5.14228E-11, 7.92788E-11,
        1.36158E-10
    ]

    assert numpy.allclose(fit_ks, test_fit_ks)

    # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_err, max_avg_err = ratefit.calc.fitting_errors(
        calc_ks, fit_ks)

    test_mean_avg_err = 0.4741670536354967
    test_max_avg_err = 0.875017109912643

    assert numpy.allclose(mean_avg_err, test_mean_avg_err)
    assert numpy.allclose(max_avg_err, test_max_avg_err)


def test__double_arrhenius_fit_dsarrfit():
    """ test ratefit.fit.arrhenius.double
    """

    # Filter the temperatures and rate constants to get valid values
    temps, calc_ks = ratefit.fit.util.get_valid_tk(
        TEMPS, RATEKS, tmin=TMIN, tmax=TMAX, bimol=False)

    # Run a double Arrhenius fit
    fit_params = ratefit.fit.arrhenius.double(
        temps, calc_ks, T_REF, 'dsarrfit', dsarrfit_path=DSARRFIT_PATH,
        a_conv_factor=1.0)

    test_fit_params = [
        1.04925e-09, 0.0415606, 16.794141398369863,
        1.2337e-13, 0.944382, 9.440829863993441
    ]

    assert numpy.allclose(fit_params, test_fit_params)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.calc.double_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        fit_params[3], fit_params[4], fit_params[5],
        T_REF, temps)

    test_fit_ks = [
        3.32464E-15, 2.87785E-14, 1.31243E-13, 4.13605E-13,
        1.02501E-12, 3.13546E-12, 7.26645E-12, 1.39697E-11,
        2.35471E-11, 3.60664E-11, 5.14209E-11, 7.92739E-11,
        1.36153E-10
    ]

    assert numpy.allclose(fit_ks, test_fit_ks)

    # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_err, max_avg_err = ratefit.calc.fitting_errors(
        calc_ks, fit_ks)

    test_mean_avg_err = 0.47471982687376374
    test_max_avg_err = 0.8711794726859022

    assert numpy.allclose(mean_avg_err, test_mean_avg_err)
    assert numpy.allclose(max_avg_err, test_max_avg_err)


if __name__ == '__main__':
    test__double_arrhenius_fit_python()
    test__double_arrhenius_fit_dsarrfit()
