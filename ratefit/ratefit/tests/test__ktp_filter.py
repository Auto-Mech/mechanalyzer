"""
Test an Arrhenius fit of T, k(T,P) rates to a single Arrhenius function
where the fits are performed using Python and dsarrfit
"""

import os
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
RATEK_FILE_NAME = 'full_rateks.csv'

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


def test__valid_tk_no_rates():
    """ test ratefit.fit.util.get_valid_tk for no rates
    """

    # Filter the temperatures and rate constants to get valid values
    # k > 0 and k != *** and tmin <= T <= tmax
    temps, calc_ks = ratefit.fit.get_valid_tk(
        TEMPS[:3], RATEKS[:3], tmin=TMIN, tmax=TMAX, bimol=False)

    # Print header for python fitting
    assert calc_ks.size == 0 and temps.size == 0


def test__valid_tk_with_rates():
    """ test ratefit.fit.util.get_valid_tk
    """

    temps, calc_ks = ratefit.fit.get_valid_tk(
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


if __name__ == '__main__':
    test__valid_tk_no_rates()
    test__valid_tk_with_rates()
