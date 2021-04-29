"""
Test the ratefit rate constant calculators for Lindemann and Troe expressions
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
TROE_FILE_NAME = 'troe.csv'

# Read csv file for data
TROE_K_DATA = _read_csv(
    os.path.join(DATA_PATH, TROE_FILE_NAME))

# Set the data for the calculations
TEMPS = numpy.array(
    [300., 600., 900., 1200., 1500.,
     1800., 2100., 2400., 2700., 3000.])
PRESSURES = numpy.array([0.1, 2.0, 5.0, 10.0])
T_REF = 1.0

A_HIGH, N_HIGH, EA_HIGH = 2.000e+12, 0.900, 4.87490
A_LOW, N_LOW, EA_LOW = 2.490e24, -2.300, 4.87490
TROE_ALPHA, TROE_T3, TROE_T1, TROE_T2 = 6.0e-1, 1.0e3, 7.0, 1.7e3

numpy.set_printoptions(precision=15)


def test__calc():
    """ test ratefit.calc.troe
    """

    # Calc high-pressure and low-pressure rates using Arrhenius fits
    # Should replace with rate constants
    highp_ks = ratefit.calc.single_arrhenius(
        A_HIGH, N_HIGH, EA_HIGH,
        T_REF, TEMPS)

    lowp_ks = ratefit.calc.single_arrhenius(
        A_LOW, N_LOW, EA_LOW,
        T_REF, TEMPS)

    # Calc the rate contants using the Lindemann functional form
    troe_ktps = ratefit.calc.troe(
        highp_ks, lowp_ks,
        TEMPS, PRESSURES,
        TROE_ALPHA, TROE_T3, TROE_T1, ts2=TROE_T2, collid_factor=1.0)

    troe_ktps1 = troe_ktps[0.1]
    troe_ktps2 = troe_ktps[2.0]
    troe_ktps3 = troe_ktps[5.0]
    troe_ktps4 = troe_ktps[10.0]

    assert numpy.allclose(troe_ktps1, TROE_K_DATA.ktp1.values, atol=0.01)
    assert numpy.allclose(troe_ktps2, TROE_K_DATA.ktp2.values, atol=0.01)
    assert numpy.allclose(troe_ktps3, TROE_K_DATA.ktp3.values, atol=0.01)
    assert numpy.allclose(troe_ktps4, TROE_K_DATA.ktp4.values, atol=0.01)


if __name__ == '__main__':
    test__calc()
