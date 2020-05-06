"""
Test the ratefit rate constant calculators
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
ARR_FILE_NAME = 'arrhenius.csv'

# Read csv file for data
ARR_K_DATA = _read_csv(
    os.path.join(DATA_PATH, ARR_FILE_NAME))

# Set the data for the calculations
TEMPS = numpy.array(
    [300., 600., 900., 1200., 1500.,
     1800., 2100., 2400., 2700., 3000.])
T_REF = 1.0
A1, N1, EA1 = 7.015e+4, 2.053, -0.3557
A2, N2, EA2 = 5.757e+12, -0.664, 0.3318

numpy.set_printoptions(precision=15)


def test__single_arrhenius():
    """ test ratefit.fxns.single_arrhenius
    """

    calc_ks = ratefit.calc.single_arrhenius(
        A1, N1, EA1,
        T_REF, TEMPS)

    ref_calc_ks = numpy.array(ARR_K_DATA.SingleArr)

    assert numpy.allclose(calc_ks, ref_calc_ks, atol=0.01)


def test__double_arrhenius():
    """ test ratefit.fxns.double_arrhenius
    """

    calc_ks = ratefit.calc.double_arrhenius(
        A1, N1, EA1,
        A2, N2, EA2,
        T_REF, TEMPS)

    ref_calc_ks = numpy.array(ARR_K_DATA.DoubleArr)

    assert numpy.allclose(calc_ks, ref_calc_ks, atol=0.01)


if __name__ == '__main__':
    test__single_arrhenius()
    test__double_arrhenius()
