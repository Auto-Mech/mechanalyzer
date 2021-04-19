"""
Test the ratefit rate constant calculators
"""

import os
import numpy
import pandas
import ratefit


def _read_csv(filename):
    """ read csv values
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


# Set path to data files
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
CHEB_FILE_NAME = 'chebyshev.csv'

# Read csv file for data
# CHEB_K_DATA = _read_csv(
#     os.path.join(DATA_PATH, CHEB_FILE_NAME))

# Set the data for the calculations
# TEMPS = numpy.array(
#     [300., 600., 900., 1200., 1500.,
#      1800., 2100., 2400., 2700., 3000.])
# PRESSURES = numpy.array([0.1, 0.9869, 2.0, 5.0])
# T_REF = 1.0
#
# numpy.set_printoptions(precision=15)

CHEB_K_DATA = ()
HIGH_PARAMS = [1.00, 0.00, 0.00]
TMIN, TMAX = (300.000, 2200.000)
PMIN, PMAX = (0.01, 20.0)
# PMIN, PMAX = (0.010, 98.702)
TEMPS = numpy.arange(TMIN, TMAX+100.0, 100.0)
PRESSURES = (0.1, 1.0, 10.0, 20.0, 40.0, 80, 98.702)
ALPHA = numpy.array([
    [8.684e+00, 7.500e-01, -7.486e-02, 1.879e-15],
    [-2.159e-01, 9.899e-02, 2.292e-02, 2.929e-17],
    [-1.557e-15, -3.331e-16, 3.324e-17, -8.346e-31],
    [2.159e-01, -9.899e-02, -2.292e-02, -2.929e-17],
    [-2.684e+00, -7.500e-01, 7.486e-02, -1.879e-15],
    [2.159e-01, -9.899e-02, -2.292e-02, -2.929e-17]
])


def test__chebyshev():
    """ test ratefit.fxns.chebyshev
    """

    cheb_ktps = ratefit.calc.chebyshev(
        ALPHA, TMIN, TMAX, PMIN, PMAX, TEMPS, PRESSURES)

    cheb_ktps1 = cheb_ktps[0.1]
    cheb_ktps2 = cheb_ktps[0.9869]
    cheb_ktps3 = cheb_ktps[2.0]
    cheb_ktps4 = cheb_ktps[5.0]

    assert numpy.allclose(cheb_ktps1,numpy.array(CHEB_K_DATA.ktp1), atol=0.01)
    assert numpy.allclose(cheb_ktps2,numpy.array(CHEB_K_DATA.ktp2), atol=0.01)
    assert numpy.allclose(cheb_ktps3,numpy.array(CHEB_K_DATA.ktp3), atol=0.01)
    assert numpy.allclose(cheb_ktps4,numpy.array(CHEB_K_DATA.ktp4), atol=0.01)


if __name__ == '__main__':
    test__chebyshev()
