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
CHEB_K_DATA = _read_csv(
    os.path.join(DATA_PATH, CHEB_FILE_NAME))

# Set the data for the calculations
# TEMPS = numpy.array(
#     [300., 600., 900., 1200., 1500.,
#      1800., 2100., 2400., 2700., 3000.])
# PRESSURES = numpy.array([0.1, 0.9869, 2.0, 5.0])
# T_REF = 1.0
#
# numpy.set_printoptions(precision=15)


def test__chebyshev():
    """ test ratefit.fxns.chebyshev
    """
    pass
    # cheb_ktps = ratefit.calc.chebyshev(
    #     ALPHA, TMIN, TMAX, PMIN, PMAX, PRESSURES, TEMPS)

    # cheb_ktps1 = cheb_ktps[0.1]
    # cheb_ktps2 = cheb_ktps[0.9869]
    # cheb_ktps3 = cheb_ktps[2.0]
    # cheb_ktps4 = cheb_ktps[5.0]

    # assertnumpy.allclose(cheb_ktps1,numpy.array(CHEB_K_DATA.ktp1), atol=0.01)
    # assertnumpy.allclose(cheb_ktps2,numpy.array(CHEB_K_DATA.ktp2), atol=0.01)
    # assertnumpy.allclose(cheb_ktps3,numpy.array(CHEB_K_DATA.ktp3), atol=0.01)
    # assertnumpy.allclose(cheb_ktps4,numpy.array(CHEB_K_DATA.ktp4), atol=0.01)


if __name__ == '__main__':
    test__chebyshev()
