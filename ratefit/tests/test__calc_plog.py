"""
Test the ratefit rate constant calculators
"""

import os
import numpy
import pandas
import ratefit


def _read_csv(filename):
    """ read csv values from file
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


# Set path to data files
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
PLOG_FILE_NAME = 'plog.csv'

# Read csv file for data
PLOG_K_DATA = _read_csv(
    os.path.join(DATA_PATH, PLOG_FILE_NAME))

# Set data for the calculations
TEMPS = numpy.array(
    [300., 600., 900., 1200., 1500.,
     1800., 2100., 2400., 2700., 3000.])
PRESSURES = numpy.array([0.1, 0.9869, 2.0, 5.0])
T_REF = 1.0

PLOG_DCT = {
    0.0296: [[2.020E+013, -1.870, 22.755]],
    0.0987: [[1.680E+018, -3.050, 24.323]],
    0.2961: [[2.500E+024, -4.630, 27.067]],
    0.9869: [[4.540E+026, -5.120, 27.572]],
    2.9607: [[7.120E+028, -5.600, 28.535]],
    9.8690: [[5.480E+029, -5.700, 28.899]]
}

numpy.set_printoptions(precision=15)


def test__plog():
    """ test ratefit.calc.plog
    """

    plog_ktps = ratefit.calc.plog(PLOG_DCT, T_REF, TEMPS, PRESSURES)

    plog_ktps1 = plog_ktps[0.1]
    plog_ktps2 = plog_ktps[0.9869]
    plog_ktps3 = plog_ktps[2.0]
    plog_ktps4 = plog_ktps[5.0]

    assert numpy.allclose(plog_ktps1, PLOG_K_DATA.ktp1.values, atol=0.01)
    assert numpy.allclose(plog_ktps2, PLOG_K_DATA.ktp2.values, atol=0.01)
    assert numpy.allclose(plog_ktps3, PLOG_K_DATA.ktp3.values, atol=0.01)
    assert numpy.allclose(plog_ktps4, PLOG_K_DATA.ktp4.values, atol=0.01)


if __name__ == '__main__':
    test__plog()
