""" test ratefit.calc.lowp_limit
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
LOWPLIM_FILE_NAME = 'lowp_limit.csv'

# Read csv file for data
LOWPLIM_K_DATA = _read_csv(
    os.path.join(DATA_PATH, LOWPLIM_FILE_NAME))

TEMPS = numpy.array(
    [300., 600., 900., 1200., 1500.,
     1800., 2100., 2400., 2700., 3000.])
PRESSURES = numpy.array([0.1, 2.0, 5.0, 10.0])

HIGHP_RATEKS = ()


def test__calc():
    """ test ratefit.calc.rate.lowp_limit
        test ratefit.calc.rate.lowp_limit_one_pressure
    """

    lowp_lim_ktps = ratefit.calc.lowp_limit(HIGHP_RATEKS, TEMPS, PRESSURES)

    lowp_lim_ktps1 = lowp_lim_ktps[0.1]
    lowp_lim_ktps2 = lowp_lim_ktps[2.0]
    lowp_lim_ktps3 = lowp_lim_ktps[5.0]
    lowp_lim_ktps4 = lowp_lim_ktps[10.0]

    assert numpy.allclose(lowp_lim_ktps1, LOWPLIM_K_DATA.ktp1.values, atol=0.01)
    assert numpy.allclose(lowp_lim_ktps2, LOWPLIM_K_DATA.ktp2.values, atol=0.01)
    assert numpy.allclose(lowp_lim_ktps3, LOWPLIM_K_DATA.ktp3.values, atol=0.01)
    assert numpy.allclose(lowp_lim_ktps4, LOWPLIM_K_DATA.ktp4.values, atol=0.01)
