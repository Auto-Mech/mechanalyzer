"""
Test the ratefit rate constant calculators
"""

import numpy as np
import pandas
import ratefit

TEMPS = np.array([300., 600., 900., 1200., 1500.,
                  1800., 2100., 2400., 2700., 3000.])
PRESSURES = np.array([0.1, 0.9869, 2.0, 5.0])
T_REF = 1.0

PLOG_DCT = {
    0.0296: [2.020E+013, -1.870, 22.755],
    0.0987: [1.680E+018, -3.050, 24.323],
    0.2961: [2.500E+024, -4.630, 27.067],
    0.9869: [4.540E+026, -5.120, 27.572],
    2.9607: [7.120E+028, -5.600, 28.535],
    9.8690: [5.480E+029, -5.700, 28.899]
}

np.set_printoptions(precision=15)


def _read_csv(filename):
    """ read csv values from file
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


def test__plog():
    """ test ratefit.fxns.plog
    """
    plog_ktps = ratefit.fxns.plog(PLOG_DCT, T_REF, PRESSURES, TEMPS)
    data = _read_csv('./data/plog.csv')
    assert np.allclose(plog_ktps[0.1], np.array(data.ktp1), atol=0.01)
    assert np.allclose(plog_ktps[0.9869], np.array(data.ktp2), atol=0.01)
    assert np.allclose(plog_ktps[2.0], np.array(data.ktp3), atol=0.01)
    assert np.allclose(plog_ktps[5.0], np.array(data.ktp4), atol=0.01)


if __name__ == '__main__':
    test__plog()
