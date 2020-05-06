"""
Test the ratefit rate constant calculators
"""

import numpy as np
import pandas
import ratefit


TEMPS = np.array([300., 600., 900., 1200., 1500.,
                  1800., 2100., 2400., 2700., 3000.])
T_REF = 1.0
A1, N1, EA1 = 7.015e+4, 2.053, -0.3557
A2, N2, EA2 = 5.757e+12, -0.664, 0.3318

np.set_printoptions(precision=15)


def _read_csv(filename):
    """ read csv values from arrhenius.csv
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


def test__single_arrhenius():
    """ test ratefit.fxns.single_arrhenius
    """
    calc_ks = ratefit.fxns.single_arrhenius(
        A1, N1, EA1,
        T_REF, TEMPS)
    data = _read_csv('./data/arrhenius.csv')
    assert np.allclose(calc_ks, np.array(data.SingleArr), atol=0.01)


def test__double_arrhenius():
    """ test ratefit.fxns.double_arrhenius
    """
    calc_ks = ratefit.fxns.double_arrhenius(
        A1, N1, EA1,
        A2, N2, EA2,
        T_REF, TEMPS)
    data = _read_csv('./data/arrhenius.csv')
    assert np.allclose(calc_ks, np.array(data.DoubleArr), atol=0.01)


if __name__ == '__main__':
    test__single_arrhenius()
    test__double_arrhenius()
