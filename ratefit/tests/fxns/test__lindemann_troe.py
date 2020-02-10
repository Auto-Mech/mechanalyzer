"""
Test the ratefit rate constant calculators for Lindemann and Troe expressions
"""

import numpy as np
import pandas
import ratefit

TEMPS = np.array([300., 600., 900., 1200., 1500.,
                  1800., 2100., 2400., 2700., 3000.])
PRESSURES = np.array([0.1, 2.0, 5.0, 10.0])
T_REF = 1.0

A_HIGH, N_HIGH, EA_HIGH = 2.000e+12, 0.900, 4.87490
A_LOW, N_LOW, EA_LOW = 2.490e24, -2.300, 4.87490
TROE_ALPHA, TROE_T3, TROE_T1, TROE_T2 = 6.0e-1, 1.0e3, 7.0, 1.7e3

np.set_printoptions(precision=15)


def _read_csv(filename):
    """ read csv values from file
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


def test__lindemann_troe():
    """ test ratefit.fxns.lindemann
        test ratefit.fxns.troe
    """
    highp_ks = ratefit.fxns.single_arrhenius(
        A_HIGH, N_HIGH, EA_HIGH,
        T_REF, TEMPS)

    lowp_ks = ratefit.fxns.single_arrhenius(
        A_LOW, N_LOW, EA_LOW,
        T_REF, TEMPS)

    lind_ktps = ratefit.fxns.lindemann(
        highp_ks, lowp_ks,
        PRESSURES, TEMPS)

    troe_ktps = ratefit.fxns.troe(
        highp_ks, lowp_ks,
        PRESSURES, TEMPS,
        TROE_ALPHA, TROE_T3, TROE_T1, TROE_T2)

    data_lind = _read_csv('./data/lindemann.csv')
    data_troe = _read_csv('./data/troe.csv')
    assert np.allclose(lind_ktps[0.1], np.array(data_lind.ktp1), atol=0.01)
    assert np.allclose(lind_ktps[2.0], np.array(data_lind.ktp2), atol=0.01)
    assert np.allclose(lind_ktps[5.0], np.array(data_lind.ktp3), atol=0.01)
    assert np.allclose(lind_ktps[10.0], np.array(data_lind.ktp4), atol=0.01)
    assert np.allclose(troe_ktps[0.1], np.array(data_troe.ktp1), atol=0.01)
    assert np.allclose(troe_ktps[2.0], np.array(data_troe.ktp2), atol=0.01)
    assert np.allclose(troe_ktps[5.0], np.array(data_troe.ktp3), atol=0.01)
    assert np.allclose(troe_ktps[10.0], np.array(data_troe.ktp4), atol=0.01)


if __name__ == '__main__':
    test__lindemann_troe()
