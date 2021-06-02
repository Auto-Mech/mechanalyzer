"""
Test the mechanalyzer.calculator.rates functions
"""

import numpy as np
import ratefit.ktpdct
from mechanalyzer.calculator import rates


PRESSURES = np.array([0.316, 1, 10, 100])
TEMPS = np.array([1000, 1500, 2000])
TEMPS2 = np.array([300, 500, 1000])  # only for Chebyshev test
LOW_P_RXN = (('N2O',), ('N2', 'O'), ('+M',))
HIGH_P_RXN = (('N2O',), ('N2', 'O'), ('(+M)',))
LOW_P_PARAMS = [1.04E+15, 0, 59810]
HIGH_P_PARAMS = [1.26E+12, 0, 62620]
HIGH_P_PARAMS2 = [1.18, 1E-30, 7900]


ARRHENIUS_RXN_PARAM_DCT = {
    LOW_P_RXN: ((LOW_P_PARAMS, None, None, None, None, None),)
}
ARRHENIUS_KTS = np.array([8.8273E+1, 2.0086E+6, 3.0299E+8])

LINDEMANN_RXN_PARAM_DCT = {
    HIGH_P_RXN: ((HIGH_P_PARAMS, LOW_P_PARAMS, None, None, None, None),)
}
LINDEMANN_10ATM_KTS = np.array([7.6096e-03, 1.3922e+02, 1.6753e+04])

TROE_RXN_PARAM_DCT = {
    HIGH_P_RXN: ((HIGH_P_PARAMS, LOW_P_PARAMS,
                  HIGH_P_PARAMS2, None, None, None),)
}
TROE_10ATM_KTS = np.array([7.76752254e-03, 1.37906401e+02, 1.62555727e+04])

PLOG_RXN_PARAM_DCT = {
    HIGH_P_RXN: (
        (LOW_P_PARAMS, None, None, None,
         {0.1: [1.04E+15, 0, 59810],
          1: [1.04E+16, 0, 59810],
          10: [1.04E+17, 0, 59810],
          100: [1.04E+18, 0, 59810]},
         None),)
}
PLOG_10ATM_KTS = np.array([2.78943384e+02, 6.34722925e+06, 9.57454718e+08])

CHEBYSHEV_RXN_PARAM_DCT = {
    HIGH_P_RXN: (
        (LOW_P_PARAMS,
         None,
         None,
         {'t_limits': [300.0, 2500.0],
          'p_limits': [1.0, 100.0],
          'alpha_elm': np.array([
             [1.0216E+01, -1.1083E+00, -1.9807E-01],  # from Bozzelli (2002)
             [7.8325E-01, 1.1609E+00, 1.1762E-01],    # C2H5 + O2 = C2H4 + HO2
             [-9.5707E-02, 1.0928E-01, 1.1551E-01],
             [-8.0290E-02, -1.0978E-01, 3.7074E-04],
             [-1.4830E-02, -6.0589E-02, -2.8056E-02],
             [6.9203E-03, -9.7259E-03, -1.3556E-02],
             [7.6648E-03, 6.6865E-03, -8.8244E-04]]),
          'a_units': 'moles'}, None, None),
    )
}
CHEBYSHEV_100ATM_KTS = np.array(
    [1.23879008e+07, 2.77379527e+08, 2.31421081e+10])

DUPLICATE_ARRHENIUS_RXN_PARAM_DCT = {
    LOW_P_RXN: (
        (LOW_P_PARAMS, None, None, None, None, None),
        (LOW_P_PARAMS, None, None, None, None, None)
    )
}
DUPLICATE_PLOG_RXN_PARAM_DCT = {
    HIGH_P_RXN: (
        (LOW_P_PARAMS, None, None, None,
            {0.1: [1.04E+15, 0, 59810],
             1: [1.04E+16, 0, 59810],
             10: [1.04E+17, 0, 59810],
             100: [1.04E+18, 0, 59810]}, None),
        (LOW_P_PARAMS, None, None, None,
            {0.1: [1.04E+15, 0, 59810],
             1: [1.04E+16, 0, 59810],
             10: [1.04E+17, 0, 59810],
             100: [1.04E+18, 0, 59810]}, None)
    )
}


def test__arrhenius():
    """ Test the Arrhenius calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        ARRHENIUS_RXN_PARAM_DCT, PRESSURES, TEMPS)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, LOW_P_RXN, 'high', 'rates')
    assert np.allclose(calc_rates, ARRHENIUS_KTS, rtol=1e-3)


def test__lindemann():
    """ Test the Lindemann calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        LINDEMANN_RXN_PARAM_DCT, PRESSURES, TEMPS)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, HIGH_P_RXN, 10.0, 'rates')
    assert np.allclose(calc_rates, LINDEMANN_10ATM_KTS, rtol=1e-3)


def test__troe():
    """ Test the Troe calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        TROE_RXN_PARAM_DCT, PRESSURES, TEMPS)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, HIGH_P_RXN, 10.0, 'rates')
    assert np.allclose(calc_rates, TROE_10ATM_KTS, rtol=1e-3)


def test__plog():
    """ Test the PLOG calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        PLOG_RXN_PARAM_DCT, PRESSURES, TEMPS)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, HIGH_P_RXN, 0.316, 'rates')
    assert np.allclose(calc_rates, PLOG_10ATM_KTS, rtol=1e-3)


def test__chebyshev():
    """ Test the Chebyshev calculator
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        CHEBYSHEV_RXN_PARAM_DCT, PRESSURES, TEMPS2)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, HIGH_P_RXN, 100.0, 'rates')
    assert np.allclose(calc_rates, CHEBYSHEV_100ATM_KTS, rtol=1e-3)


def test__dup_arrhenius():
    """ Test the Arrhenius calculator for a duplicate reaction
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        DUPLICATE_ARRHENIUS_RXN_PARAM_DCT, PRESSURES, TEMPS)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, LOW_P_RXN, 'high', 'rates')
    assert np.allclose(calc_rates, 2*ARRHENIUS_KTS, rtol=1e-3)


def test__dup_plog():
    """ Test the PLOG calculator for a duplicate reaction
    """
    rxn_ktp_dct = rates.eval_rxn_param_dct(
        DUPLICATE_PLOG_RXN_PARAM_DCT, PRESSURES, TEMPS)
    calc_rates = ratefit.ktpdct.read(rxn_ktp_dct, HIGH_P_RXN, 0.316, 'rates')
    assert np.allclose(calc_rates, 2*PLOG_10ATM_KTS, rtol=1e-3)


if __name__ == '__main__':
    # test__arrhenius()
    # test__lindemann()
    # test__troe()
    # test__plog()
    test__chebyshev()
    test__dup_arrhenius()
    test__dup_plog()
