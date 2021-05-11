""" test ratefit.fit._pdep
"""

from ratefit.fit import pressure_dependent_ktp_dct


KTP_DCT = {
    0.1:    ((500.0, 1000.0), (1.020, 1.050)),
    1.0:    ((300.0, 500.0, 1000.0), (1.0, 1.030, 1.060)),
    10.0:   ((500.0, 1000.0), (1.040, 2.080)),
    100.0:  ((500.0, 1300.0, 1000.0), (1.050, 1.50, 2.090)),
    'high': ((500.0, 1000.0), (1.060, 3.040))
}
TEMPS = (500.0,)
PVAL = 2.0


def test__assess_pdep():
    """ test ratefit.fit._pdep.pressure_dependent_ktp_dct
        test ratefit.fit._pdep.assess_pressure_dependence
    """

    ref_ktp_dct1 = {
        0.1: ((500.0, 1000.0), (1.02, 1.05)),
        1.0: ((300.0, 500.0, 1000.0), (1.0, 1.03, 1.06)),
        10.0: ((500.0, 1000.0), (1.04, 2.08)),
        100.0: ((500.0, 1300.0, 1000.0), (1.05, 1.5, 2.09)),
        'high': ((500.0, 1000.0), (1.06, 3.04))}
    ref_ktp_dct2 = {'high': ((300.0, 500.0, 1000.0), (1.0, 1.03, 1.06))}

    assert pressure_dependent_ktp_dct(KTP_DCT) == ref_ktp_dct1
    assert pressure_dependent_ktp_dct(KTP_DCT, temps=TEMPS) == ref_ktp_dct2
    assert pressure_dependent_ktp_dct(KTP_DCT, temps=TEMPS, pval=PVAL) is None
