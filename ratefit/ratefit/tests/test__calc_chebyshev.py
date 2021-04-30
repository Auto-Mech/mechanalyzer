"""
Test the ratefit rate constant calculators
"""

import numpy
import ratefit


# Pressure and Temperature Range and Parameter Setting
TMIN, TMAX = (500.0, 2000.0)
PMIN, PMAX = (0.01, 20.0)
TEMPS = numpy.arange(TMIN, TMAX+100.0, 100.0)
PRESSURES = (0.1, 1.0, 5.0, 10.0, 20.0)

ALPHA = numpy.array([
    [8.684e+00, 7.500e-01, -7.486e-02, 1.879e-15],
    [-2.159e-01, 9.899e-02, 2.292e-02, 2.929e-17],
    [-1.557e-15, -3.331e-16, 3.324e-17, -8.346e-31],
    [2.159e-01, -9.899e-02, -2.292e-02, -2.929e-17],
    [-2.684e+00, -7.500e-01, 7.486e-02, -1.879e-15],
    [2.159e-01, -9.899e-02, -2.292e-02, -2.929e-17]
])

# Set Reference Data for Comparison
REF_CHEB_KTPS = {
    0.1: numpy.array(
        [5.36149469e+05, 3.37524014e+10, 1.08278185e+07, 1.00000000e+06,
         4.67760554e+06, 7.50837664e+07, 1.00024158e+09, 6.34096052e+09,
         1.68471688e+10, 1.98701446e+10, 1.17622251e+10, 3.99642994e+09,
         8.81901757e+08, 1.40478161e+08, 1.76174826e+07, 1.86515152e+06]),
    1.0: numpy.array(
        [6.08409071e+05, 1.81382797e+11, 1.72942896e+07, 1.00000000e+06,
         6.16617966e+06, 1.81248130e+08, 4.55199960e+09, 4.78336686e+10,
         1.73639423e+11, 2.27668672e+11, 1.23451059e+11, 3.24126562e+10,
         4.78849631e+09, 4.54016299e+08, 3.08476611e+07, 1.64363098e+06]),
    5.0: numpy.array(
        [6.95918671e+05, 4.27401045e+11, 2.25326592e+07, 1.00000000e+06,
         7.15906367e+06, 3.00666831e+08, 1.11517627e+10, 1.61630221e+11,
         7.21445431e+11, 1.02306394e+12, 5.29971627e+11, 1.18905026e+11,
         1.36387243e+10, 9.26863017e+08, 4.22933849e+07, 1.43694952e+06]),
    10.0: numpy.array(
        [7.46035340e+05, 5.70314722e+11, 2.48540825e+07, 1.00000000e+06,
         7.55017560e+06, 3.63631096e+08, 1.57390262e+10, 2.59902792e+11,
         1.26441957e+12, 1.85667350e+12, 9.46804806e+11, 1.99680802e+11,
         2.06818218e+10, 1.22687767e+09, 4.75261892e+07, 1.34041907e+06]),
    20.0: numpy.array(
        [8.05396986e+05, 7.24936424e+11, 2.71535504e+07, 1.00000000e+06,
         7.90957932e+06, 4.32470186e+08, 2.16670298e+10, 4.05682226e+11,
         2.14736793e+12, 3.26728174e+12, 1.64407784e+12, 3.27077461e+11,
         3.07183062e+10, 1.59787535e+09, 5.27907769e+07, 1.24162372e+06])
}


def test__chebyshev():
    """ test ratefit.fxns.chebyshev
    """

    cheb_ktps = ratefit.calc.chebyshev(
        ALPHA, TMIN, TMAX, PMIN, PMAX, TEMPS, PRESSURES)

    assert numpy.allclose(tuple(cheb_ktps.keys()), PRESSURES)

    assert numpy.allclose(cheb_ktps[0.1], REF_CHEB_KTPS[0.1], atol=0.01)
    assert numpy.allclose(cheb_ktps[1.0], REF_CHEB_KTPS[1.0], atol=0.01)
    assert numpy.allclose(cheb_ktps[5.0], REF_CHEB_KTPS[5.0], atol=0.01)
    assert numpy.allclose(cheb_ktps[10.0], REF_CHEB_KTPS[10.0], atol=0.01)
    assert numpy.allclose(cheb_ktps[20.0], REF_CHEB_KTPS[20.0], atol=0.01)


if __name__ == '__main__':
    test__chebyshev()
