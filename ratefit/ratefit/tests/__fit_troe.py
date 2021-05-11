"""
Test an Arrhenius fit of T, k(T,P) rates to a Troe function
"""

import tempfile
import numpy
import ratefit


# Set paths to run the dsarrfit code
TROEFIT_PATH = tempfile.mkdtemp()
print('TROEFIT PATH', TROEFIT_PATH)

# Rate Constants
TEMPS = numpy.arange(300.0, 3300.0, 300.0)
KTP_DCT = {
    0.1: (TEMPS,
          numpy.array(
            [9.10804237e+12, 1.32381727e+12, 4.08085808e+11, 1.72813491e+11,
             8.71750733e+10, 4.93305843e+10, 3.03035430e+10, 1.97990317e+10,
             1.35713589e+10, 9.66617670e+09])),
    2.0: (TEMPS,
          numpy.array(
            [5.49814978e+13, 1.69633215e+13, 6.60737791e+12, 3.07653037e+12,
             1.62249315e+12, 9.39965627e+11, 5.85382158e+11, 3.85776465e+11,
             2.65960784e+11, 1.90195810e+11])),
    5.0: (TEMPS,
          numpy.array(
            [8.42788325e+13, 3.22567471e+13, 1.45617858e+13, 7.20296797e+12,
             3.90278656e+12, 2.29203205e+12, 1.43833189e+12, 9.52289745e+11,
             6.58503695e+11, 4.71883493e+11])),
    10.0: (TEMPS,
           numpy.array(
            [1.14711507e+14, 4.91180769e+13, 2.54925292e+13, 1.34565982e+13,
             7.50640433e+12, 4.47219704e+12, 2.82864928e+12, 1.88160867e+12,
             1.30503129e+12, 9.37076153e+11]))
}

# PARAMS
TROE_PARAM_FIT_LST = ('ts1', 'ts2', 'ts3', 'alpha')
HIGHP_GUESS = (8.1e-11, -0.01, 1000.0)
LOWP_GUESS = (8.1e-11, -0.01, 1000.0)

ALPHA, TS1, TS2, TS3 = 0.19, 590.0, 1.0e6, 6.0e4
FIT_TOL1, FIT_TOL2 = 1.0e-8, 1.0e-8


def test__troe_fit():
    """ test ratefit.fit.troe.std
    """

    fit_params = ratefit.fit.troe.reaction(
        KTP_DCT, TROEFIT_PATH,
        troe_param_fit_lst=TROE_PARAM_FIT_LST,
        highp_guess=HIGHP_GUESS,
        lowp_guess=LOWP_GUESS,
        alpha=ALPHA, ts1=TS1, ts2=TS2, ts3=TS3,
        fit_tol1=FIT_TOL1, fit_tol2=FIT_TOL2,
        a_conv_factor=1.0)

    print('fit_params', fit_params)


if __name__ == '__main__':
    test__troe_fit()
