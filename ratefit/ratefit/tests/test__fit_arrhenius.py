"""
Test an Arrhenius fit of T, k(T,P) rates to a single Arrhenius function
where the fits are performed using Python and dsarrfit
"""

import os
import tempfile
import numpy
import ratefit


# Set Reference Data (K(T)s from 500-2500 K, 100 K Steps)
RATEKS = numpy.array([
    1.378150410638410e+08,
    1.297836416161970e+09,
    6.769080072288560e+09,
    2.425098083508950e+10,
    6.735608561972120e+10,
    1.560888894758010e+11,
    3.164069841563090e+11,
    5.792159567474590e+11,
    9.791388652530330e+11,
    1.553294914048580e+12,
    2.340215195297280e+12,
    3.378946826539680e+12,
    4.708351078718470e+12,
    6.366577818492750e+12,
    8.390689299332350e+12,
    1.081640479955560e+13,
    1.367793986577610e+13,
    1.700791768884750e+13,
    2.083733419920020e+13,
    2.519556223651740e+13,
    3.011038339418150e+13
])

# Temperature Range and Parameter Setting
TEMPS = numpy.arange(500.0, 2600.0, 100.0)
T_REF = 1.0
NA = 6.0221409e23
NA_INV = (1.0 / NA)

# Set and build paths to run the dsarrfit code
DSARRFIT_PATH_BASE = tempfile.mkdtemp()
DSARRFIT_PATH1 = os.path.join(DSARRFIT_PATH_BASE, 'SGL')
DSARRFIT_PATH2 = os.path.join(DSARRFIT_PATH_BASE, 'DBL')
os.mkdir(DSARRFIT_PATH1)
os.mkdir(DSARRFIT_PATH2)
print(DSARRFIT_PATH1)
print(DSARRFIT_PATH2)


def test__arrhenius_fit_dsarrfit():
    """ test ratefit.fit.arrhenius.single using dsarrfit code
    """

    # Run a single Arrhenius fit
    sgl_fit_params = ratefit.fit.arrhenius.single(
        TEMPS, RATEKS, T_REF, 'dsarrfit', dsarrfit_path=DSARRFIT_PATH1)

    ref_sgl_fit_params = (26.0852, 3.76371, 7353.251918248664)

    assert numpy.allclose(sgl_fit_params, ref_sgl_fit_params)

    # Run a double Arrhenius fit
    dbl_fit_params = ratefit.fit.arrhenius.double(
        TEMPS, RATEKS, T_REF, 'dsarrfit', dsarrfit_path=DSARRFIT_PATH2)

    ref_dbl_fit_params = (75450400.0, 1.87695, 11517.021129336208,
                          673789.0, 2.40717, 10460.504113187226)

    assert numpy.allclose(dbl_fit_params, ref_dbl_fit_params)

    # Run a fit that should die to test for error setting (should be double I think)
    # runner


def test__arrhenius_fit_python():
    """ test ratefit.fit.arrhenius.double
    """

    # Run a single Arrhenius fit to build a guess for the double fit
    sgl_fit_params = ratefit.fit.arrhenius.single(
        TEMPS, RATEKS, T_REF, 'python')

    ref_sgl_fit_params = (23339999.999677252, 2.0840000000019057, 11103.999999999936)

    assert numpy.allclose(sgl_fit_params, ref_sgl_fit_params)

    # Use Single Fit Params to Run the Double Fitter
    dbl_fit_params = ratefit.fit.arrhenius.double(
        TEMPS, RATEKS, T_REF, 'python', arr1_guess=sgl_fit_params)

    ref_dbl_fit_params = (18070230.266769655, 2.0840017657840666, 11103.99642254344,
                          5269769.7413749155, 2.083993944680211, 11104.012268002933)

    assert numpy.allclose(dbl_fit_params, ref_dbl_fit_params)


if __name__ == '__main__':
    test__arrhenius_fit_dsarrfit()
    test__arrhenius_fit_python()
