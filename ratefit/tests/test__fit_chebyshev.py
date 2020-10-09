"""
Test an Arrhenius fit of T, k(T,P) rates to a single Arrhenius function
where the fits are performed using Python and dsarrfit
"""

import os
import tempfile
import numpy
import pandas
import ratefit


# def _read_csv(filename):
#     """ read csv values from arrhenius.csv
#     """
#     csv_file = open(filename, 'r')
#     data = pandas.read_csv(csv_file, comment='!', quotechar="'")
#     csv_file.close()
# 
#     return data
# 
# 
# # Set path to data files
# PATH = os.path.dirname(os.path.realpath(__file__))
# DATA_PATH = os.path.join(PATH, 'data')
# RATEK_FILE_NAME = 'rateks.csv'
# 
# # Read csv file for data
# RATEK_DATA = _read_csv(
#     os.path.join(DATA_PATH, RATEK_FILE_NAME))
# 
# # Set the rate constant data
# TEMPS = RATEK_DATA.TEMPS.values
# RATEKS = RATEK_DATA.rateks.values

# Set NA and the T0 value in the (T/T0)^n term in the Arrhenius expr.
NA = 6.0221409e23
NA_INV = (1.0 / NA)
TDEG = 6
PDEG = 4

TMIN, TMAX = (300.000, 2200.000)
PMIN, PMAX = (0.010, 98.702)
TEMPS = numpy.arange(TMIN, TMAX, 100.0)
PRESSURES = (0.1, 1.0, 10.0, 20.0)

RATEKS = {
 0.1:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
numpy.array([5.36149469e+05, 9.97979437e+08, 1.22231205e+06, 3.16691208e+06,
       6.07522933e+07, 8.62460463e+08, 5.27826470e+09, 1.43145909e+10,
       2.03018305e+10, 1.78439264e+10, 1.11745217e+10, 5.54400107e+09,
       2.35585025e+09, 9.07197360e+08, 3.29643819e+08, 1.16327304e+08,
       4.06873664e+07, 1.43082405e+07, 5.10931923e+06, 1.86515152e+06])),
 1:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
 numpy.array([6.08409071e+05, 3.29464630e+09, 1.29188256e+06, 3.86300416e+06,
       1.39653166e+08, 3.77763148e+09, 3.77506955e+10, 1.39199261e+11,
       2.27880646e+11, 2.02594174e+11, 1.15986346e+11, 4.88073967e+10,
       1.66646999e+10, 4.96423316e+09, 1.36002040e+09, 3.55962747e+08,
       9.14721736e+07, 2.35314021e+07, 6.14383793e+06, 1.64363098e+06])),
 10:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
 numpy.array([7.46035340e+05, 7.75068486e+09, 1.35762901e+06, 4.45236841e+06,
       2.68933325e+08, 1.26310114e+10, 1.95384632e+11, 9.58250266e+11,
       1.80984024e+12, 1.64895566e+12, 8.81539127e+11, 3.22821335e+11,
       9.10634707e+10, 2.15903869e+10, 4.58714006e+09, 9.15161978e+08,
       1.77362652e+08, 3.42179850e+07, 6.68748385e+06, 1.34041907e+06])),
 20:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
 numpy.array([8.05396986e+05, 9.38351170e+09, 1.37652636e+06, 4.59549763e+06,
       3.16415767e+08, 1.72298258e+10, 3.00821082e+11, 1.60072866e+12,
       3.15592329e+12, 2.90411737e+12, 1.52719941e+12, 5.39053260e+11,
       1.44349518e+11, 3.21298356e+10, 6.35842272e+09, 1.17527171e+09,
       2.10303924e+08, 3.73883339e+07, 6.72786756e+06, 1.24162372e+06])),
 40:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
 numpy.array([8.75609139e+05, 1.10168031e+10, 1.39496387e+06, 4.71892054e+06,
       3.66356183e+08, 2.29349334e+10, 4.49776775e+11, 2.59156548e+12,
       5.33332535e+12, 4.96272211e+12, 2.57208449e+12, 8.77096804e+11,
       2.23525873e+11, 4.68293404e+10, 8.65420688e+09, 1.48567949e+09,
       2.46042278e+08, 4.03992924e+07, 6.70768270e+06, 1.14206209e+06])),
 80:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
  numpy.array([9.58650400e+05, 1.25432664e+10, 1.41291616e+06, 4.82082588e+06,
       4.17427299e+08, 2.97911950e+10, 6.53065400e+11, 4.06642107e+12,
       8.73484414e+12, 8.22860771e+12, 4.21124046e+12, 1.39062420e+12,
       3.38128348e+11, 6.68480958e+10, 1.15657986e+10, 1.84866620e+09,
       2.84020221e+08, 4.31685406e+07, 6.62747314e+06, 1.04313314e+06])),
 98.702:  (numpy.arange(TMIN, TMAX+100.0, 100.0),
    numpy.array([9.86707971e+05, 1.29675352e+10, 1.41825752e+06, 4.84722339e+06,
       4.32895027e+08, 3.20934576e+10, 7.26992115e+11, 4.63260659e+12,
       1.00809999e+13, 9.53449737e+12, 4.86278376e+12, 1.59093559e+12,
       3.81553306e+11, 7.41568869e+10, 1.25829404e+10, 1.96914411e+09,
       2.95864770e+08, 4.39478436e+07, 6.59159521e+06, 1.01347109e+06]))}



def test__chebyshev_fit():
    """ test ratefit.fit.troe.std
    """

    # Run a Chebyshev fit
    alpha, trange, prange = ratefit.fit.chebyshev.kfit(
        TEMPS, RATEKS, tdeg=TDEG, pdeg=PDEG)
    print('alpha\n', alpha)
    # assert numpy.allclose(fit_params, ref_fit_params)

    # # Calculate fitted rate constants using the fitted parameters
    fit_ktps = ratefit.calc.chebyshev(
        alpha, TMIN, TMAX, PMIN, PMAX, TEMPS, PRESSURES)
    # print('fit ks', fit_ks)
    # assert numpy.allclose(fit_ks, ref_fit_ks)

    # # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_errs, max_avg_errs = [], []
    for pressure in PRESSURES:
        rate_kts = RATEKS[pressure][1]
        fit_kts = numpy.array(fit_ktps[pressure])
        mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
            rate_kts, fit_kts)
        mean_avg_errs.append(mean_avg_err)
        max_avg_errs.append(max_avg_err)
        print('rate kts')
        print(rate_kts)
        print('fit kts')
        print(fit_kts)
        print('errs')
        print(mean_avg_err, max_avg_err)

    # ref_mean_avg_err = 2.504321595549729
    # ref_max_avg_err = 7.827700474963309

    # assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    # assert numpy.allclose(max_avg_err, ref_max_avg_err)


if __name__ == '__main__':
    test__chebyshev_fit()
