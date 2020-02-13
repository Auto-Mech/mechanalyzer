"""
Test an Arrhenius fit of T, k(T,P) rates to a single Arrhenius function
where the fits are performed using Python
"""

import ratefit

NA = 6.0221409e23
NA_INV = (1.0 / NA)

# Obtain list of temperatures and rate constants from initial pair list
PAIRS = [
    ['300', '***'],
    ['400', '***'],
    ['500', '-30.0'],
    ['600', '93.5678'],
    ['700', '2287.02'],
    ['800', '27103.4'],
    ['900', '193348'],
    ['1000', '955781'],
    ['1100', '3.60E+06'],
    ['1200', '1.10E+07'],
    ['1300', '2.85E+07'],
    ['1400', '6.51E+07'],
    ['1500', '1.34E+08'],
    ['1600', '2.52E+08'],
    ['1700', '4.43E+08'],
    ['1800', '7.34E+08'],
    ['1900', '1.16E+09'],
    ['2000', '1.74E+09']
    # ['2100', '2.53E+09'],
    # ['2200', '3.56E+09'],
    # ['2300', '4.86E+09'],
    # ['2400', '6.48E+09'],
    # ['2500', '8.45E+09'],
    # ['2600', '1.08E+10'],
    # ['2700', '1.36E+10'],
    # ['2800', '1.68E+10'],
    # ['2900', '2.06E+10'],
    # ['3000', '2.48E+10']
]
TEMPS = [pair[0] for pair in PAIRS]
RATE_CONSTANTS = [pair[1] for pair in PAIRS]
BIMOL = False

# Set the T0 value in the (T/T0)^n term in the Arrhenius expr.
T_REF = 1.0


def test__troe_fit():
    """ test ratefit.fit.arrhenius.single
    """

    # Filter the temperatures and rate constants to get valid values
    # k > 0 and k != *** and tmin <= T <= tmax
    tmin = None
    tmax = None
    temps, calc_ks = ratefit.fit.get_valid_tk(
        TEMPS, RATE_CONSTANTS, BIMOL, tmin=tmin, tmax=tmax)
    print('Fit Range =', [tmin, tmax])
    print('Calc ks =', calc_ks)

    # Print header for python fitting
    print('\n\nTroe fits')

    # Run the TroeFit program to fit the rate constants
    fit_params = ratefit.fit.arrhenius.single(
        temps, calc_ks, T_REF, 'python')
    print('\nFit Parameters:')
    print('A =', fit_params[0])
    print('n =', fit_params[1])
    print('Ea =', fit_params[2])

    # Calculate fitted rate constants using the fitted parameters
    fit_ks = ratefit.fxns.single_arrhenius(
        fit_params[0], fit_params[1], fit_params[2],
        T_REF, temps)

    # Print the fitted rate constants and errors
    print('\nComparison of Calculated vs. Fitted Rate Constants:')
    print('Temp (K)  Calc ks      Fit ks')
    for i, _ in enumerate(temps):
        print('{0:6.1f}    {1:1.5E}  {2:1.5E}'.format(
            temps[i], calc_ks[i], fit_ks[i]))

    # Calculate the sum-of-square errors and mean-average-errors
    sse, mean_avg_err, max_avg_err = ratefit.err.calc_sse_and_mae(
        calc_ks, fit_ks)
    print('\nSSE =', sse)
    print('Mean Avg. Err = ', mean_avg_err)
    print('Max Avg. Err = ', max_avg_err)

    # Print header for dsarrfit fitting
    print('\n\n\nSingle Arrhenius Fit with dsarrfit:')

    # Run a single Arrhenius fit
    fit_params2 = ratefit.fit.arrhenius.single(
        temps, calc_ks, T_REF, 'dsarrfit', dsarrfit_path='.',
        a_conv_factor=1.0)
    print('\nFit Parameters:')
    print('A =', fit_params2[0])
    print('n =', fit_params2[1])
    print('Ea =', fit_params2[2])

    # Calculate fitted rate constants using the fitted parameters
    fit_ks2 = ratefit.fxns.single_arrhenius(
        fit_params2[0], fit_params2[1], fit_params2[2],
        T_REF, temps)

    # Print the fitted rate constants and errors
    print('\nComparison of Calculated vs. Fitted Rate Constants:')
    print('Temp (K)  Calc ks      Fit ks')
    for i, _ in enumerate(temps):
        print('{0:6.1f}    {1:1.5E}  {2:1.5E}'.format(
            temps[i], calc_ks[i], fit_ks2[i]))

    # Calculate the sum-of-square errors and mean-average-errors
    sse2, mean_avg_err2, max_avg_err2 = ratefit.err.calc_sse_and_mae(
        calc_ks, fit_ks2)
    print('\nSSE =', sse2)
    print('Mean Avg. Err = ', mean_avg_err2)
    print('Max Avg. Err = ', max_avg_err2)


if __name__ == '__main__':
    test__single_troe_fit()
