"""
Test what happens if no rates are found for a given pressure
"""

import ratefit

NA = 6.0221409e23
NA_INV = (1.0 / NA)

# Obtain list of temperatures and rate constants from initial pair list
PAIRS = [
    ['300', '***'],
    ['400', '***'],
    ['500', '-30.0'],
]
TEMPS = [pair[0] for pair in PAIRS]
RATE_CONSTANTS = [pair[1] for pair in PAIRS]

# Set the T0 value in the (T/T0)^n term in the Arrhenius expr.
T_REF = 1.0


def test__single_arrhenius_fit():
    """ test ratefit.fit.arrhenius.single
    """

    # Filter the temperatures and rate constants to get valid values
    # k > 0 and k != *** and tmin <= T <= tmax
    temps, calc_ks = ratefit.fit.util.get_valid_tk(
        TEMPS, RATE_CONSTANTS, tmin=None, tmax=None)

    # Print header for python fitting
    if calc_ks.size == 0:
        fit_params = [1.00, 0.00, 0.00]
    print(fit_params)
    

if __name__ == '__main__':
    test__single_arrhenius_fit()

