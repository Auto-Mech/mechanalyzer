"""
  Fit the rate constants read from the MESS output to
  Arrhenius expressions
"""


import numpy
import ratefit
import chemkin_io
from ratefit.fit import arrhenius
from ratefit.fit._util import set_a_conversion_factor


# def pes(ktp_dct, inp_temps, reaction, mess_path,
def pes(ktp_dct, reaction, mess_path,
        t_ref=1.0, tdeg=6, pdeg=4,
        fit_tolerance=20.0):
    """ Read the rates for each channel and perform the fits
    """

    # Initialize the a conversion factor
    a_conv_factor = set_a_conversion_factor(reaction)

    # Obtain the fit paramts for the 1-atm rate constants, if available
    if 1 in ktp_dct.keys():
        a_conv_factor = ''
        [temps, rate_constants] = ktp_dct[1]
        one_atm_params = arrhenius.single(
            temps, rate_constants, t_ref, 'python',
            dsarrfit_path=mess_path, a_conv_factor=a_conv_factor)
    else:
        one_atm_params = [1.0, 0.0, 0.0]

    # If there are rates for more than one pressure
    pressures = tuple(pressure for pressure in ktp_dct.keys()
                      if pressure != 'high')

    # check existence of rates at all conditions
    num_kts = []
    for pressure in pressures:
        rate_kts = ktp_dct[pressure][1] * a_conv_factor
        num_kts.append(len(rate_kts))
    fit_viable = True
    if len(set(num_kts)) != 1:
        print(
            'Different number of k(T) values at different pressures...')
        fit_viable = False

    if fit_viable:

        fit_temps = list(set(inp_temps) & set(ktp_dct[pressures[0]][0]))
        fit_temps.sort()
        # print('fit_temps', fit_temps)

        # Fit rate constants to Chebyshev polynomial
        alpha, trange, prange = ratefit.fit.chebyshev.reaction(
            fit_temps, ktp_dct, tdeg=tdeg, pdeg=pdeg,
            a_conv_factor=a_conv_factor)
        tmin, tmax = trange
        pmin, pmax = prange

        # Calculate the fitted rate constants
        fit_ktps = ratefit.calc.chebyshev(
            alpha, tmin, tmax, pmin, pmax, fit_temps, pressures)

        # Calculate errors
        err_dct, temp_dct = {}, {}
        # ioprinter.debug_message('reaction in rate fit test:', reaction)
        for pressure in pressures:
            rate_kts = ktp_dct[pressure][1] * a_conv_factor
            fit_kts = numpy.array(fit_ktps[pressure])
            mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
                rate_kts, fit_kts)

            # Add to the dct and lst
            err_dct[pressure] = [mean_avg_err, max_avg_err]

        # look at err and same num k(T) at each P
        if max((vals[1] for vals in err_dct.values())) > fit_tolerance:
            print(
                'Errors from Chebyshev fit too large (see string)...')
            fit_viable = False

        # Write the Chemkin strings
        chemkin_str = chemkin_io.writer.reaction.chebyshev(
            reaction, one_atm_params, alpha, tmin, tmax, pmin, pmax)
        chemkin_str += '\n'
        chemkin_str += chemkin_io.writer.reaction.fit_info(
            pressures, temp_dct, err_dct)

    if not fit_viable:
        # Print message and reset string to empty to trigger Arrhenius
        print('Chemkin string from Chebyshev fit')
        chemkin_str = ''

    return chemkin_str
