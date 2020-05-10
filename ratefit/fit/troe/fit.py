""" fit rate constants to Arrhenius expressions
"""

import os
from ratefit.fit.troe import troefit_io


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def std_form(kpt_dct, troe_param_fit_lst, troefit_path,
             highp_a=8.1e-11, highp_n=-0.01, highp_ea=1000.0,
             lowp_a=8.1e-11, lowp_n=-0.01, lowp_ea=1000.0,
             alpha=0.19, ts1=590, ts2=1.e6, ts3=6.e4,
             fit_tol1=1.0e-8, fit_tol2=1.0e-8,
             a_conv_factor=1.0):
    """ call the troefit program to get the troefits
    """

    # Write the input file for the ratefit code
    ratefit_inp_str = troefit_io.write_input(
        kpt_dct, troe_param_fit_lst,
        highp_a=highp_a, highp_n=highp_n, highp_ea=highp_ea,
        lowp_a=lowp_a, lowp_n=lowp_n, lowp_ea=lowp_ea,
        alpha=alpha, ts1=ts1, ts2=ts2, ts3=ts3,
        fit_tol1=fit_tol1, fit_tol2=fit_tol2)
    troefit_inp_file = os.path.join(troefit_path, 'troefit.dat')
    with open(troefit_inp_file, 'w') as troefit_infile:
        troefit_infile.write(ratefit_inp_str)

    # Run the ratefit program
    troefit_io.run_troefit(troefit_path)

    troename = 'troefit.out'
    troefit_out_file = os.path.join(troefit_path, troename)
    with open(troefit_out_file, 'r') as troefit_outfile:
        troefit_out_str = troefit_outfile.read()

    # Parse the ratefit files for the Arrhenius fit parameters
    fit_params = troefit_io.read_params(
        troefit_out_str, a_conv_factor)

    return fit_params
