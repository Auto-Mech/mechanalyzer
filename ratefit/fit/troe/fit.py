""" fit rate constants to Arrhenius expressions
"""

import os
import numpy as np
from scipy.optimize import leastsq
from ratefit.fit.arrhenius import troefit_io


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def single(temps, rate_constants, t_ref, method,
           a_guess=8.1e-11, n_guess=-0.01, ea_guess=2000.0,
           troefit_path=None, a_conv_factor=1.00):
    """ call the troefit program to get the troefits
    """

    # Write the input file for the ratefit code
    ratefit_inp_str = troefit_io.write_input(
        temps, rate_constants, a_guess, n_guess, ea_guess)
    troefit_inp_file = os.path.join(troefit_path, 'arrfit.dat')
    with open(troefit_inp_file, 'w') as arrfit_infile:
        arrfit_infile.write(ratefit_inp_str)

    # Run the ratefit program
    troefit_io.run_troefit(troefit_path)

    arrname = 'arrfit.out'
    troefit_out_file = os.path.join(troefit_path, arrname)
    with open(troefit_out_file, 'r') as arrfit_outfile:
        arrfit_out_str = arrfit_outfile.read()

    # Parse the ratefit files for the Arrhenius fit parameters
    fit_params = troefit_io.read_params(
        arrfit_out_str, fit_type, a_conv_factor)

    return fit_params
