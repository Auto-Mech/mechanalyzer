""" fit rate constants to Arrhenius expressions
"""

import os
import numpy as np
from scipy.optimize import leastsq
from ratefit.fit.arrhenius import dsarrfit_io


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def single(temps, rate_constants, t_ref, method,
           a_guess=8.1e-11, n_guess=-0.01, ea_guess=2000.0,
           dsarrfit_path=None, a_conv_factor=1.00):
    """ call the single arrhenius fitter
    """

    if method == 'python' or len(rate_constants) <= 3:
        fit_params = _single_arrhenius_numpy(
            temps, rate_constants, t_ref, a_conv_factor)
    elif method == 'dsarrfit':
        assert dsarrfit_path is not None
        fit_params = _dsarrfit(
            temps, rate_constants, a_guess, n_guess, ea_guess,
            'single', dsarrfit_path, a_conv_factor)
    else:
        raise NotImplementedError

    return fit_params


def double(temps, rate_constants, t_ref, method,
           a_guess=8.1e-11, n_guess=-0.01, ea_guess=2000.0,
           dsarrfit_path=None, a_conv_factor=1.00):
    """ call the double arrhenius fitter
    """

    if len(rate_constants) <= 3:
        fit_params = _single_arrhenius_numpy(
            temps, rate_constants, t_ref, a_conv_factor)
    elif method == 'dsarrfit':
        assert dsarrfit_path is not None
        fit_params = _dsarrfit(
            temps, rate_constants, a_guess, n_guess, ea_guess,
            'double', dsarrfit_path, a_conv_factor)
    elif method == 'python':
        fit_params = _double_arrhenius_scipy(
            temps, rate_constants, t_ref, a_guess, n_guess, ea_guess)
    else:
        raise NotImplementedError

    return fit_params


def _single_arrhenius_numpy(temps, rate_constants, t_ref, a_conv_factor=1.):
    """ Fit a set of T-dependent rate constants to a single Arrhenius
        functional expression using numpy.
        :param numpy.ndarray temps: temperatures
        :param numpy.ndarray rate_constants: rate constants
        :param float t_ref: reference temperature (K)
        :param float a_conv_factor: Conversion factor for A parameter
        :return fit_params: A, n, Ea fitting parameters for function
        :rtype: list
    """

    # temps = temps[0]
    # rate_constants = rate_constants[0]

    # consider several cases depending on the number of valid rate constants
    # no k is positive, so return all zeros
    if rate_constants.size == 0:
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0

    # if num(k) > 0 is 1: set A = k
    elif rate_constants.size == 1:
        a_fit, n_fit, ea_fit = rate_constants[0], 0.0, 0.0

    # if num(k) > 0 is 2,3: fit A and Ea
    elif rate_constants.size in (2, 3):
        # Build vectors and matrices used for the fitting
        a_vec = np.ones(len(temps))
        ea_vec = (-1.0 / RC) * (1.0 / temps)
        coeff_mat = np.array([a_vec, ea_vec], dtype=np.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = np.log(rate_constants)
        # Perform the least-squares fit
        theta = np.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = np.exp(theta[0]), 0.0, theta[1]

    # if num(k) > 0 is more than 3: fit A, n, and Ea
    elif rate_constants.size > 3:
        # Build vectors and matrices used for the fitting
        a_vec = np.ones(len(temps))
        n_vec = np.log(temps / t_ref)
        ea_vec = (-1.0 / RC) * (1.0 / temps)
        coeff_mat = np.array([a_vec, n_vec, ea_vec], dtype=np.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = np.log(rate_constants)
        # Perform the least-squares fit
        theta = np.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = np.exp(theta[0]), theta[1], theta[2]

    # Pack the parameters into a list
    fit_params = [a_fit, n_fit, ea_fit]
    fit_params[0] *= a_conv_factor

    return fit_params


def _double_arrhenius_scipy(temps, rate_constants, t_ref,
                            sgl_a, sgl_n, sgl_ea):
    """ Fit a set of T-dependent rate constants to a double Arrhenius
        functional expression using scipy.
        :param numpy.ndarray temps: temperatures
        :param numpy.ndarray rate_constants: rate constants
        :param float sgl_a: seed guess value for A parameters
        :param float sgl_n: seed guess value for n parameters
        :param float sgl_ea: seed guess value for Ea parameters
        :return fit_params: A1, A2, n1, n2, Ea1 and Ea2
                            fitting parameters for function
        :rtype: list
    """

    # Build a guess vector
    guess_params = [(sgl_a / 2.0), (sgl_n + 0.1), sgl_ea,
                    (sgl_a / 2.0), (sgl_n - 0.1), sgl_ea]

    # Perform a new least-squares fit
    plsq = leastsq(_mod_arr_residuals, guess_params,
                   args=(rate_constants, temps, t_ref),
                   ftol=1.0E-9, xtol=1.0E-9, maxfev=100000)
    fit_params = list(plsq[0])

    return fit_params


def _mod_arr_residuals(guess_params, rate_constant, temp, t_ref):
    """ Subroutine computes the residual used by the nonlinear solver
        when determing a double Arrhenius fit using scipy.
        :param numpy.ndarray guess_params: Guess Parameters for function
        :param float rate_constant: Rate constant value
        :param float temp: Temperature 
        :param float t_ref: Reference temperature
        :return err: error in the residual
        :rtype: numpy.ndarray
    """

      # compute the fitted rate constant
    k_fit1 = np.exp(
        np.log(guess_params[0]) +
        guess_params[1] * np.log(temp/t_ref) -
        guess_params[2]/(RC * temp)
    )
    k_fit2 = np.exp(
        np.log(guess_params[3]) +
        guess_params[4] * np.log(temp/t_ref) -
        guess_params[5]/(RC * temp)
    )
    k_fit = k_fit1 + k_fit2

    # calculate error
    err = np.log10(rate_constant) - np.log10(k_fit)

    return err

 
def _dsarrfit(temps, rate_constants,
              a_guess, n_guess, ea_guess,
              fit_type, dsarrfit_path, a_conv_factor):
    """ Routine calls the dsarrfit code to fit a set of rate constants
        to a single or double Arrhenius functional expression
        :param numpy.ndarray temps: temperatures
        :param numpy.ndarray rate_constants: rate constants
        :param float a_guess: seed guess value for A parameters
        :param float n_guess: seed guess value for n parameters
        :param float ea_guess: seed guess value for Ea parameters
        :param string fit_type: var signaling for a single or double fit
        :param string dsarrfit_path: path to run dsarrfit
        :param float a_conv_factor: Conversion factor for A parameter
        :return fit_params: fitting parameters for function
        :rtype: list
    """
    """

    # Write the input file for the ratefit code
    ratefit_inp_str = dsarrfit_io.write_input(
        temps, rate_constants, a_guess, n_guess, ea_guess)
    dsarrfit_inp_file = os.path.join(dsarrfit_path, 'arrfit.dat')
    with open(dsarrfit_inp_file, 'w') as arrfit_infile:
        arrfit_infile.write(ratefit_inp_str)

    # Run the ratefit program
    dsarrfit_io.run_dsarrfit(dsarrfit_path)

    # Read the output of the single and double fit
    arrname = 'arrfit.out'
    dsarrfit_out_file = os.path.join(dsarrfit_path, arrname)
    with open(dsarrfit_out_file, 'r') as arrfit_outfile:
        arrfit_out_str = arrfit_outfile.read()

    # Parse the ratefit files for the Arrhenius fit parameters
    fit_params = dsarrfit_io.read_params(
        arrfit_out_str, fit_type, a_conv_factor)

    return fit_params
