""" fit rate constants to Arrhenius expressions
"""

import os
from ratefit.fit.troe import troefit_io


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def reaction(kpt_dct, troe_param_fit_lst, troefit_path,
             highp_a=8.1e-11, highp_n=-0.01, highp_ea=1000.0,
             lowp_a=8.1e-11, lowp_n=-0.01, lowp_ea=1000.0,
             alpha=0.19, ts1=590.0, ts2=1.0e6, ts3=6.0e4,
             fit_tol1=1.0e-8, fit_tol2=1.0e-8,
             a_conv_factor=1.0):
    """ Fits a set of T,P-dependent rate constants [k(T,P)]s to a
        Troe functional expression using troefit fitting code.

        :param kpt_dct: k(T,P) values
        :type kpt_dct: dict[pressure:temps]
        :param troe_param_fit_lst: Troe parameters to include in fitting
        :param troe_param_fit_lst: list(str)
        :param float highp_a: seed guess value for A parameters at high-P
        :type highp_a: float
        :param float highp_n: seed guess value for n parameters at high-P
        :type highp_n: float
        :param float highp_ea: seed guess value for Ea parameters at high-P
        :type highp_ea: float
        :param float lowp_a: seed guess value for A parameters at low-P
        :type lowp_a: float
        :param float lowp_n: seed guess value for n parameters at low-P
        :type lowp_n: float
        :param float lowp_ea: seed guess value for Ea parameters at low-P
        :type lowp_ea: float
        :param alpha: Troe alpha parameter
        :type alpha: float
        :param ts3: Troe T3 parameter
        :type ts3: float
        :param ts1: Troe T1 parameter
        :type ts1: float
        :param ts2: Troe T2 parameter
        :type ts2: float
        :param fit_tol1: tolerance for fitting ???
        :type fit_tol1: float
        :param fit_tol2: tolerance for fitting ???
        :type fit_tol2: float
        :param troe_path: path to run troe
        :type troe_path: str
        :param a_conv_factor: Conversion factor for A parameter
        :type a_conv_factor: float
        :return fit_params: fitting parameters for function
        :rtype: list
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
