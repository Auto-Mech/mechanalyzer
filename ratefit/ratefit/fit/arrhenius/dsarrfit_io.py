""" interface to dsarrfit code by sjk
"""

import os
import subprocess
from mako.template import Template


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def write_input(temps, rate_constants, fit_type='single',
                arr1_guess=(8.1e-11, -0.01, 2000.0), arr2_guess=()):
    """ Write the dsarrfit input file.

        :param list temps: temperatures (K)
        :type temps: list(float)
        :param rate_constants: T-dependent rate constants ()
        :type rate_constants: list(float)
        :param a_guess: Guess value for pre-exponential A parameter
        :type a_guess: float
        :param n_guess: Guess value for n parameter
        :type n_guess: float
        :param ea_guess: Guess value for activation energy Ea parameter
        :type ea_guess: float
        :return dsarrfit_str: string for the dsarrfit input file
        :rtype: string
    """

    # Format the string to tell the code what to optimize
    assert fit_type in ('single, double, single-double')
    if fit_type == 'single':
        arr_runs = 's'
        opt_str = '1 1 1 0 0 0'
    elif fit_type == 'double':
        arr_runs = 'd'
        opt_str = '1 1 1 1 1 1'
    else:
        arr_runs = 'sd'
        opt_str = '1 1 1 1 1 1'

    # Set the param guess lines
    param_lines = '  {} {} {}'.format(*arr1_guess)
    if arr2_guess:
        param_lines += '\n  {} {} {}'.format(*arr2_guess)
        num_param_lines = 2
    else:
        num_param_lines = 1

    # Format the lines with temps and rate constants
    assert len(temps) == len(rate_constants)
    num_tk = len(temps)
    tk_str = ''
    for i, _ in enumerate(temps):
        tk_str += '{0:<10.1f}{1:<1.3E}\n'.format(temps[i], rate_constants[i])

    # Build the fill value dictionary
    fit_keys = {
        'arr_runs': arr_runs,
        'arr_param_opt_idxs': opt_str,
        'num_param_lines': num_param_lines,
        'param_lines': param_lines,
        'num_tk': num_tk,
        'tk_lines': tk_str
        }

    # Set template name and path for the global keywords section
    template_file_name = 'dsarrfit.mako'
    template_file_path = os.path.join(SRC_PATH, template_file_name)

    # Build global section string
    dsarrfit_str = Template(filename=template_file_path).render(**fit_keys)

    return dsarrfit_str


def run_dsarrfit(path):
    """ Run the dsarrfit executable.

        :param path: Path where the dsarrfit input file exists
        :type path: str
    """

    # Go to path
    start_path = os.getcwd()
    os.chdir(path)

    # Run the executable
    exe_cmd = 'dsarrfit.x_cfg'
    # exe_cmd = 'ratefit/external/dsarrfit/debug/dsarrfit.x_cfg'
    try:
        subprocess.check_call([exe_cmd])
    except subprocess.CalledProcessError:
        print('dsarrfit failed for', path)

    # Return to starting dir
    os.chdir(start_path)


def read_params(output_str, fit_type, a_conv_factor=1.000, ea_conv_factor=RC):
    """ Parse the output of the dsarrfit code for the fitting
        parameters for either a single or double Arrhenius fit.
        :param output_str: output of dsarrfit code
        :type output_str: str
        :param fit_type:
        :type fit_type: str
        :param a_conv_factor: Conversion factor for A parameter
        :type a_conv_factor: float
    """

    assert fit_type in ('single', 'double')

    # Loop over the lines and find the resulting fit params line
    lines = output_str.splitlines()
    lines.reverse()
    params_str = ''
    if fit_type == 'single':
        for line in lines:
            if line.startswith(' results for iteration'):
                params_str = lines[lines.index(line)-3]
                break
    elif fit_type == 'double':
        for line in lines:
            if line.startswith(' results from sum of two modified arrhenius'):
                params_str = lines[lines.index(line)-3]
                break

    # Assess status of fits (single always assumed True for now)
    single_fit_success = True
    double_fit_success = False
    if fit_type == 'double' and params_str:
        double_fit_success = True

    # Grab the fitting parameters
    # Multiple A by given conversion factor and Ea/R term by R to get Ea
    if fit_type == 'single' and single_fit_success:
        fit_params = [float(param) for param in params_str.split()]
        fit_params[0] *= a_conv_factor
        fit_params[2] *= ea_conv_factor
    elif fit_type == 'double' and double_fit_success:
        fit_params = [float(param) for param in params_str.split()]
        fit_params[0] *= a_conv_factor
        fit_params[2] *= ea_conv_factor
        fit_params[3] *= a_conv_factor
        fit_params[5] *= ea_conv_factor
    elif fit_type == 'double' and not double_fit_success:
        fit_params = []

    return fit_params
