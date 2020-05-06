""" interface to dsarrfit code by sjk
"""

import os
import subprocess
from mako.template import Template


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def write_input(temps, rate_constants,
                a_guess=8.1e-11,
                n_guess=-0.01,
                ea_guess=1000.0):
    """ write the dsarrfit input file
    """

    # Format the lines with temps and rate constants
    assert len(temps) == len(rate_constants)
    num_tk = len(temps)
    tk_str = ''
    for i, _ in enumerate(temps):
        tk_str += '{0:<10.1f}{1:<1.3E}\n'.format(temps[i], rate_constants[i])

    # Build the fill value dictionary
    fit_keys = {
        'a_guess': a_guess,
        'n_guess': n_guess,
        'ea_guess': ea_guess,
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
    """ run arrfit code
    """

    # Go to path
    start_path = os.getcwd()
    os.chdir(path)

    # Run the executable
    exe_cmd = os.path.join(SRC_PATH, 'dsarrfit', 'dsarrfit.x_cfg')
    # exe_cmd = 'dsarrfit.x_cfg'
    subprocess.check_call([exe_cmd])

    # Return to starting dir
    os.chdir(start_path)


def read_params(output_string, fit, conv_factor=1.000):
    """ obtain information from the arrfit output
    """

    assert fit in ('single', 'double')

    # Loop over the lines and find the resulting fit params line
    lines = output_string.splitlines()
    lines.reverse()
    params_str = ''
    if fit == 'single':
        for line in lines:
            if line.startswith(' results for iteration'):
                params_str = lines[lines.index(line)-3]
                break
    elif fit == 'double':
        for line in lines:
            if line.startswith(' results from sum of two modified arrhenius'):
                params_str = lines[lines.index(line)-3]
                break

    # Assess status of fits (single always assumed True for now)
    single_fit_success = True
    double_fit_success = False
    if fit == 'double' and params_str:
        double_fit_success = True

    # Grab the fitting parameters
    # Multiple A by given conversion factor and Ea/R term by R to get Ea
    if fit == 'single' and single_fit_success:
        fit_params = [float(param) for param in params_str.split()]
        fit_params[0] *= conv_factor
        fit_params[2] *= RC
    elif fit == 'double' and double_fit_success:
        fit_params = [float(param) for param in params_str.split()]
        fit_params[0] *= conv_factor
        fit_params[2] *= RC
        fit_params[3] *= conv_factor
        fit_params[5] *= RC
    elif fit == 'double' and not double_fit_success:
        fit_params = []

    return fit_params
