""" interface to troefit code by sjk
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
    """ write the troefit input file
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
    template_file_name = 'troefit.mako'
    template_file_path = os.path.join(SRC_PATH, template_file_name)

    # Build global section string
    troefit_str = Template(filename=template_file_path).render(**fit_keys)

    return troefit_str


def run_troefit(path):
    """ run arrfit code
    """
    start_path = os.getcwd()
    os.chdir(path)
    # Set the full path to the troefit executable
    exe_path = os.path.join(SRC_PATH, 'troefit', 'troefit.x')
    # Run the executable
    subprocess.check_call([exe_path])
    os.chdir(start_path)


def read_params(output_string, conv_factor=1.000):
    """ obtain information from the arrfit output
    """

    # Loop over the lines and find the resulting fit params line
    lines = output_string.splitlines()
    lines.reverse()
    kinf_str, k0_str, ptroe_str = '', '', ''
    for line in lines:
        if line.startswith(' results for iteration'):
            kinf_str = lines[lines.index(line)-3]
            k0_str = lines[lines.index(line)-4]
            ptroe_str = lines[lines.index(line)-5]
            break

    # Assess status of fits (single always assumed True for now)
    fit_success = bool(kinf_str and k0_str and ptroe_str)

    # Grab the fitting parameters
    # Multiple A by given conversion factor and Ea/R term by R to get Ea
    if fit_success:
        high_params = [float(param) for param in kinf_str.split()]
        low_params = [float(param) for param in k0_str.split()]
        troe_params = [float(param) for param in ptroe_str.split()]
        temp_params = troe_params[:3]
        fcent = troe_params[3]
        fit_params = [high_params, low_params, temp_params, fcent]
    else:
        fit_params = []

    return fit_params


if __name__ == '__main__':
    with open('ex/0/arrfit.out', 'r') as f:
        string = f.read()
    print(read_params(string, conv_factor=1.000))
