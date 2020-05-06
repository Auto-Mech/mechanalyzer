""" interface to troefit code by sjk
"""

import os
import subprocess
from mako.template import Template


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K

# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))


def write_input(kpt_dct, troe_param_fit_lst,
                highp_a=8.1e-11, highp_n=-0.01, highp_ea=1000.0,
                lowp_a=8.1e-11, lowp_n=-0.01, lowp_ea=1000.0,
                alpha=0.19, ts1=590, ts2=1.e6, ts3=6.e4,
                fit_tol1=1.0e-8, fit_tol2=1.0e-8):
    """ write the troefit input file
    """

    # Write the fitting parameters at the top of the file
    fit_tol_str = '{0:<12.10f} {1:<12.10f}'.format(fit_tol1, fit_tol2)

    # Set the indices for what is fitted
    fit_idxs = ['1', '1', '1', '1', '1', '1', None, None, None, None]
    fit_idxs[6] = '1' if 'alpha' in troe_param_fit_lst else '0'
    fit_idxs[7] = '1' if 'ts1' in troe_param_fit_lst else '0'
    fit_idxs[8] = '1' if 'ts2' in troe_param_fit_lst else '0'
    fit_idxs[9] = '1' if 'ts3' in troe_param_fit_lst else '0'
    fit_idx_str = ' '.join(fit_idxs)
    fit_option_idx = 1

    # Write the parameters strings
    highp_params_str = '{0:>8.5E} {1:>8.5E} {2:>8.5E}'.format(
        highp_a, highp_n, highp_ea)
    lowp_params_str = '{0:>8.5E} {1:>8.5E} {2:>8.5E}'.format(
        lowp_a, lowp_n, lowp_ea)
    troe_params_str = '{0:>8.5E} {1:>8.5E} {2:>8.5E} {3:>8.5E}'.format(
        alpha, ts1, ts2, ts3)

    # Write the rate constants
    num_tk = len(kpt_dct)
    kpt_str = ''
    for temp, pk_arr in kpt_dct.items():
        [pressures, rate_constants] = pk_arr
        kpt_str += '{0:<8.2f}{1:<4d}\n'.format(temp, len(pressures))
        for rate, pressure in zip(rate_constants, pressures):
            kpt_str += '{0:<14.5E}{1:<14.5f}\n'.format(rate, pressure)

    # Build the fill value dictionary
    fit_keys = {
        'fit_tol_str': fit_tol_str,
        'fit_idx_str': fit_idx_str,
        'fit_option_idx': fit_option_idx,
        'highp_params_str': highp_params_str,
        'lowp_params_str': lowp_params_str,
        'troe_params_str': troe_params_str,
        'num_tk': num_tk,
        'kpt_str': kpt_str
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

    # Go to path
    start_path = os.getcwd()
    os.chdir(path)

    # Run the executable
    exe_cmd = 'troefit.x'
    subprocess.check_call([exe_cmd])

    # Return to starting dir
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
    _ = conv_factor
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
