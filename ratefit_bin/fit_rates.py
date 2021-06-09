#!/usr/bin/env python
"""
  Script to fit rates for a standalone
"""

import os
import json
import argparse
import ratefit


def _read_json(path, file_name):
    """ read a json file, send to pathtools/mechanalyzer parser
    """

    json_path = os.path.join(path, file_name)
    if os.path.exists(json_path):
        with open(json_path, 'r') as fobj:
            json_dct = json.load(fobj)
    else:
        json_dct = None

    return json_dct


# Set useful global variables
CWD = os.getcwd()

# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-m', '--mess', default='rate.out',
                 help='MESS ouput name (rate.out)')
PAR.add_argument('-l', '--label', default='label.inp',
                 help='label dct name (label.inp)')
# PAR.add_argument('-o', '--output', default='output.dat',
#                  help='fitted rates out (output.dat)')
PAR.add_argument('-p', '--pes', default='pes.ckin',
                 help='pes formulat (pes.ckin)')
OPTS = vars(PAR.parse_args())

# Set paths
mess_path = CWD
# mess_path = os.path.join(CWD, OPTS['mess'])
ckin_path = os.path.join(CWD, OPTS['pes'])

# Read label dct
label_dct = _read_json(CWD, OPTS['label'])

# Fit the rates
ckin_dct = ratefit.fit.fit_ktp_dct(
    mess_path=mess_path,
    inp_fit_method='arrhenius',  # set to plog (change to pdep method)
    # pdep_dct=ratefit_dct['pdep_fit'],  # use defaults for now
    # arrfit_dct=ratefit_dct['arrfit_fit'],  # use defaults for now
    # chebfit_dct=ratefit_dct['chebfit_fit'],  # use defaults for now
    # troefit_dct=ratefit_dct['troefit_fit'],  # use defaults
    label_dct=label_dct,
    # fit_temps=pes_mod_dct[pes_mod]['rate_temps'], # read from MESS out
    # fit_pressures=pes_mod_dct[pes_mod]['pressures'],  # read from MESS out
    # fit_tunit=pes_mod_dct[pes_mod]['temp_unit'], # read MESS out
    # fit_punit=pes_mod_dct[pes_mod]['pressure_unit']  # read MESS out
)

# Write the header part
# ckin_dct.update({
#     'header': writer.ckin.model_header((spc_mod,), spc_mod_dct)
# })
CKIN_STR = ''
for _, rstring in ckin_dct.items():
    CKIN_STR += rstring + '\n\n'
with open(ckin_path, 'w') as cfile:
    cfile.write(CKIN_STR)
