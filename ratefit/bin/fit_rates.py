"""
  Script to fit rates for a standalone
"""

import os
import sys
import json
import ratefit

PES_FORMULA = 'PES'

mess_out_name = sys.argv[0]
label_dct_inp_name = sys.argv[1]

# Build a label dictionary for bulding the rates
if label_dct_inp_name:
    # Read a label dictionary from a file
    label_dct_file = os.path.join(os.getcwd(), label_dct_inp_name)
    with open(label_dct_inp_name) as fobj:
        data = fobj.read()
    label_dct = json.loads(data)
else:
    label_dct = {}

# Fit the rates
ratefit.fit.fit_ktp_dct(
    mess_path, pes_formula,
    # pdep_dct={},
    # arrfit_dct={},
    # troe_dct={},
    label_dct=None,
    fit_temps=None, fit_pressure=None,
    fit_tunit='K', fit_punit='atm')
