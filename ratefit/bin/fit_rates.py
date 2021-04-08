"""
  Script to fit rates for a standalone
"""

import os
import sys
import json
import ratefit

mess_inp_name = sys.argv[1]
mess_out_name = sys.argv[2]
label_dct_inp_name = sys.argv[3]

# Maybe just read the fit methods from the file

# Build a label dictionary for bulding the rates
if label_dct_inp_name:
    # Read a label dictionary from a file
    label_dct_file = os.path.join(os.getcwd(), label_dct_inp_name)
    with open(label_dct_inp_name) as fobj:
        data = fobj.read()
    label_dct = json.loads(data)
else:
    # Build a label dictionary using labels from a MESS input file
    if mess_inp_name:
    mess_inp_file = os.path.join(os.getcwd(), mess_inp_name)
        with open(mess_inp_file) as fobj:
            data = fobj.read()
        labels = mess_io.reader.labels(mess_inp_str, read_fake=False)
        label_dct = dict(zip(label_dct, label_dct))
        
# Call the fitting functions
mess_out_file = os.path.join(os.getcwd(), mess_out_name)

ckin_str_dct = ktproutines.fit.fit_rates(
    inp_temps=pes_model_dct[pes_model]['rate_temps'],  # init to none, just use mess vals if so
    inp_pressures=pes_model_dct[pes_model]['pressures'], # ^^^
    inp_tunit=pes_model_dct[pes_model]['tunit'],  # ^^^         
    inp_punit=pes_model_dct[pes_model]['punit'],  # ^^^
    pes_formula=pes_formula,
    label_dct=label_dct,
    # es_info=es_info,  # remove, just build before call in mechdriver
    # pf_model=pf_model,  # remove, just build before call in mechdriver
    mess_path=mess_path,
    inp_fit_method=fit_method,
    pdep_fit=pdep_fit,
    arrfit_thresh=arrfit_thresh)


