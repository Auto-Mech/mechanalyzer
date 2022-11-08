""" Script for running a comparison of rate constants between mechanisms
"""

import os
import sys
import numpy
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import mechanalyzer.plotter._util as util
import mechanalyzer.parser.new_spc as spc_parser
import mechanalyzer.parser.ckin_ as ckin_parser
from ioformat import pathtools

# INPUTS
# Filenames
MECH_FILES = [
    # 'rates_noabs_11_8_22.ckin',
    # 'allrates_11_24_20.ckin',
    'rates_shortb.ckin',
    'rmg_short.ckin',
]
THERM_FILES = [
    # '1dhr_10_27_22.ckin',
    '1dhr_10_27_22.ckin',
    'rmg.ckin',
]
CSV_FILES = [
    # 'anl_new.csv',
    'anl_new.csv',
    'rmg_new.csv',
]
out_fname = 'ratesb.pdf'
mech_nicknames = [
    'ANL',
    'RMG',
]
out_txt_fname = 'ordering.txt'  # filename for the ordered text file

# Conditions
TEMPS_LST = [numpy.linspace(500, 1500, 11)]
pressures = [1, 10, 100]

# Options
SORT_METHOD = 'ratios'  # either 'ratios' or None
rev_rates = True
remove_loners = True
write_file = False


# RUN FUNCTIONS; DON'T CHANGE THIS

# Get the job path and load the dcts
if len(sys.argv) > 1:
    JOB_PATH = sys.argv[1]
    print(f'The job path is {JOB_PATH}')
elif len(sys.argv) == 1:
    JOB_PATH = os.getcwd()
    print(f'No job path input; using the current directory, {JOB_PATH}')
rxn_ktp_dcts = ckin_parser.load_rxn_ktp_dcts(
    MECH_FILES, JOB_PATH, TEMPS_LST, pressures)
spc_therm_dcts = ckin_parser.load_spc_therm_dcts(
    THERM_FILES, JOB_PATH, TEMPS_LST[0])  # NOTE: taking first entry
spc_dcts = spc_parser.load_mech_spc_dcts(CSV_FILES, JOB_PATH)

# Get the algn_rxn_ktp_dct
TEMPS = TEMPS_LST[0]  # function receives a single Numpy array of temps
algn_rxn_ktp_dct = compare.get_algn_rxn_ktp_dct(
    rxn_ktp_dcts, spc_therm_dcts, spc_dcts, TEMPS, rev_rates=rev_rates,
    remove_loners=remove_loners, write_file=write_file)

# Run the plotter
figs = plot_rates.build_plots(
    algn_rxn_ktp_dct,
    mech_names=mech_nicknames,
    ratio_sort=bool(SORT_METHOD == 'ratios'))
util.build_pdf(figs, filename=out_fname, path=JOB_PATH)

# Write the ordered text file
FSTR = compare.write_ordered_str(algn_rxn_ktp_dct, dct_type='rxn')
pathtools.write_file(FSTR, JOB_PATH, out_txt_fname)

