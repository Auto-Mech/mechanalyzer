""" Script for running a comparison of thermodynamic properties between mechanisms
"""

import os
import sys
import numpy
from mechanalyzer.calculator import compare
import mechanalyzer.plotter.thermo as plot_thermo
import mechanalyzer.plotter._util as util
import mechanalyzer.parser.new_spc as spc_parser
import mechanalyzer.parser.ckin_ as ckin_parser
from ioformat import pathtools

# INPUTS
# Filenames
# Git merge complained yet again about filename changes 

OUTPUT_FILENAME = 'rmg_hr_hrf_10_25.pdf'
#OUTPUT_FILENAME = 'hr_hrf_rmg.pdf'
#OUT_TXT_FNAME = 'hr_ordered.txt'  # filename for ordered text file
OUT_TXT_FNAME = 'ordered.txt'  # filename for ordered text file

THERMO_FILENAMES = [
    'rmg.ckin',
    '1dhr_10_25_22.ckin',
    '1dhrf_10_25_22.ckin',
]
SPC_CSV_FILENAMES = [
    'rmg_new.csv',
    'anl_new.csv',
    'anl_new.csv',
]
MECH_NAMES = [
    'RMG',
    'ANL_hr',
    'ANL_hrf',
]

# Conditions
TEMPS = numpy.linspace(300, 1500, 13)

# Options
# SORT = False
SORT = True
SORT_INSTR = 'g'  # either 'h', 'cp', 's', 'g', or None
SORT_TEMP = None  # can be (1) None to sort by max difference or (2) a number
REMOVE_LONERS = False
WRITE_FILE = False  # this currently does nothing
PRINT_MISSING = True  # print spcs that are in mech(s) but not in spc.csv

# RUN FUNCTIONS
# Fix temps to include the sort_temps if it doesn't already
if SORT_TEMP is not None and SORT_TEMP not in TEMPS:
    TEMPS = numpy.append(TEMPS, SORT_TEMP)

# Get the job path and load the dcts
if len(sys.argv) > 1:
    JOB_PATH = sys.argv[1]
    print(f'The job path is {JOB_PATH}')
elif len(sys.argv) == 1:
    JOB_PATH = os.getcwd()
    print(f'No job path input; using the current directory, {JOB_PATH}')
SPC_THERM_DCTS = ckin_parser.load_spc_therm_dcts(THERMO_FILENAMES, JOB_PATH,
                                                 TEMPS)
SPC_DCTS = spc_parser.load_mech_spc_dcts(SPC_CSV_FILENAMES, JOB_PATH)

# Get the algn_spc_therm_dct
ALGN_SPC_THERM_DCT = compare.get_algn_spc_therm_dct(
    SPC_THERM_DCTS, SPC_DCTS, remove_loners=REMOVE_LONERS,
    write_file=WRITE_FILE)

# Get the combined spc_dct (used for including SMILES and InChis)
COMB_SPC_DCT = compare.get_mult_comb_mech_spc_dct(SPC_DCTS)

# Run the plotter
FIGS, SORT_ALGN_SPC_THERM_DCT = plot_thermo.build_plots(
    ALGN_SPC_THERM_DCT, spc_dct=COMB_SPC_DCT, mech_names=MECH_NAMES,
    sort=SORT, sort_instr=SORT_INSTR, sort_temp=SORT_TEMP)
util.build_pdf(FIGS, filename=OUTPUT_FILENAME, path=JOB_PATH)

# Write the ordered text file
FSTR = compare.write_ordered_str(SORT_ALGN_SPC_THERM_DCT, dct_type='therm',
    comb_mech_spc_dct=COMB_SPC_DCT, print_missing=PRINT_MISSING)
pathtools.write_file(FSTR, JOB_PATH, OUT_TXT_FNAME)
# Saving this for later...writes a comparison of species
#comp_str = compare.write_comparison(ALGN_SPC_THERM_DCT, dct_type='therm')
