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

# INPUTS
# Filenames
THERMO_FILENAMES = [
    'anl.therm',
    'rmg.ckin',
]
SPC_CSV_FILENAMES = [
    'anl.csv',
    'rmg.csv',
]
MECH_NAMES = [
    'ANL',
    'RMG',
]
OUTPUT_FILENAME = 'anl_rmg_therm_g.pdf'

# Conditions
TEMPS = numpy.linspace(300, 1500, 13)

# Options
# SORT = False
SORT = True
SORT_INSTR = 'g'  # either 'h', 'cp', 's', 'g', or None
SORT_TEMP = None  # can be (1) None to sort by max difference or (2) a number
REMOVE_LONERS = True
WRITE_FILE = False  # this currently does nothing

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
FIGS = plot_thermo.build_plots(
    ALGN_SPC_THERM_DCT, spc_dct=COMB_SPC_DCT, mech_names=MECH_NAMES,
    sort=SORT, sort_instr=SORT_INSTR, sort_temp=SORT_TEMP)
util.build_pdf(FIGS, filename=OUTPUT_FILENAME, path=JOB_PATH)
