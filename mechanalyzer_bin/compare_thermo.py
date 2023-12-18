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

OUTPUT_FILENAME = 'nuig_m2-m5.pdf'
#OUTPUT_FILENAME = 'gold.pdf'
#OUTPUT_FILENAME = 'hr_hrf_rmg.pdf'

THERMO_FILENAMES = [
    #'NUIGMech1.2.THERM',
    #'frank.therm',
    'NUIGMech1.2.THERM',
    #'to_us.ckin',
    'm2.ckin',
    'm3.ckin',
    'm4.ckin',
    'm5_all.ckin',
]
SPC_CSV_FILENAMES = [
    #'canon_NUIG_species.csv',
    #'frank_species.csv',
    'canon_NUIG_species.csv',
    #'yao_spc.csv',
    'canon_2.csv',
    'canon_2.csv',
    'canon_2.csv',
    'canon_2.csv',
    #'canon_3.csv',
    #'canon_3.csv',
    #'canon_3.csv',
    #'canon_species.csv',
    #'canon_species.csv',
]
MECH_NAMES = [
    #'NUIG',
    #'Goldsmith',
    'NUIGMech1.2',
    #'Yao',
    'M2',
    'M3',
    'M4',
    'M5',
]

# Conditions
TEMPS = numpy.linspace(300, 1500, 13)

# Options
# SORT = False
SORT = True
SORT_INSTR = 'lnq'  # either 'h', 'cp', 's', 'g', 'lnq', or None
SORT_TEMP = None  # can be (1) None to sort by max difference or (2) a number
REMOVE_LONERS = True
WRITE_FILE = False  # this currently does nothing
PRINT_MISSING = True  # print spcs that are in mech(s) but not in spc.csv
OUT_FILENAME = 'comparetool.out'  #

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
pathtools.write_file(FSTR, JOB_PATH, OUT_FILENAME)
# Saving this for later...writes a comparison of species
#comp_str = compare.write_comparison(ALGN_SPC_THERM_DCT, dct_type='therm')
