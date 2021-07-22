""" Script for running a comparison of thermodynamic properties between mechanisms
"""

import sys
import numpy
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.thermo as plot_thermo
import mechanalyzer.plotter._util as util
import mechanalyzer.parser.spc as spc_parser  
import mechanalyzer.parser.ckin_ as ckin_parser

# INPUTS
# Filenames
thermo_filenames = ['NUIGMech1.2.new.THERM', 'NUIGMech1.2.old.THERM']
spc_csv_filenames = ['NUIG_species.csv', 'NUIG_species.csv']

output_filename = 'compare_thermo.pdf'
mech_nicknames = ['ANL', 'NUIG']

# Conditions
temps = numpy.linspace(300,1000, 15)

# options
sort = True
sort_instr = 's' # either 'h', 'cp', 's', 'g', or None
remove_loners = True
write_file = False


# RUN FUNCTIONS
# Load dcts
assert len(sys.argv) == 2, (
    'There should be one command line input; namely, the job path')
JOB_PATH = sys.argv[1]
spc_therm_dcts = ckin_parser.load_spc_therm_dcts(thermo_filenames, JOB_PATH, 
                                                 temps)
spc_dcts = spc_parser.load_spc_dcts(spc_csv_filenames, JOB_PATH)

# Get the algn_spc_therm_dct
algn_spc_therm_dct = compare.get_algn_spc_therm_dct(
    spc_therm_dcts, spc_dcts, remove_loners=remove_loners,
    write_file=write_file)

# Get the combined spc_dct (used for including SMILES and InChis)
comb_spc_dct = compare.get_mult_comb_spc_dct(spc_dcts)

# Run the plotter
figs = plot_thermo.build_plots(
    algn_spc_therm_dct, spc_dct=comb_spc_dct, mech_names=mech_nicknames,
    sort=sort, sort_instr=sort_instr)
util.build_pdf(figs, filename=output_filename, path=JOB_PATH)

