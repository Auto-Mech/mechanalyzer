""" Script for running a comparison of thermodynamic properties between mechanisms
"""

import sys
import numpy 
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.thermo as plot_thermo
import mechanalyzer.plotter._util as util
import mechanalyzer.parser.mech as mech_parser

# INPUTS
# Filenames
thermo_filenames = ['glarborg.therm', 'stagni.therm']
spc_csv_filenames = ['nh3_species.csv', 'nh3_species.csv',]

output_filename = 'test_thermo.pdf'
mech_nicknames = ['glarborg', 'stagni']

# Conditions
temps = numpy.linspace(400, 2000, 17)

# options
sort_method = 'h' # either 'h', 'cp', 's', 'g', or None
remove_loners = True
write_file = False


# RUN FUNCTIONS
# Load dcts
assert len(sys.argv) == 2, (
    'There should be one command line input; namely, the job path')

JOB_PATH = sys.argv[1]
spc_therm_dcts = compare.load_spc_therm_dcts_chemkin(thermo_filenames, JOB_PATH, temps)
spc_dcts = compare.load_spc_dcts(spc_csv_filenames, JOB_PATH)

# Get the algn_spc_therm_dct 
algn_spc_therm_dct = compare.get_algn_spc_therm_dct(
    spc_therm_dcts, spc_dcts, remove_loners=remove_loners,
    write_file=write_file)

# Run the plotter    
figs = plot_thermo.build_plots(
    algn_spc_therm_dct, path=JOB_PATH, filename=output_filename,
    mech_names=mech_nicknames, ratio_sort=ratio_sort)
util.build_pdf(figs, filename=output_filename, path=JOB_PATH)

