""" Script for running a comparison of rate constants between mechanisms
"""

import sys
import numpy
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import mechanalyzer.plotter._util as util
import mechanalyzer.parser.spc as spc_parser
import mechanalyzer.parser.ckin_ as ckin_parser

# INPUTS
# Filenames
mech_filenames = ['Tran_v33_short.ckin', 'Tran_v42_short.ckin']
thermo_filenames = ['Tran_v33_short.ckin', 'Tran_v42_short.ckin']
spc_csv_filenames = ['Tran_species_dumb.csv', 'Tran_species_dumb.csv']

output_filename = '33_42_short.pdf'
mech_nicknames = ['v33', 'v42']

# Conditions
temps = numpy.linspace(400, 700, 31)
pressures = numpy.array([1, 10])

# options
sort_method = 'ratios' # either 'ratios' or None
rev_rates = True
remove_loners = True
write_file = False


# RUN FUNCTIONS
# Load dcts
assert len(sys.argv) == 2, (
    'There should be one command line input; namely, the job path')

JOB_PATH = sys.argv[1]
rxn_ktp_dcts = ckin_parser.load_rxn_ktp_dcts(
    mech_filenames, JOB_PATH, temps, pressures)
spc_therm_dcts = ckin_parser.load_spc_therm_dcts(
    thermo_filenames, JOB_PATH, temps)
spc_dcts = spc_parser.load_spc_dcts(spc_csv_filenames, JOB_PATH)

# Get the algn_rxn_ktp_dct
algn_rxn_ktp_dct = compare.get_algn_rxn_ktp_dct(
    rxn_ktp_dcts, spc_therm_dcts, spc_dcts, temps, rev_rates=rev_rates,
    remove_loners=remove_loners, write_file=write_file
)

if sort_method == 'ratios':
    ratio_sort = True
else:
    ratio_sort = False

# Run the plotter
figs = plot_rates.build_plots(algn_rxn_ktp_dct, mech_names=mech_nicknames,
                              ratio_sort=ratio_sort)
util.build_pdf(figs, filename=output_filename, path=JOB_PATH)


