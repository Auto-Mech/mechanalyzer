""" Script for running a comparison of rate constants between mechanisms
"""

import os
import sys
import numpy
import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import mechanalyzer.plotter._util as util
import mechanalyzer.parser.spc as spc_parser
import mechanalyzer.parser.ckin_ as ckin_parser
import ratefit

# INPUTS
# Filenames
mech_filenames = ['Hong_MECH_no_troe.ckin',]
thermo_filenames = ['Hong_THERM.ckin',]
spc_csv_filenames = ['H2_O2_species.csv',]

output_filename = 'fit_comp.pdf'
mech_nicknames = ['original', 'refit']

# Conditions
temps = numpy.linspace(400, 2000, 17)
pressures = numpy.array([1, 10])

# options
sort_method = 'ratios' # either 'ratios' or None
rev_rates = True
remove_loners = True
write_file = False

# RUN FUNCTIONS
# Get the job path and load the dcts
if len(sys.argv) > 1:
    JOB_PATH = sys.argv[1]
    print(f'The job path is {JOB_PATH}')
elif len(sys.argv) == 1:
    JOB_PATH = os.getcwd()
    print(f'No job path input; using the current directory, {JOB_PATH}')

# Get the therm and spc dcts for mechanism aligning 
spc_therm_dcts = ckin_parser.load_spc_therm_dcts(
    thermo_filenames, JOB_PATH, temps)
spc_dcts = spc_parser.load_spc_dcts(spc_csv_filenames, JOB_PATH)

# Read the initial ktp dct and generate a new one with the fitter 
rxn_ktp_dcts = ckin_parser.load_rxn_ktp_dcts(
    mech_filenames, JOB_PATH, temps, pressures)
rxn_ktp_dct = rxn_ktp_dcts[0]

new_rxn_ktp_dct = {}
for rxn, ktp_dct in rxn_ktp_dct.items():
    # Dummy variables
    fake_path = os.getcwd() 
    reaction = 'W1=P1'  # to multiply by 1
    new_rxn_ktp_dct[rxn] = ratefit.fit.arrhenius.pes(
        ktp_dct, reaction, fake_path,
        dbltol=15.0,
        dblcheck='max')

new_rxn_ktp_dcts = [new_rxn_ktp_dct, rxn_ktp_dct]
new_spc_therm_dcts = [spc_therm_dcts[0], spc_therm_dcts[0]]
new_spc_dcts = [spc_dcts[0], spc_dcts[0]]

# Get the algn_rxn_ktp_dct
algn_rxn_ktp_dct = compare.get_algn_rxn_ktp_dct(
    new_rxn_ktp_dcts, new_spc_therm_dcts, new_spc_dcts, temps, rev_rates=rev_rates,
    remove_loners=remove_loners, write_file=write_file
)

print(algn_rxn_ktp_dct)
print(len(algn_rxn_ktp_dct))

if sort_method == 'ratios':
    ratio_sort = True
else:
    ratio_sort = False

# Run the plotter
figs = plot_rates.build_plots(algn_rxn_ktp_dct, mech_names=mech_nicknames,
                              ratio_sort=ratio_sort)
util.build_pdf(figs, filename=output_filename, path=JOB_PATH)
