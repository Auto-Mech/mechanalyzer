import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import numpy as np
import mechanalyzer.parser.mech as mech_parser
import sys

# INPUTS
# Filenames
mech_filenames = ['Tran_v33_short.ckin', 'Tran_v42_short.ckin']
thermo_filenames = ['Tran_v33_short.ckin', 'Tran_v42_short.ckin']
spc_csv_filenames = ['Tran_species_dumb.csv', 'Tran_species_dumb.csv']

output_filename = '33_42_short.pdf'
mech_nicknames = ['v33', 'v42']

# Conditions
temps = np.linspace(400, 700, 31)
pressures = np.array([1, 10])

# options
sort_method = 'ratios' # sorting; either 'rates', 'ratios', or None
rev_rates = True
remove_loners = True
write_file = False


# RUN FUNCTIONS
# Load dcts
assert len(sys.argv) == 2, (
    'There should be one input specified on the command line, namely the job path'
)
JOB_PATH = sys.argv[1]
rxn_ktp_dcts = compare.load_rxn_ktp_dcts_chemkin(mech_filenames, JOB_PATH, temps, pressures)
spc_thermo_dcts = compare.load_spc_thermo_dcts_chemkin(thermo_filenames, JOB_PATH, temps)
spc_ident_dcts = compare.load_spc_ident_dcts(spc_csv_filenames, JOB_PATH)

# Get the aligned_rxn_ktp_dct 
aligned_rxn_ktp_dct = compare.get_aligned_rxn_ktp_dct(
    rxn_ktp_dcts, spc_thermo_dcts, spc_ident_dcts, temps, rev_rates=rev_rates,
    remove_loners=remove_loners, write_file=write_file
)

# Sort as indicated in the inputs
if sort_method == 'rates':
    SORT_STR = ['molecularity', 'rxn_max_vals', 'rxn_max_ratio', 'rxn_class_broad', 0]
    ISOLATE_SPCS = []
    mech_info = mech_parser.build_dct(spc_ident_dcts[0], aligned_rxn_ktp_dct)
    sorted_idx, _, _ = mech_parser.sort_mechanism(mech_info, spc_dct_full, SORT_STR, ISOLATE_SPCS)
    aligned_rxn_ktp_dct = mech_parser.reordered_mech(aligned_rxn_ktp_dct, sorted_idx)
    ratio_sort = False
elif sort_method == 'ratios':
    ratio_sort = True
else: 
    ratio_sort = False

# Run the plotter    
plot_rates.build_plots(aligned_rxn_ktp_dct, path=JOB_PATH, 
                       filename=output_filename,
                       mech_names=mech_nicknames, ratio_sort=ratio_sort)

