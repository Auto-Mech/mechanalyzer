import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import numpy as np
import mechanalyzer.parser.mech as mech_parser
import sys

# INPUTS
# Filenames
mech_filenames = ['chomech_v13dCLEANEDUP.txt', 'chomech_v11a_ascii_cleaned.inp']
thermo_filenames = ['therm_v11c_ascii_cleaned.dat', 'therm_v11c_ascii_cleaned.dat']
spc_csv_filenames = ['species_stereo2.csv', 'species_stereo2.csv']  

# Conditions
temps = np.linspace(300, 3000, 28)
pressures = np.array([1, 10, 100])

# Sorting; either 'rates', 'ratios', or None
sort_method = 'ratios'

# Load dcts
JOB_PATH = sys.argv[1]
rxn_ktp_dcts = compare.load_rxn_ktp_dcts_chemkin(mech_filenames, JOB_PATH, temps, pressures)
spc_thermo_dcts = compare.load_spc_thermo_dcts_chemkin(thermo_filenames, JOB_PATH, temps)
spc_ident_dcts = compare.load_spc_ident_dcts(spc_csv_filenames, JOB_PATH)

# Get the aligned_rxn_ktp_dct 
aligned_rxn_ktp_dct = compare.get_aligned_rxn_ktp_dct(
    rxn_ktp_dcts, spc_thermo_dcts, spc_ident_dcts, temps, rev_rates=True,
    remove_loners=False, write_file=True
)

# Sort as indicated in the inputs
if sort_method == 'rates':
    SORT_STR = ['molecularity', 'rxn_max_vals', 'rxn_max_ratio', 'rxn_class_broad', 0]
    ISOLATE_SPECIES = []
    mech_info = mech_parser.build_dct(spc_ident_dcts[0], aligned_rxn_ktp_dct)
    sorted_idx, _, _ = mech_parser.sort_mechanism(mech_info, spc_dct_full, SORT_STR, ISOLATE_SPECIES)
    aligned_rxn_ktp_dct = mech_parser.reordered_mech(aligned_rxn_ktp_dct, sorted_idx)
    ratio_sort = False
elif sort_method == 'ratios':
    ratio_sort = True
else: 
    ratio_sort = False

# Run the plotter    
plot_rates.build_plots(aligned_rxn_ktp_dct, mech_names=['v13', 'v11'], ratio_sort=ratio_sort)

