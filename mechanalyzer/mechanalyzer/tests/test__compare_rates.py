import mechanalyzer.calculator.compare as compare
import mechanalyzer.plotter.rates as plot_rates
import numpy as np
import sys

# Input filenames
mech_filenames = ['data/mech1.txt', 'data/mech2.txt']
thermo_filenames = ['data/thermo1.txt', 'data/thermo2.txt']
spc_csv_filenames = ['data/spc1.csv', 'data/spc2.csv']

# Input calculation conditions
temps = np.array([500, 1000, 1500])
#temps = np.linspace(500, 2000, 31)
pressures = np.array([1, 10, 100])

# Load dcts
JOB_PATH = sys.argv[1]
rxn_ktp_dcts = compare.load_rxn_ktp_dcts_chemkin(mech_filenames, JOB_PATH, temps, pressures)
spc_thermo_dcts = compare.load_spc_thermo_dcts_chemkin(thermo_filenames, JOB_PATH, temps)
spc_ident_dcts = compare.load_spc_ident_dcts(spc_csv_filenames, JOB_PATH)

# Get the aligned_rxn_ktp_dct 
aligned_rxn_ktp_dct = compare.get_aligned_rxn_ktp_dct(
    rxn_ktp_dcts, spc_thermo_dcts, spc_ident_dcts, temps, rev_rates=True, 
    remove_loners=True, write_file=True
)

