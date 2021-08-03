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
# thermo_filenames = ['fa_07.THERM', 'NUIGMech1.2.THERM']
# spc_csv_filenames = ['C5-C8_species.csv', 'NUIG_stereo.csv']
# thermo_filenames = ['nc5_fa_07.THERM', 'NUIG_nc5.THERM']
# spc_csv_filenames = ['nc5_species.csv', 'NUIG_nc5_species.csv']
# thermo_filenames = ['NUIG_nc5.THERM', 'nc5_fa_07.THERM', 'nc5_f.THERM']
# spc_csv_filenames = ['NUIG_nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv']
# thermo_filenames = ['NUIG_nc5.THERM', 'nc5_fa_07.THERM', 'nc5_f.THERM', 'nc5_f_wbt.THERM']
# spc_csv_filenames = ['NUIG_nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv']
# thermo_filenames = ['NUIG_nc5.THERM', 'nc5_1-5_fa_07.THERM', 'nc5_1-5_f.THERM', 'nc5_1-5_1dhr.THERM']
# spc_csv_filenames = ['NUIG_nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv']
# thermo_filenames = ['NUIG_nc5.THERM', '1-5_wbt_fa_07.THERM', '1-5_wbt_f.THERM', '1-5_wbt_1dhr.THERM']
# spc_csv_filenames = ['NUIG_nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv']
# THERMO_FILENAMES = ['NUIG_nc5.THERM', 'nc5_1-5_f.THERM', '1-5_wbt_f.THERM', '1-5_b3sp_f.THERM']
# THERMO_FILENAMES = ['NUIG_nc5.THERM', 'nc5_1-5_1dhr.THERM', '1-5_wbt_1dhr.THERM', '1-5_b3sp_1dhr.THERM']
# SPC_CSV_FILENAMES = ['NUIG_nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv', 'nc5_species.csv']
# THERMO_FILENAMES = ['NUIG_nc3.THERM', '1-12_wbs_fa.THERM', '1-12_wbs_f.THERM', '1-12_wbs_hr.THERM']
# SPC_CSV_FILENAMES = ['NUIG_nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv']
# THERMO_FILENAMES = ['NUIG_nc3.THERM', 'conf1.THERM', 'conf2.THERM', 'conf3.THERM', 'conf4.THERM', 'conf5.THERM']
# THERMO_FILENAMES = ['NUIG_nc3.THERM', 'all_therm.ckin_0']
# SPC_CSV_FILENAMES = ['NUIG_nc3_species.csv', 'nc3_species.csv']
THERMO_FILENAMES = ['NUIG_nc3.THERM', 'all_therm.ckin_0', 'all_therm.ckin_1', 'all_therm.ckin_2', 'all_therm.ckin_3', 'all_therm.ckin_4']
SPC_CSV_FILENAMES = ['NUIG_nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv']
# THERMO_FILENAMES = ['NUIG_nc3.THERM', '1-12_b3sp_fa.THERM', '1-12_b3sp_f.THERM', '1-12_b3sp_hr.THERM']
# THERMO_FILENAMES = ['NUIG_nc3.THERM', '1-12_wbs_hr.THERM', '1-12_b3sp_hr.THERM']
# SPC_CSV_FILENAMES = ['NUIG_nc3_species.csv', 'nc3_species.csv', 'nc3_species.csv']

OUTPUT_FILENAME = 'compare_thermo.pdf'
# MECH_NICKNAMES = ['NUIG', 'ANL_wbs_hr', 'ANL_b3_hr']
# MECH_NICKNAMES = ['NUIG', 'ANL_wbs_fa', 'ANL_wbs_f', 'ANHL_wbs_hr']
# MECH_NICKNAMES = ['NUIG', 'ANL_b3sp_fa', 'ANL_b3sp_f', 'ANHL_b3sp_hr']
# MECH_NICKNAMES = ['NUIG', 'wbs_hr']
MECH_NICKNAMES = ['NUIG', 'wbs_hr', 'conf2', 'conf3', 'conf4', 'conf5']
#mech_nicknames = ['NUIG', 'ANL_fa_07', 'ANL_f', 'ANL_hr']

# Conditions
TEMPS = numpy.linspace(300, 1000, 15)

# options
SORT = False
# sort = True
SORT_INSTR = 'h' # either 'h', 'cp', 's', 'g', or None
REMOVE_LONERS = False
WRITE_FILE = False


# RUN FUNCTIONS
# Load dcts
assert len(sys.argv) == 2, (
    'There should be one command line input; namely, the job path')
if len(sys.argv) > 1: 
    JOB_PATH = sys.argv[1]
else:
    JOB_PATH = os.getcwd()
print('file test:', THERMO_FILENAMES, JOB_PATH)
SPC_THERM_DCTS = ckin_parser.load_spc_therm_dcts(THERMO_FILENAMES, JOB_PATH,
                                                 TEMPS)
SPC_DCTS = spc_parser.load_spc_dcts(SPC_CSV_FILENAMES, JOB_PATH)

# Get the algn_spc_therm_dct
ALGN_SPC_THERM_DCT = compare.get_algn_spc_therm_dct(
    SPC_THERM_DCTS, SPC_DCTS, remove_loners=REMOVE_LONERS,
    write_file=WRITE_FILE)

# Get the combined spc_dct (used for including SMILES and InChis)
COMB_SPC_DCT = compare.get_mult_comb_spc_dct(SPC_DCTS)

# Run the plotter
FIGS = plot_thermo.build_plots(
    ALGN_SPC_THERM_DCT, spc_dct=COMB_SPC_DCT, mech_names=MECH_NICKNAMES,
    sort=SORT, sort_instr=SORT_INSTR)
util.build_pdf(FIGS, filename=OUTPUT_FILENAME, path=JOB_PATH)
