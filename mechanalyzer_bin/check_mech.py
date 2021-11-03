import numpy as np
import sys
import ioformat.pathtools as fileio
from chemkin_io.parser import mechanism as parser_mech
from chemkin_io.parser import reaction as parser_rxn
from mechanalyzer.calculator import rates as calc_rates
from mechanalyzer.builder import checker

# INPUTS
mech_filename = 'chomech_v11a_ascii_cleaned.inp'
temps = np.linspace(300, 3000, 28)
pressures = np.array([1, 10, 100])
k_thresholds = [1e11, 1e15, 1e22]
rxn_num_threshold = 2
output_filename = 'mech_check.txt'

# Load dcts
JOB_PATH = sys.argv[1]
mech_str = fileio.read_file(JOB_PATH, mech_filename)
ea_units, a_units = parser_mech.reaction_units(mech_str)
rxn_block_str = parser_mech.reaction_block(mech_str)
RXN_PARAM_DCT = parser_rxn.param_dct(rxn_block_str, ea_units, a_units)
RXN_KTP_DCT = calc_rates.eval_rxn_param_dct(RXN_PARAM_DCT, pressures, temps)

output_str = checker.run_all_checks(
    RXN_PARAM_DCT, RXN_KTP_DCT, k_thresholds, rxn_num_threshold)
fileio.write_file(output_str, JOB_PATH, output_filename)
