import numpy as np
import sys
import ioformat.pathtools as fileio
import mechanalyzer.parser.ckin_ as ckin_parser
from mechanalyzer.builder import checker

# INPUTS
MECH_FILENAME = 'mechanism.dat'
TEMPS = [np.linspace(300, 3000, 28)]
PRESSURES = [1, 10, 100]
K_THRESHOLDS = [1e11, 1e15, 1e22]
RXN_NUM_THRESHOLD = 2
OUT_FILENAME = 'mech_check.txt'

# Load dcts
JOB_PATH = sys.argv[1]
RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(MECH_FILENAME, JOB_PATH)
RXN_KTP_DCT = ckin_parser.load_rxn_ktp_dct(MECH_FILENAME, JOB_PATH,
                                           TEMPS, PRESSURES)

output_str = checker.run_all_checks(
    RXN_PARAM_DCT, RXN_KTP_DCT, K_THRESHOLDS, RXN_NUM_THRESHOLD)
fileio.write_file(output_str, JOB_PATH, OUT_FILENAME)
