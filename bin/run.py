""" Run some analysis function
"""

import os
import numpy as np
import plot
import branching
import stability
from _format import read_file


# Set files
PATH = os.path.dirname(os.path.realpath(__file__))
CKIN_FILES = [
    # './abs/l4_fin/full.ckin',
    # './abs/l3_ckin/full.ckin',
    # './abs/l2_ckin/full.ckin',
    # './abs/l1_ckin/full.ckin',
    # './mig/l4_fin/full.ckin',
    # './mig/l3_refit/full.ckin',
    # './mig/l2_ckin/full.ckin',
    # './mig/l1_ckin/full.ckin',
    # './think/test.dat',
    './think/think.mech',
    './think/aramco2.mech'
    # './think/think_small.mech',
    # './think/aramco2_small.mech'
]
CSV_FILE1 = './think/think.csv'
CSV_FILE2 = './think/aramco2.csv'
# CSV_FILE = './abs/species.csv'
# CSV_FILE = 'mig/species.csv'

# Set temperatures and pressures
T_REF = 1.0
# TEMPS = np.array([1500.0])
TEMPS = np.arange(500.0, 2025.0, 25.00)
# TEMPS = np.arange(300.0, 3025.0, 25.00)
# TEMPS = np.arange(500.0, 1600.0, 100.00)
# TEMPS = np.array([300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0])
# PRESSURES = [0.1, 1.0, 10.0, 100.0]
PRESSURES = [0.1, 1.0, 10.0, 100.0, 'high']
# PRESSURES = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 'high']
# PRESSURES = [1.0]
# PRESSURES = ['high']
ALLOW_RCTS = []
# ALLOW_RCTS = ['H', 'OH', 'CH3']
MLBLS = ['TK', 'AR']
RATE_PREFIX = '.'
THERMO_PREFIX = '.'

# Combine
# comb.combine_mech_files(CKIN_FILES)

# Rates and Branching
# branching.calc_multimech_rates2(
#     TEMPS, PRESSURES, T_REF, CKIN_FILES, PATH,
#     rtyp='ignore', allow_rcts=ALLOW_RCTS)

# Stability
# stability.read_radical_stabilities2(PRESSURES, CKIN_FILES, PATH)

# single Plotter
# MECH1_STR = read_file(os.path.join(PATH, CKIN_FILES[0]))
# plot.sm_rates(MECH1_STR, T_REF, TEMPS, PRESSURES, dir_prefix=RATE_PREFIX,
#               ignore_reverse=True, remove_bad_fits=True)

# multi Plotter
MECH1_STR = read_file(os.path.join(PATH, CKIN_FILES[0]))
MECH2_STR = read_file(os.path.join(PATH, CKIN_FILES[1]))
# MECH2_STR = ''
with open(CSV_FILE1, 'r') as csvfile:
    CSV_STR1 = csvfile.read()
with open(CSV_FILE2, 'r') as csvfile:
    CSV_STR2 = csvfile.read()
plot.rates(MECH1_STR, MECH2_STR, T_REF, TEMPS, PRESSURES,
           mech1_csv_str=CSV_STR1, mech2_csv_str=CSV_STR2,
           dct_idx='inchi', dir_prefix=RATE_PREFIX,
           mech_labels=MLBLS)
