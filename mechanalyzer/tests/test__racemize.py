""" This script tests the merging of a stereo-specific submechanism (e.g., one
    calculated with AutoMech) with a non-stereo mechanism from the literature 
"""

import os
import tempfile
import numpy
from mechanalyzer.parser import new_spc as spc_parser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.builder import merge_ste
from mechanalyzer.builder import strip_ste
from mechanalyzer.builder import racemize
from mechanalyzer.calculator.rates import check_p_t
from chemkin_io.writer import mechanism
from ioformat import pathtools
from autoreact import params
from automol import amchi

# Set Paths to test/data directory and output directory
DAT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
TMP_OUT = tempfile.mkdtemp()

PRESSURES = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100]
TEMPS = [numpy.arange(360, 1500, 60)]

# Filenames
SPC_CSV = 'species.csv'
#CKIN = 'amech.ckin'
#CKIN = 'amech_short.ckin'
CKIN = 'amech_140423.ckin'
THERM = 'amech.therm'

# Load things
MECH_SPC_DCT = spc_parser.load_mech_spc_dct(SPC_CSV, DAT_PATH, canon_ent=True)
RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(CKIN, DAT_PATH)
SPC_NASA7_DCT = ckin_parser.load_spc_nasa7_dct(THERM, DAT_PATH)

# Get isomer sets, racemic_sets, and racemic rxn_param_dct
iso_sets = racemize.find_iso_sets(MECH_SPC_DCT)
rac_sets, rac_names, rac_mech_spc_dct = racemize.get_rac_sets(
    iso_sets, MECH_SPC_DCT)
rac_rxn_param_dct = racemize.get_rac_rxn_param_dct(
    rac_sets, rac_names, RXN_PARAM_DCT)

# Run the check for unbalanced reactions
racemize.check_bal_rxns(rac_rxn_param_dct, rac_mech_spc_dct)
breakpoint()

# Get the lumped rxn parameter dictionary
lump_rxn_param_dct = racemize.lump(rac_rxn_param_dct, TEMPS, PRESSURES)

# Write to file
mech_str = mechanism.write_chemkin_file(
    mech_spc_dct=rac_mech_spc_dct, rxn_param_dct=lump_rxn_param_dct)
pathtools.write_file(mech_str, DAT_PATH, 'racemize.out')
#csv_str = 

