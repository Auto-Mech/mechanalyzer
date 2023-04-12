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

# Filenames
SPC_CSV = 'species.csv'
CKIN = 'amech.ckin'
THERM = 'amech.therm'

# Load things
MECH_SPC_DCT = spc_parser.load_mech_spc_dct(SPC_CSV, DAT_PATH, canon_ent=True)
RXN_PARAM_DCT = ckin_parser.load_rxn_param_dct(CKIN, DAT_PATH)
SPC_NASA7_DCT = ckin_parser.load_spc_nasa7_dct(THERM, DAT_PATH)

spc = 'C4CEHY-2m1jatAB0'
spc_dct = MECH_SPC_DCT[spc]
inchi = spc_dct['canon_enant_ich']
racem = amchi.racemic(inchi)

spcs = ['C4CEHY-Y2VjatAB0', 'C4CEHY-Y2VjatAB1', 'C4CEHY-Y2VjatAA0',
        'C4CEHY-Y2VjatAA1']
for spc in spcs:
    spc_dct = MECH_SPC_DCT[spc]
    inchi = spc_dct['canon_enant_ich']
    racem = amchi.racemic(inchi)
    print('inchi: ', inchi)
    print('racem: ', racem)

for spc in MECH_SPC_DCT:
    if spc[-1] in ('E', 'Z'):
        print('E or Z spc: ', spc)

CANON_ENT = False
#mech_spc_dct_strpd, mech_spc_dct_no_ste = strip_ste.strip_mech_spc_dct(
#    MECH_SPC_DCT, canon_ent=CANON_ENT)

iso_sets = racemize.find_iso_sets(MECH_SPC_DCT, canon_ent=CANON_ENT)
print('iso_sets')
[print(iso_set) for iso_set in iso_sets]

rac_sets = racemize.get_rac_sets(iso_sets, MECH_SPC_DCT)

rxns_by_idx = racemize.get_rxns_by_idx(rac_sets, RXN_PARAM_DCT)

for rxn_by_idx, params_lst in rxns_by_idx.items():
    if len(params_lst) > 1:
        print('rxn_by_idx: ', rxn_by_idx)

print(rxns_by_idx[((41,),(79,),(None,))])
print(rxns_by_idx[((79,),(41,),(None,))])

print('more than one ste species')
for rxn in rxns_by_idx:
    rcts, prds, _ = rxn
    num_ste = 0  # number of chiral 
    for rct in rcts:
        if len(rac_sets[rct]) > 1:
            num_ste += 1
    for prd in prds:
        if len(rac_sets[prd]) > 1:
            num_ste += 1
    if num_ste > 1:
        print(rxn)

racemize.lump(rxns_by_idx)

breakpoint()

# Write the mechanism to a Chemkin file
#mech_str = mechanism.write_chemkin_file(
#    rxn_param_dct=new_rxn_param_dct,
#    mech_spc_dct=new_mech_spc_dct,
#    spc_nasa7_dct=new_spc_nasa7_dct)
#pathtools.write_file(mech_str, DAT_PATH, 'merge_ste.out')

