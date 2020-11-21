"""
Read the mechanism file
"""

import os
import sys
import mechanalyzer
import chemkin_io
import pandas as pd
import numpy as np
import copy

CWD = os.getcwd()

# possibly turn this into run.dat file?

SPC_NAME = sys.argv[1]
MECH_NAME = sys.argv[2]
SORT_NAME = sys.argv[3]


############ input reading ####################

# Read species file
with open(os.path.join(CWD, SPC_NAME), 'r') as file_obj:
    SPC_STR = file_obj.read() 

# Read input mechanism file
with open(os.path.join(CWD, MECH_NAME), 'r') as file_obj:
    MECH_STR = file_obj.read() 

# Read sorting options
SORT_STR = list(np.genfromtxt(SORT_NAME,dtype=str,comments='#'))

# (0) Build spc dct
spc_dct = mechanalyzer.parser.spc.build_spc_dct(SPC_STR,'csv')

# (1) Build pes dct, rxn block
[formulas_dct,formulas, rct_names, prd_names, rxn_names] = mechanalyzer.parser.pes.read_mechanism_file(MECH_STR,'chemkin',spc_dct)
#print(rct_names)
# extract rxn block and build species dct
block_str = chemkin_io.parser.mechanism.reaction_block(MECH_STR)
rxn_param_dct = chemkin_io.parser.reaction.param_dct(block_str)
#print(block_str)
#print(rxn_param_dct)
# the keys are the tuples with the reactants and product names
# for consistency: replace the keys with rct and prd names re-ordered
rxn_param_dct = dict(zip(list(zip(rct_names, prd_names)),rxn_param_dct.values()))

# (2) Modify the pes dct
# call a class in pes.py: store information about the PES
srt_mch = mechanalyzer.parser.pes.SORT_MECH(formulas_dct,formulas,rct_names,prd_names,rxn_names,spc_dct)
# sort according to the selected criteria
srt_mch.sort(SORT_STR)
new_idx,cmts_dct = srt_mch.return_mech_df()

# (3) Write the pes_dct to a mech str
# associate the ordered list of tuples the the rxn dictionary

new_val = list(map(rxn_param_dct.get,new_idx))
rxn_param_dct = dict(zip(new_idx,new_val))

# (4) Write the file
el_block = chemkin_io.parser.mechanism.element_block(MECH_STR)
elem_tuple = chemkin_io.parser.species.names(el_block)
chemkin_io.writer.mechanism.write_mech_file(elem_tuple,spc_dct,rxn_param_dct,filename='sorted_mech.txt',comments=cmts_dct)