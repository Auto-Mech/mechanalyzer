"""
Read the mechanism file
"""

import os
import mechanalyzer


CWD = os.getcwd()

# MODIFY THIS SECTION WITH INPUT NAMES AND SORTING OPTIONS

SPC_NAME = 'LLNL_iC8H18_species.csv'
MECH_NAME = 'mechanism.dat'
SORTMECH_NAME = 'sorted_mech_IC8_cl_species.txt'
ISOLATE_SPECIES = ['IC8','IC8-1R','IC8-3R','IC8-4R','IC8-5R'] # LIST OF SPECIES TO BE INCLUDED IN THE MECH; IF EMPTY: PROCESS THE FULL MECH
# ISOLATE_SPECIES = [] # leave empty if you want to process the full mechanism
SORT_STR = ['RXN_CLASS_GRAPH','SPECIES',1] # LIST WITH SORTING CRITERIA IN HIERARCHICAL ORDER
                # THE LAST ELEMENT IS THE N OF CRITERIA TO BE USED FOR THE HEADERS - ALSO 0 IS AVAILABLE
# AVAILABLE:
# SPECIES: GROUPS THE REACTIONS ACCORDING TO THE SPECIES LISTED IN ISOLATE_SPECIES; if ISOLATE_SPECIES = [], does nothing
# RXN_CLASS_GRAPH: FIND REACTION CLASS WITH GRAPH CLASSIFICATION IN AUTOMOL
# PES: ORDER BY STOICHIOMETRY
# SUBPES: ORDER BY CONNECTED CHANNELS IN 1 PES
# R1: ORDER BY FIRST REACTANT (IN BIMOL REACTIONS: ALWAYS ORDERED ACCORDING TO HIGHER N OF ATOMS)
# MULT_R1: ORDERED BY MULTIPLICITY OF THE FIRST REACTANT
# numC: ORDERED BY NUMBER OF CARBON ATOMS IN THE MOLECULE/TOT IN THE BIMOL REACTANTS
# specieslist: ALL REACTIONS INVOLVING THE SPECIES IN THE LIST

############ input reading ####################

# READ FILE
SPC_STR,MECH_STR = mechanalyzer.parser.mech.readfiles(os.path.join(CWD,SPC_NAME),os.path.join(CWD,MECH_NAME))

# BUILD DICTIONARIES AND MECH INFORMATION
spc_dct,mech_info,rxn_param_dct = mechanalyzer.parser.mech.build_dct(SPC_STR,MECH_STR)

# SORTING: sort the mech and build the sorted rxn param dct
sorted_idx,cmts_dct = mechanalyzer.parser.mech.sort_mechanism(mech_info,spc_dct,SORT_STR,ISOLATE_SPECIES)
rxn_param_dct_sorted = mechanalyzer.parser.mech.reordered_mech(rxn_param_dct,sorted_idx,cmts_dct)

# WRITE THE NEW MECHANISM
mechanalyzer.parser.mech.write_reordered_mech(MECH_STR,spc_dct,rxn_param_dct_sorted,cmts_dct,SORTMECH_NAME)

