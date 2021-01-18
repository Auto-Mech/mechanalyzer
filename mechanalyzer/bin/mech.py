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
# LIST OF SPECIES TO BE INCLUDED IN THE MECH; IF EMPTY: PROCESS THE FULL MECH
ISOLATE_SPECIES = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
# ISOLATE_SPECIES = [] # leave empty if you want to process the full mechanism
# LIST WITH SORTING CRITERIA IN HIERARCHICAL ORDER
SORT_STR = ['RXN_CLASS_GRAPH', 'SPECIES', 1]
# THE LAST ELEMENT IS THE N OF CRITERIA TO BE USED FOR THE HEADERS - ALSO 0 IS AVAILABLE
# AVAILABLE:
# SPECIES: GROUPS THE REACTIONS ACCORDING TO THE SPECIES LISTED IN ISOLATE_SPECIES;
# if ISOLATE_SPECIES = [], does nothing
# RXN_CLASS_GRAPH: FIND REACTION CLASS WITH GRAPH CLASSIFICATION IN AUTOMOL
# PES: ORDER BY STOICHIOMETRY
# SUBPES: ORDER BY CONNECTED CHANNELS IN 1 PES
# R1: ORDER BY FIRST REACTANT (IN BIMOL REACTIONS: ALWAYS ORDERED ACCORDING TO HIGHER N OF ATOMS)
# MULT: ORDERED BY MULTIPLICITY OF  REACTANTS mR1*mR2..
# numC: ORDERED BY NUMBER OF CARBON ATOMS IN THE MOLECULE/TOT IN THE BIMOL REACTANTS
# specieslist: ALL REACTIONS INVOLVING THE SPECIES IN THE LIST

############ input reading ####################

# READ FILE AND BUILD DICTIONARIES
spc_dct, rxn_param_dct, elem_tuple = mechanalyzer.parser.mech.readfiles(
    os.path.join(CWD, SPC_NAME), os.path.join(CWD, MECH_NAME))

# BUILD  MECH INFORMATION
mech_info = mechanalyzer.parser.mech.build_dct(spc_dct, rxn_param_dct)

# SORTING: sort the mech and build the sorted rxn param dct
sorted_idx, cmts_dct, spc_dct = mechanalyzer.parser.mech.sort_mechanism(
    mech_info, spc_dct, SORT_STR, ISOLATE_SPECIES)
rxn_param_dct_sorted = mechanalyzer.parser.mech.reordered_mech(
    rxn_param_dct, sorted_idx, cmts_dct)

# WRITE THE NEW MECHANISM
mechanalyzer.parser.mech.write_reordered_mech(
    elem_tuple, spc_dct, rxn_param_dct_sorted, cmts_dct, SORTMECH_NAME)
