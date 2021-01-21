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
# LIST OF SPECIES TO BE INCLUDED IN THE MECH; IF EMPTY: PROCESS THE FULL MECH
ISOLATE_SPECIES = ['IC8', 'SUBMECH']
# IF 'SUBMECH' ACTIVE IN ISOLATE_SPECIES: ISOLATE THE CORRESPONDING
# PYROLYSIS/OXIDATION SUBMECH (NB AVAILABLE ONLY FOR 1 SPECIES)
# ISOLATE_SPECIES = [] # leave empty if you want to process the full mechanism
# LIST WITH SORTING CRITERIA IN HIERARCHICAL ORDER
SORT_STR = ['SUBMECH', 'molecularity', 'RXN_CLASS_BROAD', 'SPECIES', 1]
# THE LAST ELEMENT IS THE N OF CRITERIA TO BE USED FOR THE HEADERS - ALSO 0 IS AVAILABLE
# AVAILABLE:
# SPECIES: GROUPS THE REACTIONS ACCORDING TO THE SPECIES LISTED IN ISOLATE_SPECIES;
# if ISOLATE_SPECIES = [], does nothing
# RXN_CLASS_BROAD: REACTION TYPE BY BROADER CRITERIA
# RXN_CLASS_GRAPH: FIND REACTION CLASS WITH GRAPH CLASSIFICATION IN AUTOMOL
# PES: ORDER BY STOICHIOMETRY
# SUBPES: ORDER BY CONNECTED CHANNELS IN 1 PES
# R1: ORDER BY FIRST REACTANT (IN BIMOL REACTIONS: ALWAYS ORDERED ACCORDING TO HIGHER N OF ATOMS)
# MULT: ORDERED BY PRODUCT OF MULTIPLICITIES OF THE REACTANTS
# MOLECULARITY: N OF REACTANTS
# SUBMECH: SUBMECH OF THE FUEL INDICATED (FUEL, RAD, RO2, RO4..)

############ input reading ####################

# READ FILE AND BUILD DICTIONARIES
spc_dct_full, rxn_param_dct, elem_tuple = mechanalyzer.parser.mech.readfiles(
    os.path.join(CWD, SPC_NAME), os.path.join(CWD, MECH_NAME))

# BUILD  MECH INFORMATION
mech_info = mechanalyzer.parser.mech.build_dct(spc_dct_full, rxn_param_dct)

# SORTING: sort the mech and build the sorted rxn param dct
sorted_idx, cmts_dct, spc_dct = mechanalyzer.parser.mech.sort_mechanism(
    mech_info, spc_dct_full, SORT_STR, ISOLATE_SPECIES)
rxn_param_dct_sorted = mechanalyzer.parser.mech.reordered_mech(
    rxn_param_dct, sorted_idx)

# WRITE THE NEW MECHANISM
mechanalyzer.parser.mech.write_mech(
    elem_tuple, spc_dct, rxn_param_dct_sorted, SORTMECH_NAME, comments=cmts_dct)

if ISOLATE_SPECIES:
    rxn_param_dct_rest = mechanalyzer.parser.util.filter_keys(
        rxn_param_dct, rxn_param_dct_sorted)
    mechanalyzer.parser.mech.write_mech(
        elem_tuple, spc_dct_full, rxn_param_dct_rest, 'mechanism_rest.txt')
