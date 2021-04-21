"""
Read the mechanism file
"""

import os
import automol
import chemkin_io
import mechanalyzer


CWD = os.getcwd()

# MODIFY THIS SECTION WITH INPUT NAMES AND SORTING OPTIONS

SPC_NAME = 'LLNL_iC8H18_species.csv'

MECH_NAME = 'mechanism.dat'
SORTMECH_NAME = 'sorted_mech_IC8_cl_species.txt'
# LIST OF SPECIES TO BE INCLUDED IN THE MECH; IF EMPTY: PROCESS THE FULL MECH
ISOLATE_SPECIES = ['IC8', 'submech']
# IF 'SUBMECH' ACTIVE IN ISOLATE_SPECIES: ISOLATE THE CORRESPONDING
# PYROLYSIS/OXIDATION SUBMECH (NB AVAILABLE ONLY FOR 1 SPECIES)
# ISOLATE_SPECIES = [] # leave empty if you want to process the full mechanism
# LIST WITH SORTING CRITERIA IN HIERARCHICAL ORDER
SORT_STR = ['submech', 'molecularity', 'rxn_class_broad', 'species', 1]
# LAST ELEMENT IS NUM OF CRITERIA FOR THE HEADERS - ALSO 0 IS AVAILABLE
# AVAILABLE:
# SPECIES: GROUPS REACTIONS ACCORDING TO SPECIES LISTED IN ISOLATE_SPECIES;
# if ISOLATE_SPECIES = [], does nothing
# RXN_CLASS_BROAD: REACTION TYPE BY BROADER CRITERIA
# RXN_CLASS_GRAPH: FIND REACTION CLASS WITH GRAPH CLASSIFICATION IN AUTOMOL
# PES: ORDER BY STOICHIOMETRY
# SUBPES: ORDER BY CONNECTED CHANNELS IN 1 PES
# R1: ORDER BY FIRST REACTANT
#    (IN BIMOL REACTIONS: ALWAYS ORDERED ACCORDING TO HIGHER N OF ATOMS)
# MULT: ORDERED BY PRODUCT OF MULTIPLICITIES OF THE REACTANTS
# MOLECULARITY: N OF REACTANTS
# SUBMECH: SUBMECH OF THE FUEL INDICATED (FUEL, RAD, RO2, RO4..)

# INPUT READING #

# READ FILE AND BUILD DICTIONARIES

# BUILD  MECH INFORMATION
#     os.path.join(CWD, SPC_NAME), os.path.join(CWD, MECH_NAME))
spc_dct_full = mechanalyzer.parser.spc.build_spc_dct(spc_str, spc_type)
rxn_param_dct, mech_info, elem_tuple = mechanalyzer.parser.mech.parse_mechanism(
    mech_str, mech_type, spc_dct_full)

# SORTING: sort the mech and build the sorted rxn param dct
sorted_idx, cmts_dct, spc_dct = mechanalyzer.parser.mech.sort_mechanism(
    mech_info, spc_dct_full, SORT_STR, ISOLATE_SPECIES)
rxn_param_dct_sorted = mechanalyzer.parser.mech.reordered_mech(
    rxn_param_dct, sorted_idx)

# WRITE THE NEW MECHANISM
spc_dct = mechanalyzer.parser.spc.order_species_by_atomcount(spc_dct)
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elem_tuple, spc_dct=spc_dct, rxn_param_dct=rxn_param_dct_sorted,
    comments=None)

if ISOLATE_SPECIES:
    # write a second mechanism file of the isolated species
    rxn_param_dct_rest = automol.util._dict_.filter_keys(
        rxn_param_dct, rxn_param_dct_sorted)
    mech_str2 = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=elem_tuple, spc_dct=spc_dct_full,
        rxn_param_dct=rxn_param_dct_rest,
        comments=None)

# write the files
# name = 'mechanism_rest.txt'
