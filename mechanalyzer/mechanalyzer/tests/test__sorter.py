"""
test the sorter
"""

import os
import mechanalyzer


CWD = os.getcwd()
CWD_RESULTS = os.path.join(CWD, 'sorter_results')
if not os.path.isdir(CWD_RESULTS):
    os.makedirs(CWD_RESULTS)

# FUNCTIONS TO PROCESS DIFFERENT MECHANISMS IN data/ WITH DIFFERENT SORTING OPTIONS


def test__readwrite_thirdbody():
    """
    check read/write of a small set of rxns involving third bodies
    """
    SPC_NAME = os.path.join(CWD, 'data', 'NUIG_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    MECH_REST_NAME = os.path.join(CWD_RESULTS, 'NUIG_mech_rest.txt')
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'NUIG_test_readwrite_thirdbody.txt')
    ISOLATE_SPECIES = []
    SORT_STR = ['pes', 0]  # ARRANGE BY PES- NO HEADERS INCLUDED

    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)


def test__sortby_mult():
    """
    sort by multiplicity of the reaction
    """
    SPC_NAME = os.path.join(CWD, 'data', 'C10H10_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'C10H10_HP_mech.dat')
    MECH_REST_NAME = os.path.join(CWD_RESULTS, 'C10H10_mech_rest.txt')
    SORTMECH_NAME = os.path.join(CWD_RESULTS, 'C10H10_test_sortby_mult.txt')
    ISOLATE_SPECIES = []
    SORT_STR = ['mult', 0]  # NO HEADERS INCLUDED

    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)


def test__sortby_molec_R1():
    """
    sort by first (heavier) reactant and molecularity of the reaction
    """
    SPC_NAME = os.path.join(CWD, 'data', 'C10H10_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'C10H10_HP_mech.dat')
    MECH_REST_NAME = os.path.join(CWD_RESULTS, 'C10H10_mech_rest.txt')
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'C10H10_test_sortby_molec_R1.txt')
    ISOLATE_SPECIES = []
    SORT_STR = ['r1', 'molecularity', 0]  # NO HEADERS INCLUDED

    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)


def test__sortby_rxnclass():
    """
    sort by reaction class: both "broad" (based on multiplicity, reactants/products..)
    and "graph" (based on graph classification - warning, CPU intensive)
    prior to rxn class, the mech is also subdivided into PESs
    """
    SPC_NAME = os.path.join(CWD, 'data', 'C10H10_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'C10H10_Pdep_mech.dat')
    MECH_REST_NAME = os.path.join(CWD_RESULTS, 'C10H10_mech_rest.txt')
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'C10H10_test_sortby_rxnclass.txt')
    ISOLATE_SPECIES = []
    # HEADER INDICATING THE REACTION CLASS
    SORT_STR = ['pes', 'rxn_class_broad', 'rxn_class_graph', 1]

    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)


def test__sortby_species_subpes():
    """
    select a species subset from a mechanism and extract all reactions they are involved to
    Within the reaction subset, classify according to subpes (or potentially any other criteria)
    """
    SPC_NAME = os.path.join(CWD, 'data', 'LLNL_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    MECH_REST_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_IC8_SPECIES_mech_rest.txt')
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_species_subpes_IC8.txt')
    ISOLATE_SPECIES = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
    SORT_STR = ['species', 'subpes', 1]  # HEADER INDICATING THE SPECIES SUBSET
    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)
    # NOW ORDER JUST BY SUBPES
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_subpes_IC8.txt')
    SORT_STR = ['subpes', 0]  # NO HEADER
    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)


def test__sortby_submech():
    """
    sort by fuel submechanism: extract reactions of
    fuel, fuel radicals, R+O2, R+O4
    then order by subpes
    """
    SPC_NAME = os.path.join(CWD, 'data', 'LLNL_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    MECH_REST_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_IC8_submech_mech_rest.txt')
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_submech_IC8.txt')
    ISOLATE_SPECIES = ['IC8']
    SORT_STR = ['submech', 1]  # HEADER INDICATING THE SPECIES SUBSET
    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)
    # NOW ORDER JUST BY SUBPES
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_submech_subpes_broadclass_IC8.txt')
    SORT_STR = ['submech', 'subpes','rxn_class_broad', 1]  # NO HEADER
    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)

def test__sortby_submech_class():
    """
    sort by fuel submechanism: extract reactions of
    fuel, fuel radicals, R+O2, R+O4
    then order by subpes
    """
    SPC_NAME = os.path.join(CWD, 'data', 'LLNL_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    MECH_REST_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_IC8_submech_class_mech_rest.txt')
    SORTMECH_NAME = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_submech_class_IC8.txt')
    ISOLATE_SPECIES = ['IC8', 'submech']
    SORT_STR = ['submech', 'rxn_class_broad','rxn_class_graph', 1]  # NO HEADER
    sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME,
              MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR)


def sort_main(CWD, SPC_NAME, MECH_NAME, SORTMECH_NAME, MECH_REST_NAME, ISOLATE_SPECIES, SORT_STR):
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
    spc_dct = mechanalyzer.parser.spc.order_species_by_atomcount(spc_dct)
    chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=elem_tuple, spc_dct=spc_dct, rxn_param_dct=rxn_param_dct_sorted,
        comments=cmts_dct)

    # IF YOU DERIVED A MECH SUBSET: SAVE ALL THE REMAINING REACTIONS IN ANOTHER FILE
    if ISOLATE_SPECIES:
        rxn_param_dct_rest = mechanalyzer.parser.util.filter_keys(
            rxn_param_dct, rxn_param_dct_sorted)
        mechanalyzer.parser.mech.write_mech(
            elem_tuple, spc_dct_full, rxn_param_dct_rest, MECH_REST_NAME)


if __name__ == '__main__':
 #   test__readwrite_thirdbody()
 #   test__sortby_species_subpes()
    test__sortby_submech()
 #   test__sortby_mult()
 #   test__sortby_molec_R1()
 #   test__sortby_rxnclass()
#    test__sortby_submech_class()
