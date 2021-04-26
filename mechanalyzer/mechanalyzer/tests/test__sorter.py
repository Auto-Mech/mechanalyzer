""" test mechanalyzer.parser.sort for different mechanisms in 'data/'
    using different sorting options
"""

import os
import automol.util.dict_
import chemkin_io
import mechanalyzer


# Set Paths to test/data directory and output directory
CWD = os.path.dirname(os.path.realpath(__file__))
CWD_RESULTS = os.path.join(CWD, 'sorter_results')
if not os.path.isdir(CWD_RESULTS):
    os.makedirs(CWD_RESULTS)

# Set types for parsing mechanisms
SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'


def test__readwrite_thirdbody():
    """ test mechanalyzer.parser.sort

        Checks read/write of a small set of rxns involving third bodies
    """

    spc_name = os.path.join(CWD, 'data', 'NUIG_species.csv')
    mech_name = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    mech_rest_name = os.path.join(CWD_RESULTS, 'NUIG_mech_rest.txt')
    sortmech_name = os.path.join(
        CWD_RESULTS, 'NUIG_test_readwrite_thirdbody.txt')
    isolate_species = []
    sort_str = ['pes', 0]  # ARRANGE BY PES- NO HEADERS INCLUDED

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def test__sortby_mult():
    """ test mechanalyzer.parser.sort

        Sort by multiplicity of the reaction
    """

    spc_name = os.path.join(CWD, 'data', 'C10H10_species.csv')
    mech_name = os.path.join(CWD, 'data', 'C10H10_HP_mech.dat')
    mech_rest_name = os.path.join(CWD_RESULTS, 'C10H10_mech_rest.txt')
    sortmech_name = os.path.join(CWD_RESULTS, 'C10H10_test_sortby_mult.txt')
    isolate_species = []
    sort_str = ['mult', 0]  # NO HEADERS INCLUDED

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def test__sortby_molec_r1():
    """ test mechanalyzer.parser.sort

        Sort by first (heavier) reactant and molecularity of the reaction
    """

    spc_name = os.path.join(CWD, 'data', 'C10H10_species.csv')
    mech_name = os.path.join(CWD, 'data', 'C10H10_HP_mech.dat')
    mech_rest_name = os.path.join(CWD_RESULTS, 'C10H10_mech_rest.txt')
    sortmech_name = os.path.join(
        CWD_RESULTS, 'C10H10_test_sortby_molec_R1.txt')
    isolate_species = []
    sort_str = ['r1', 'molecularity', 0]  # NO HEADERS INCLUDED

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def test__sortby_rxnclass():
    """ test mechanalyzer.parser.sort

        sort by reaction class:
            both "broad" (based on multiplicity, reactants/products..)
        and "graph" (based on graph classification - warning, CPU intensive)
        prior to rxn class, the mech is also subdivided into PESs
    """

    spc_name = os.path.join(CWD, 'data', 'C10H10_species.csv')
    mech_name = os.path.join(CWD, 'data', 'C10H10_Pdep_mech.dat')
    mech_rest_name = os.path.join(CWD_RESULTS, 'C10H10_mech_rest.txt')
    sortmech_name = os.path.join(
        CWD_RESULTS, 'C10H10_test_sortby_rxnclass.txt')
    isolate_species = []

    # HEADER INDICATING THE REACTION CLASS
    sort_str = ['pes', 'rxn_class_broad', 'rxn_class_graph', 1]

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def test__sortby_species_subpes():
    """ test mechanalyzer.parser.sort

        Select a species subset from a mechanism and
        extract all reactions they are involved to
        Within the reaction subset, classify according
        to subpes (or potentially any other criteria)
    """

    spc_name = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_name = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    mech_rest_name = os.path.join(
        CWD_RESULTS, 'LLNL_IC8_SPECIES_mech_rest.txt')
    sortmech_name = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_species_subpes_IC8.txt')
    isolate_species = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
    sort_str = ['species', 'subpes', 1]  # HEADER INDICATING THE SPECIES SUBSET

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)

    # NOW ORDER JUST BY SUBPES
    sortmech_name = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_subpes_IC8.txt')
    sort_str = ['subpes', 0]  # NO HEADER
    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def test__sortby_submech():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes
    """

    spc_name = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_name = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    mech_rest_name = os.path.join(
        CWD_RESULTS, 'LLNL_IC8_submech_mech_rest.txt')

    sortmech_name = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_submech_IC8.txt')
    isolate_species = ['IC8']
    sort_str = ['submech', 1]  # HEADER INDICATING THE SPECIES SUBSET

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)

    # NOW ORDER JUST BY SUBPES
    sortmech_name = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_submech_subpes_broadclass_IC8.txt')
    # NO HEADER
    sort_str = ['submech', 'subpes', 'rxn_class_broad', 1]

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def test__sortby_submech_class():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes
    """

    spc_name = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_name = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    mech_rest_name = os.path.join(
        CWD_RESULTS, 'LLNL_IC8_submech_class_mech_rest.txt')

    sortmech_name = os.path.join(
        CWD_RESULTS, 'LLNL_test_sortby_submech_class_IC8.txt')
    isolate_species = ['IC8', 'submech']
    # NO HEADER
    sort_str = ['submech', 'rxn_class_broad', 'rxn_class_graph', 1]

    _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str)


def _sort_main(spc_name, mech_name, sortmech_name,
               mech_rest_name, isolate_species, sort_str):
    """ Function that conducts the sorting process for all of the above tests
    """

    # Read the files
    with open(spc_name, 'r') as spc_obj:
        spc_str = spc_obj.read()
    with open(mech_name, 'r') as mech_obj:
        mech_str = mech_obj.read()

    # Build mech information
    spc_dct_full = mechanalyzer.parser.spc.build_spc_dct(spc_str, SPC_TYPE)
    rxn_param_dct, mech_info, elems = mechanalyzer.parser.mech.parse_mechanism(
        mech_str, MECH_TYPE, spc_dct_full)

    # Sorting: sort the mech and build the sorted rxn param dct
    # sorted_idx, cmts_dct, spc_dct = mechanalyzer.parser.mech.sort_mechanism(
    #    mech_info, spc_dct_full, sort_str, isolate_species)
    srt_mch = mechanalyzer.parser.mech.sorting(
        mech_info, spc_dct_full, sort_str, isolate_species)
    sorted_idx, cmts_dct, spc_dct = mechanalyzer.parser.mech.sorted_mech(
        srt_mch)
    rxn_param_dct_sorted = mechanalyzer.parser.mech.reordered_mech(
        rxn_param_dct, sorted_idx)

    # Write the new mechanism
    spc_dct_ord = mechanalyzer.parser.spc.order_species_by_atomcount(spc_dct)
    mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
        elem_tuple=elems, spc_dct=spc_dct_ord,
        rxn_param_dct=rxn_param_dct_sorted,
        comments=cmts_dct)
    with open(sortmech_name, 'w') as mech1_obj:
        mech1_obj.write(mech_str)

    # If isolated species provided, save remaining reactions in another file
    if isolate_species:
        rxn_param_dct_rest = automol.util.dict_.filter_keys(
            rxn_param_dct, rxn_param_dct_sorted)
        mech_rest_str = chemkin_io.writer.mechanism.write_chemkin_file(
            elem_tuple=elems, spc_dct=spc_dct_full,
            rxn_param_dct=rxn_param_dct_rest,
            comments=cmts_dct)
        with open(mech_rest_name, 'w') as mech2_obj:
            mech2_obj.write(mech_rest_str)


if __name__ == '__main__':
    test__readwrite_thirdbody()
    test__sortby_species_subpes()
    test__sortby_submech()
    test__sortby_mult()
    test__sortby_molec_r1()
    test__sortby_rxnclass()
    test__sortby_submech_class()
