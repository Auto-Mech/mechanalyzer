""" test mechanalyzer.parser.sort for different mechanisms in 'data/'
    using different sorting options
"""

import os
import tempfile
import numpy as np
from mechanalyzer.builder import sorter
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import spc as sparser


# Set Paths to test/data directory and output directory
CWD = os.path.dirname(os.path.realpath(__file__))
TMP_OUT = tempfile.mkdtemp()

# Set types for parsing mechanisms
SPC_TYPE = 'csv'
MECH_TYPE = 'chemkin'

# Test data
BIG_ARRAY = np.array([1e15, 1e15, 1e15])
MIDDLE_ARRAY = np.array([1e14, 1e14, 1e14])
LITTLE_ARRAY = np.array([1e13, 1e13, 1e13])

AL_KTP_DCT = {
    (('H2', 'O'), ('OH', 'H'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([0.157572885e+134, 2.79926202e+143, 1.72670689e+149])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([6.57572885e+134, 8.79926202e+143, 4.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}],
    (('H', 'O2'), ('OH', 'O'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([4.6548277231154764e+45,
                      8.556998184634325e+52, 4.662500917095324e+56])),
         1: (np.array([500, 1000, 1500]),
             np.array([4.6548277231154764e+45,
                       8.556998184634325e+52, 4.662500917095324e+56])),
         10: (np.array([500, 1000, 1500]),
              np.array([4.6548277231154764e+45,
                        8.556998184634325e+52, 4.662500917095324e+56]))}],
    (('H2', 'O'), ('OH', 'OH'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        None],
    (('H', 'O'), ('OH',), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([1.420849319576619e+96,
                      2.8405686431169553e+77, 3.4922934313599517e+72])),
         1: (np.array([500, 1000, 1500]),
             np.array([5.8295576381190475e+100,
                       2.3308958102634265e+82, 4.2985260083885116e+77])),
         10: (np.array([500, 1000, 1500]),
              np.array([2.3917907260059445e+105,
                        1.912671707993609e+87, 5.2908858341829314e+82]))}],
    (('H', 'O'), ('OH',), ('(+M)',)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([9.813202359645695e+109,
                      1.569488025355258e+92, 6.512342336821681e+87])),
         1: (np.array([500, 1000, 1500]),
             np.array([4.0262276922599165e+114,
                       1.2878805345625882e+97, 8.015784887656628e+92])),
         10: (np.array([500, 1000, 1500]),
              np.array([1.6519081983453455e+119,
                        1.0568008449314422e+102, 9.866312924289953e+97]))}],
    (('H2', 'O(S)'), ('OH', 'H'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
         1: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139]))}],
    (('H2', 'O2'), ('HO2V', 'H'), (None,)): [None, {
        'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
        1: (np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
        10: (np.array([500, 1000, 1500]),
             np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}]}


def __sort_with_input():
    """ sort by using the auxlilary input files to specify parameters
    """

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'NUIG_species.csv')
    mech_path = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    sort_path = os.path.join(CWD, 'data', 'sort.dat')

    spc_str, mech_str, sort_str = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism
    isolate_spc, sort_lst = mparser.parse_sort(sort_str)

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test__readwrite_thirdbody():
    """ test mechanalyzer.parser.sort

        Checks read/write of a small set of rxns involving third bodies
    """

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'NUIG_species.csv')
    mech_path = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by PES - No Headers Included
    isolate_spc = []
    sort_lst = ['pes', 0]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test__sortby_mult():
    """ test mechanalyzer.parser.sort

        Sort by multiplicity of the reaction
    """

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'C10H10_species.csv')
    mech_path = os.path.join(CWD, 'data', 'C10H10_HP_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by multiplicity - No Headers Included
    isolate_spc = []
    sort_lst = ['mult', 0]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test__sortby_molec_r1():
    """ test mechanalyzer.parser.sort

        Sort by first (heavier) reactant and molecularity of the reaction
    """

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'C10H10_species.csv')
    mech_path = os.path.join(CWD, 'data', 'C10H10_HP_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by R1-molecularity - No Headers Included
    isolate_spc = []
    sort_lst = ['r1', 'molecularity', 0]  # NO HEADERS INCLUDED

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test_sortby_rxnclass():
    """ test mechanalyzer.parser.sort

        sort by reaction class:
            both "broad" (based on multiplicity, reactants/products..)
        and "graph" (based on graph classification - warning, CPU intensive)
        prior to rxn class, the mech is also subdivided into PESs
    """

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_classes.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by reaction class
    isolate_spc = []
    sort_lst = ['pes', 'rxn_class_broad', 'rxn_class_graph', 1]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test__sortby_species_subpes():
    """ test mechanalyzer.parser.sort

        Select a species subset from a mechanism and
        extract all reactions they are involved to
        Within the reaction subset, classify according
        to subpes (or potentially any other criteria)
    """

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by subpes with Headers for specues subset
    isolate_spc = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
    sort_lst = ['species', 'subpes', 1]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)

    # Now sort mechanism by subpes no header
    isolate_spc = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
    sort_lst = ['subpes', 0]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test__sortby_submech():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes
    """

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['IC8']
    sort_lst = ['submech', 1]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)

    # Sort again
    isolate_spc = ['IC8']
    sort_lst = ['submech', 'subpes', 'rxn_class_broad', 1]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test_sortby_submech_class():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes
    """

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism
    isolate_spc = ['IC8']
    sort_lst = ['submech', 'rxn_class_broad', 'rxn_class_graph', 1]

    param_dct_sort, _, _, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    print(param_dct_sort)


def test_sort_ktp():
    """ test mechanalyzer.parser.sort

        sort ktp dictionary according to highest rate values/ratios
    """

    # Read mechanism files into strings
    spc_paths = [
        os.path.join(CWD, 'data', 'spc2.csv'),
        os.path.join(CWD, 'data', 'spc1B.csv')]
    mech_path = None
    sort_path = None

    spc_str, _, _ = _read_files(spc_paths[1], mech_path, sort_path)

    # Build spc and mech information
    spc_dct_full = sparser.build_spc_dct(spc_str, SPC_TYPE)
    mech_info = mparser.mech_info(AL_KTP_DCT, spc_dct_full)

    # Sort the mechanism
    isolate_spc = []
    sort_lst = ['molecularity', 'rxn_max_vals',
                'rxn_max_ratio', 'rxn_class_broad', 0]

    srt_mch = sorter.sorting(
        mech_info, spc_dct_full, sort_lst, isolate_spc)
    sorted_idx, _, _ = srt_mch.return_mech_df()
    al_ktp_dct_sorted = sorter.reordered_mech(AL_KTP_DCT, sorted_idx)

    print(al_ktp_dct_sorted)


def test__build_sorted_pesdct():
    """ sort by subpes for a subset of species and
        get the corresponding subpes dictionaries
    """

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism and build pes
    isolate_spc = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
    sort_lst = ['subpes', 0]

    pes_dct = sorter.sorted_pes_dct(
        spc_str, mech_str, isolate_spc, sort_lst)

    print(pes_dct)


# Helper function
def _read_files(spc_path, mech_path, sort_path):
    """ read file names
    """

    spc_str, mech_str, sort_str = '', '', ''

    if spc_path is not None:
        with open(spc_path) as fobj:
            spc_str = fobj.read()
    if mech_path is not None:
        with open(mech_path) as fobj:
            mech_str = fobj.read()
    if sort_path is not None:
        with open(sort_path) as fobj:
            sort_str = fobj.read()

    return spc_str, mech_str, sort_str
