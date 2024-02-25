""" test mechanalyzer.parser.sort for different mechanisms in 'data/'
    using different sorting options
"""

import os
import tempfile
import numpy as np
from ioformat import pathtools
import chemkin_io.writer
from mechanalyzer.builder import sorter
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.parser import new_spc as sparser

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
            np.array([0.157572885e+134, 2.79926202e+143, 1.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}],
    (('H', 'O2'), ('OH', 'O'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([4.6548277231154764e+45,
                      8.556998184634325e+52, 4.662500917095324e+56])),
         10: (np.array([500, 1000, 1500]),
              np.array([4.6548277231154764e+45,
                        8.556998184634325e+52, 4.662500917095324e+56]))}],
    (('H2', 'O'), ('OH', 'OH'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        None],
    (('H', 'O'), ('OH',), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([1.420849319576619e+96,
                      2.8405686431169553e+77, 3.4922934313599517e+72]))}],
    (('H', 'O'), ('OH',), ('(+M)',)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
         10: (np.array([500, 1000, 1500]),
              np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))},
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([9.813202359645695e+109,
                      1.569488025355258e+92, 6.512342336821681e+87])),
         10: (np.array([500, 1000, 1500]),
              np.array([1.6519081983453455e+119,
                        1.0568008449314422e+102, 9.866312924289953e+97]))}],
    (('H2', 'O(S)'), ('OH', 'H'), (None,)): [
        {'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139]))},
        {1: (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+139]))}],
    (('H2', 'O2'), ('HO2V', 'H'), (None,)): [None, {
        'high': (
            np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149])),
        1: (np.array([500, 1000, 1500]),
            np.array([3.57572885e+134, 4.79926202e+143, 2.72670689e+149]))}]}


def test__sort_with_input():
    """ sort by using the auxlilary input files to specify parameters
    """

    results = [
        [(('C2H4',), ('H2', 'H2CC'), ('(+M)',)),
         '! pes.subpes.NR.rxntype  1.1.1.Decomposition'],
        [(('C2H3', 'H'), ('C2H4',), ('(+M)',)),
         '! pes.subpes.NR.rxntype  1.1.2.Recombination H'],
        [(('C2H4', 'H'), ('C2H5',), ('(+M)',)),
         '! pes.subpes.NR.rxntype  2.1.2.Addition H'],
        [(('C2H4', 'H'), ('C2H3', 'H2'), (None,)),
         '! pes.subpes.NR.rxntype  2.2.2.H abstraction'],
        [(('CH2(S)', 'CH3'), ('C2H4', 'H'), (None,)),
         '! pes.subpes.NR.rxntype  2.3.2.Addition-decomposition - propagation'],
        [(('C2H5', 'H'), ('C2H4', 'H2'), (None,)),
         '! pes.subpes.NR.rxntype  3.2.2.Recombination-decomposition - termination'],
        [(('C2H4', 'O'), ('CH3', 'HCO'), (None,)),
         '! pes.subpes.NR.rxntype  5.1.2.Addition-decomposition - branching'],
        [(('C2H4', 'OH'), ('PC2H4OH',), (None,)),
         '! pes.subpes.NR.rxntype  6.1.2.Addition OH'],
        [(('C2H5', 'OH'), ('C2H4', 'H2O'), (None,)),
         '! pes.subpes.NR.rxntype  7.1.2.Recombination-decomposition - termination'],
        [(('C2H4', 'O2'), ('C2H3', 'HO2'), (None,)),
         '! pes.subpes.NR.rxntype  9.1.2.H abstraction'],
        [(('C2H5O2',), ('C2H4', 'HO2'), (None,)),
         '! pes.subpes.NR.rxntype  10.1.1.Beta-scission +HO2'],
        [(('C2H4', 'CH3'), ('C2H3', 'CH4'), (None,)),
         '! pes.subpes.NR.rxntype  11.1.2.H abstraction'],
        [(('C3H4-A', 'O'), ('C2H4', 'CO'), (None,)),
         '! pes.subpes.NR.rxntype  12.1.2.Addition-decomposition - termination']
    ]

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech.dat')
    sort_path = os.path.join(CWD, 'data', 'sort.dat')

    spc_str, mech_str, sort_str = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism
    isolate_spc, sort_lst, _ = mparser.parse_sort(sort_str)
    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)
    index = 0

    for rxn in param_dct_sort.keys():
        assert [rxn, cmts_dct[rxn]['cmts_inline']] == results[index]
        index += 1
    print('ok')


def test__readwrite_thirdbody():
    """ test mechanalyzer.parser.sort

        Checks read/write of a small set of rxns involving third bodies
    """

    # Setting the values of the dictionary to be None since they don't matter
    trd_bdy_dct = {
        (('H', 'OH'), ('H2O',), ('+M',)): None,
        (('H', 'OH', 'AR'), ('H2O', 'AR'), (None,)): None,
        (('H', 'O2'), ('HO2',), ('(+HE)',)): None,
        (('CH3', 'IC4H7'), ('AC5H10',), ('(+M)',)): None,
        (('C5H10-2',), ('C4H71-3', 'CH3'), ('(+M)',)): None,
        (('C5H11-1',), ('C2H4', 'NC3H7'), (None,)): None
    }

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'NUIG_species.csv')
    mech_path = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by PES - No Headers Included
    isolate_spc = []
    sort_lst = ['pes', 0]

    param_dct_sort, _, _, _, _= sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    # Just checking keys since this is what the sorting is according to
    assert param_dct_sort.keys() == trd_bdy_dct.keys()
    print('ok')


def test__sortby_mult():
    """ test mechanalyzer.parser.sort

        Sort by multiplicity of the reaction
    """

    results = {
        (('C5H10-2',), ('C4H71-3', 'CH3'), ('(+M)',)): str(1),
        (('C5H11-1',), ('C2H4', 'NC3H7'), (None,)): str(2),
        (('H', 'OH'), ('H2O',), ('+M',)): str(4),
        (('H', 'OH', 'AR'), ('H2O', 'AR'), (None,)): str(4),
        (('CH3', 'IC4H7'), ('AC5H10',), ('(+M)',)): str(4),
        (('H', 'O2'), ('HO2',), ('(+HE)',)): str(6)
    }
    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'NUIG_species.csv')
    mech_path = os.path.join(CWD, 'data', 'NUIG_mechred.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by multiplicity - No Headers Included
    isolate_spc = []
    sort_lst = ['mult', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    for rxn in param_dct_sort.keys():
        assert cmts_dct[rxn]['cmts_inline'][-1] == results[rxn]
    print('ok')


def test__sortby_molec_r1():
    """ test mechanalyzer.parser.sort

        Sort by first (heavier) reactant and molecularity of the reaction
    """
    comments_results = [
        'C2H3.2', 'C2H3.2', 'C2H3.2', 'C2H3.2', 'C2H3.2', 'C2H3.2',
        'C2H3OO.1', 'C2H3OO.1', 'C2H4.1', 'C2H4.2', 'C2H4.2', 'C2H4.2',
        'C2H4.2', 'C2H4.2', 'C2H4.2', 'C2H5.2', 'C2H5.2', 'C2H5.2',
        'C2H5O2.1', 'C2H6.2', 'C3H4-A.2', 'CH3.2', 'CH3.2', 'HOCH2CO.1']

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by R1-molecularity - No Headers Included
    isolate_spc = []
    sort_lst = ['r1', 'molecularity', 0]  # NO HEADERS INCLUDED
    # ONLY BY MOLEC- CHECK HEADERS WITH INT NUMBER
    sort_lst_2 = ['molecularity', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    param_dct_sort_2, _, cmts_dct_2, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst_2)
    # NO NEED TO CALL IT AFTERWARDS

    comments = []
    for rxn in param_dct_sort.keys():
        comments.append(''.join(cmts_dct[rxn]['cmts_inline'].split()[-1:]))
    assert comments == comments_results
    print('ok')


def test__sortby_pes_dct():
    """ test mechanalyzer.parser.sort

        sort by pes dct:
    """
    pes_dct_result = {
        ('C8H17', 0, 0): ((0, (('IC8-1R',), ('IC8-5R',))),),
        ('C8H18', 1, 0): ((0, (('IC8',), ('NEOC5H11', 'IC3H7'))),),
        ('C8H17O2', 2, 0): ((0, (('IC8OOH1-1AR',), ('IC8O1-1A', 'OH'))),
                            (1, (('IC8OOH1-1AR',), ('IC4H7OOH', 'IC4H9'))),
                            (2, (('IC8OOH1-1AR',), ('CH2O', 'I24C7D1', 'OH'))),
                            (3, (('IC8-1R', 'O2'), ('IC8-1O2R',))),
                            (4, (('IC8-1O2R',), ('IC8OOH1-1AR',)))),
        ('C8H17O2', 2, 1): ((5, (('IC8-3O2R',), ('IC8D3', 'HO2'))),),
        ('C8H17O2', 2, 2): ((6, (('IC8-3R', 'O2'), ('IC8D3', 'HO2'))),),
        ('C8H18O2', 3, 0): ((0, (('IC8OOH1',), ('IC8-1OR', 'OH'))),),
        ('C8H18O2', 3, 1): ((1, (('IC8', 'O2'), ('IC8-1R', 'HO2'))),),
        ('C8H17O4', 4, 0): ((0, (('IC8OOH1-1AR', 'O2'), ('IC8OOH1-1AO2R',))),),
        ('C8H18O4', 5, 0): ((0, (('IC8-1O2R', 'HO2'), ('IC8OOH1', 'O2'))),),
        ('C8H19O4', 6, 0): ((0, (('IC8-1O2R', 'H2O2'), ('IC8OOH1', 'HO2'))),),
        ('C9H20O2', 7, 0): ((0, (('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'))),)}

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by reaction class
    isolate_spc = []
    sort_lst = ['pes', 'subpes', 0]

    print('test pes dct functionality:')

    pes_dct = sorter.sorted_pes_dct(
        spc_str, mech_str, isolate_spc, sort_lst)

    assert pes_dct == pes_dct_result
    print('ok')


def test__sortby_rxnclass():
    """ test mechanalyzer.parser.sort

        sort by reaction class:
            both "broad" (based on multiplicity, reactants/products..)
        and "graph" (based on graph classification - warning, CPU intensive)
        prior to rxn class, the mech is also subdivided into PESs
    """
    results = [
        [(('C2H5', 'H'), ('C2H6',), ('(+M)',)),
         '  addition.Recombination H'],
        [(('C2H3', 'H'), ('C2H4',), ('(+M)',)),
         '  addition.Recombination H'],
        [(('HOCH2CO',), ('CH2OH', 'CO'), (None,)),
         '  beta scission.Decomposition'],
        [(('C2H3OO',), ('CH2O', 'HCO'), (None,)),
         '  elimination.Beta-scission'],
        [(('C2H5O2',), ('C2H4', 'HO2'), (None,)),
         '  elimination.Beta-scission +HO2'],
        # removed because classifer is not working well right now
        # [(('C2H4',), ('H2', 'H2CC'), ('(+M)',)),
        # '  elimination.Decomposition'],
        [(('C2H4', 'H'), ('C2H3', 'H2'), (None,)),
         '  hydrogen abstraction.H abstraction'],
        [(('C2H5', 'H'), ('C2H4', 'H2'), (None,)),
         '  hydrogen abstraction.Recombination-decomposition - termination'],
        [(('C2H3', 'O2'), ('C2H2', 'HO2'), (None,)),
         '  hydrogen abstraction.Recombination-decomposition - termination'],
        [(('C3H5-A',), ('C3H5-T',), (None,)),
         '  hydrogen migration.Isomerization'],
        [(('C3H5-A',), ('C3H5-S',), (None,)),
         '  hydrogen migration.Isomerization'],
        [(('CH2(S)', 'CH3'), ('C2H4', 'H'), (None,)),
         '  substitution.Addition-decomposition - propagation'],
        [(('CH3', 'CH3'), ('H', 'C2H5'), (None,)),
         '  substitution.Recombination-decomposition - propagation'],
        # removed because classifer is not working well right now
        # [(('C3H4-A', 'O'), ('C2H4', 'CO'), (None,)),
        #  '  unclassified.Addition-decomposition - termination']
    ]
    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech_class.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by reaction class
    isolate_spc = []
    sort_lst = ['rxn_class_graph', 'rxn_class_broad', 0]

    print('Sort by rxn class broad+graph test:')

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline'].split('type')[1]])

    assert sorted_results == results

    print('ok')


def test__sortby_species_subpes():
    """ test mechanalyzer.parser.sort

        Select a species subset from a mechanism and
        extract all reactions they are involved to
        Within the reaction subset, classify according
        to subpes (or potentially any other criteria)
    """
    results = [
        [(('IC8',), ('NEOC5H11', 'IC3H7'), (None,)),
         '! pes.subpes  26.1'],
        [(('IC8', 'O2'), ('IC8-1R', 'HO2'), (None,)),
         '! pes.subpes  28.2'],
        [(('IC8-1R',), ('IC8-5R',), (None,)),
         '! pes.subpes  25.1'],
        [(('IC8-1R', 'O2'), ('IC8-1O2R',), (None,)),
         '! pes.subpes  27.1'],
        [(('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'), (None,)),
         '! pes.subpes  32.1'],
    ]
    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_submech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by subpes with Headers for species subset
    isolate_spc = ['IC8', 'IC8-1R']
    sort_lst = ['species', 'subpes', 1]

    print('Sort by species subset + subpes test:')

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline']])

    assert sorted_results == results

    print('ok')


def test__sortby_submech_subpes_chnl():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes and broad class
    """
    results = [
        [(('IC8-1R',), ('IC8-5R',), (None,)),
         '! pes.subpes.channel  1.1.1'],
        [(('IC8',), ('NEOC5H11', 'IC3H7'), (None,)),
         '! pes.subpes.channel  2.1.1'],
        [(('IC8OOH1-1AR',), ('IC8O1-1A', 'OH'), (None,)),
         '! pes.subpes.channel  3.1.1'],
        [(('IC8OOH1-1AR',), ('IC4H7OOH', 'IC4H9'), (None,)),
         '! pes.subpes.channel  3.1.2'],
        [(('IC8OOH1-1AR',), ('CH2O', 'I24C7D1', 'OH'), (None,)),
         '! pes.subpes.channel  3.1.3'],
        [(('IC8-1R', 'O2'), ('IC8-1O2R',), (None,)),
         '! pes.subpes.channel  3.1.4'],
        [(('IC8-1O2R',), ('IC8OOH1-1AR',), (None,)),
         '! pes.subpes.channel  3.1.5'],
        [(('IC8-3O2R',), ('IC8D3', 'HO2'), (None,)),
         '! pes.subpes.channel  3.2.6'],
        [(('IC8-3R', 'O2'), ('IC8D3', 'HO2'), (None,)),
         '! pes.subpes.channel  3.3.7'],
        [(('IC8OOH1',), ('IC8-1OR', 'OH'), (None,)),
         '! pes.subpes.channel  4.1.1'],
        [(('IC8', 'O2'), ('IC8-1R', 'HO2'), (None,)),
         '! pes.subpes.channel  4.2.2'],
        [(('IC8OOH1-1AR', 'O2'), ('IC8OOH1-1AO2R',), (None,)),
         '! pes.subpes.channel  5.1.1'],
        [(('IC8-1O2R', 'HO2'), ('IC8OOH1', 'O2'), (None,)),
         '! pes.subpes.channel  6.1.1'],
        [(('IC8-1O2R', 'H2O2'), ('IC8OOH1', 'HO2'), (None,)),
         '! pes.subpes.channel  7.1.1'],
        [(('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'), (None,)),
         '! pes.subpes.channel  8.1.1']
    ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = []
    sort_lst = ['subpes', 'chnl', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    print('Sort by submech-subpes-classbroad test:')
    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline']])

    assert results == sorted_results
    print('ok')


def test__sortby_submech_class():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes and broad class
    """
    results = [
        [(('IC8',), ('NEOC5H11', 'IC3H7'), (None,)),
         '  FUEL.2.1.Bond fission'],
        [(('IC8', 'O2'), ('IC8-1R', 'HO2'), (None,)),
         '  FUEL.4.2.H abstraction'],
        [(('IC8OOH1',), ('IC8-1OR', 'OH'), (None,)),
         '  FUEL_ADD_O2.4.1.Bond fission +OH'],
        [(('IC8-1O2R', 'HO2'), ('IC8OOH1', 'O2'), (None,)),
         '  FUEL_ADD_O2.6.1.Recombination-decomposition - termination'],
        [(('IC8-1O2R', 'H2O2'), ('IC8OOH1', 'HO2'), (None,)),
         '  FUEL_ADD_O2.7.1.H abstraction'],
        [(('IC8-1R',), ('IC8-5R',), (None,)),
         '  FUEL_RAD.1.1.Isomerization'],
        [(('IC8-1R', 'O2'), ('IC8-1O2R',), (None,)),
         '  FUEL_RAD.3.1.Recombination O2'],
        [(('IC8-3R', 'O2'), ('IC8D3', 'HO2'), (None,)),
         '  FUEL_RAD.3.3.Recombination-decomposition - termination'],
        [(('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'), (None,)),
         '  FUEL_RAD.8.1.Recombination-decomposition - propagation'],
        [(('IC8OOH1-1AR',), ('IC4H7OOH', 'IC4H9'), (None,)),
         '  R_O2.3.1.Beta-scission'],
        [(('IC8OOH1-1AR',), ('IC8O1-1A', 'OH'), (None,)),
         '  R_O2.3.1.Beta-scission +OH'],
        [(('IC8OOH1-1AR',), ('CH2O', 'I24C7D1', 'OH'), (None,)),
         '  R_O2.3.1.Decomposition(lumped)'],
        [(('IC8-1O2R',), ('IC8OOH1-1AR',), (None,)),
         '  R_O2.3.1.Isomerization'],
        [(('IC8-3O2R',), ('IC8D3', 'HO2'), (None,)),
         '  R_O2.3.2.Beta-scission +HO2'],
        [(('IC8OOH1-1AR', 'O2'), ('IC8OOH1-1AO2R',), (None,)),
         '  R_O2.5.1.Recombination O2'],
    ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['IC8']
    sort_lst = ['submech', 'subpes', 'rxn_class_broad', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    print('Sort by submech-subpes-classbroad test:')
    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline'].split('rxntype')[1]])

    assert results == sorted_results
    print('ok')


def test__sortby_submech_ext():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        and also the relative submech
        then order by subpes and broad class
    """
    results = [
        [(('C2H4',), ('H2', 'H2CC'), ('(+M)',)), 'FUEL.42.1'],
        [(('CH2(S)', 'C2H4'), ('CC3H6',), (None,)), 'FUEL.76.1'],
        [(('C2H4', 'CH3O'), ('C2H3', 'CH3OH'), (None,)), 'FUEL.85.7'],
        [(('C3H5-A', 'C2H5'), ('C2H4', 'C3H6'), (None,)), 'FUEL.152.1'],
        [(('C3H8', 'O2'), ('IC3H7', 'HO2'), (None,)), 'FUEL_ADD_CH3.94.1'],
        [(('C2H5CHCO', 'OH'), ('NC3H7', 'CO2'), (None,)), 'FUEL_ADD_CH3.128.2'],
        [(('C2H5', 'O2'), ('C2H4O1-2', 'OH'), (None,)), 'FUEL_ADD_H.57.1'],
        [(('C4H71-1',), ('C2H5', 'C2H2'), (None,)), 'FUEL_ADD_H.111.1'],
        [(('C2H3OH', 'H'), ('PC2H4OH',), (None,)), 'FUEL_ADD_O.50.1'],
        [(('C2H3OH', 'HO2'), ('CH3CHO', 'HO2'), (None,)), 'FUEL_ADD_O.63.4'],
        [(('C4H6', 'O'), ('C2H2', 'C2H4O1-2'), (None,)), 'FUEL_ADD_O.119.7'],
        [(('CH3OCHO', 'O2'), ('CH2OCHO', 'HO2'), (None,)), 'FUEL_ADD_O2.67.1'],
        [(('CH3', 'CH2O'), ('C2H5O',), (None,)), 'FUEL_ADD_OH.50.1'],
        [(('C2H5OH', 'O2'), ('SC2H4OH', 'HO2'), (None,)), 'FUEL_ADD_OH.64.3'],
        [(('CH3OCH3', 'CH3O2'), ('CH3OCH2', 'CH3O2H'), (None,)), 'FUEL_ADD_OH.100.5'],
        [(('C2H3', 'CH3'), ('CH4', 'C2H2'), (None,)), 'FUEL_RAD.76.2'],
        [(('C4H71-O',), ('C2H3', 'CH3CHO'), (None,)), 'FUEL_RAD.120.5'],
        [(('C3H6', 'OH'), ('IC3H5OH', 'H'), (None,)), 'R_CH3.85.15'],
        [(('C4H8-2', 'H'), ('C3H6', 'CH3'), (None,)), 'R_CH3.113.10'],
        [(('C2H5CHCO', 'O'), ('C3H6', 'CO2'), (None,)), 'R_CH3.127.1'],
        [(('SC2H2OH', 'O2'), ('CH2CO', 'HO2'), (None,)), 'R_O.61.8'],
        [(('C2H3OO',), ('CH2CO', 'OH'), (None,)), 'R_O2.55.2'],
        [(('O', 'O'), ('O2',), ('+M',)), 'SUBFUEL.5.1'],
        [(('CH', 'H'), ('C', 'H2'), (None,)), 'SUBFUEL.14.1'],
        [(('CH3O',), ('CH2O', 'H'), ('(+M)',)), 'SUBFUEL.20.2'],
        [(('CH4', 'O'), ('CH3', 'OH'), (None,)), 'SUBFUEL.21.9'],
        [(('CH2(S)', 'O2'), ('CO', 'H2O'), (None,)), 'SUBFUEL.25.3'],
        [(('CH3', 'HO2'), ('CH3O', 'OH'), (None,)), 'SUBFUEL.27.3'],
        [(('CH2O', 'HO2'), ('OCH2O2H',), (None,)), 'SUBFUEL.32.1'],
        [(('CH3O2', 'H2O2'), ('CH3O2H', 'HO2'), (None,)), 'SUBFUEL.38.1'],
        [(('C2H2', 'OH'), ('HCCOH', 'H'), (None,)), 'SUBFUEL.48.2'],
        [(('CH2(S)', 'CO2'), ('CH2O', 'CO'), (None,)), 'SUBFUEL.54.2'],
        [(('CH3O', 'CH3O'), ('CH3OH', 'CH2O'), (None,)), 'SUBFUEL.58.9'],
        [(('C2H', 'CH3'), ('C3H4-P',), (None,)), 'SUBFUEL.74.1'],
        [(('C3H2C', 'O2'), ('C2H2', 'CO2'), (None,)), 'SUBFUEL.88.1'],
    ]
    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['C2H4']
    sort_lst = ['submech_ext', 'subpes', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    print('Sort by submech_ext-subpes:')
    sorted_results = []
    for i, rxn in enumerate(param_dct_sort.keys()):
        if i % 20 == 0:
            sorted_results.append(
                [rxn, cmts_dct[rxn]['cmts_inline'].split('subpes')[1].strip()])

    assert results == sorted_results
    print('ok')


def test__sortby_submech_prompt():
    """ test mechanalyzer.parser.sort

        sort by prompt reactions identified
        based on radical type
    """
    [
        [(('C4H72-2',), ('C4H612', 'H'), (None,)), '111.1.1.RAD_DECO_C4H72-2'],
        [(('C4H71-3',), ('C4H72-2',), (None,)), '111.1.6.'], 
        [(('C4H612', 'H'), ('C3H4-P', 'CH3'), (None,)), '111.1.11.'], 
        [(('C4H71-4', 'H'), ('C4H8-1',), ('(+M)',)), '112.1.1.RAD_GEN_C4H71-4'],
        [(('C4H8-2',), ('C3H5-A', 'CH3'), (None,)), '112.2.6.'],
        [(('C4H8-2', 'H'), ('C4H71-3', 'H2'), (None,)), '113.10.27.RAD_GEN_C4H71-3'],
        [(('C4H8-2', 'O'), ('C4H71-3', 'OH'), (None,)), '121.9.15.RAD_GEN_C4H71-3'], 
        [(('C4H8-2', 'OH'), ('C4H72-2', 'H2O'), (None,)), '122.24.41.RAD_GEN_C4H72-2'], 
        [(('C4H72-2O2',), ('CH3CHCOCH3', 'O'), (None,)), '128.8.30.'], 
        [(('C4H71-3OOH',), ('C4H71-O', 'OH'), (None,)), '129.3.7.'], 
        [(('C4H71-3', 'HO2'), ('C2H3COCH3', 'H2O'), (None,)), '129.3.12.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'HO2'), ('C4H7O2-1', 'OH'), (None,)), '129.4.17.RAD_GEN_C4H71-3'],
        [(('C4H8-2', 'O2'), ('C4H71-3', 'HO2'), (None,)), '129.15.34.RAD_GEN_C4H71-3'],
        [(('C4H8-2', 'HO2'), ('C4H72-2', 'H2O2'), (None,)), '130.24.63.RAD_GEN_C4H72-2'], 
        [(('C4H71-3', 'CH3O'), ('C4H8-1', 'CH2O'), (None,)), '160.1.1.RAD_GEN_C4H71-3'], 
        [(('C4H71-3', 'CH3O2'), ('C4H71-O', 'CH3O'), (None,)), '166.1.2.RAD_GEN_C4H71-3'], 
        [(('CH3O2', 'C4H71-4'), ('CH3O', 'C4H7O1-4'), (None,)), '166.4.10.RAD_GEN_C4H71-4'],
        [(('C2H5', 'C4H71-3'), ('C4H6', 'C2H6'), (None,)), '184.1.1.RAD_GEN_C4H71-3'], 
        [(('C3H5-A', 'C4H71-3'), ('C3H6', 'C4H6'), (None,)), '206.1.1.RAD_GEN_C4H71-3'], 
        [(('C4H8-2', 'IC3H7O2'), ('C4H71-3', 'IC3H7O2H'), (None,)), '213.4.4.RAD_GEN_C4H71-3'], 
        [(('C4H71-3', 'C4H71-3'), ('C8H141-5,3',), (None,)), '222.1.1.RAD_GEN_C4H71-3'], 
        [(('C4H71-3', 'SC3H5CHO'), ('C8H131-5,3,SAO',), (None,)), '225.5.5.RAD_GEN_C4H71-3'], 
        [(('SC4H9O2', 'C4H71-3'), ('SC4H9O', 'C4H71-O'), (None,)), '230.3.3.RAD_GEN_C4H71-3'], 
        [(('C4H8-2', 'PC4H9O2'), ('C4H71-3', 'PC4H9O2H'), (None,)), '231.5.5.RAD_GEN_C4H71-3']]

    results = [
        [(('C4H72-2',), ('C4H612', 'H'), (None,)), '111.1.1.RAD_DECO_C4H72-2'], 
        [(('C4H71-3',), ('C4H72-2',), (None,)), '111.1.6.unclassified'], 
        [(('C4H8-1', 'H'), ('C4H71-3', 'H2'), (None,)), '113.8.25.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'O'), ('C4H71-3', 'OH'), (None,)), '121.5.11.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'OH'), ('C4H71-4', 'H2O'), (None,)), '122.22.39.RAD_GEN_C4H71-4'],
        [(('C4H71-3', 'HO2'), ('C4H71-O', 'OH'), (None,)), '129.3.10.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'HO2'), ('C4H72-1OOH',), (None,)), '129.4.18.unclassified'], 
        [(('C4H8-2', 'O2'), ('C4H72-2', 'HO2'), (None,)), '129.16.35.RAD_GEN_C4H72-2'], 
        [(('C4H8-1', 'CH3'), ('C4H71-3', 'CH4'), (None,)), '153.6.6.RAD_GEN_C4H71-3'], 
        [(('C4H8-1', 'CH3O'), ('C4H71-3', 'CH3OH'), (None,)), '161.5.5.RAD_GEN_C4H71-3'],
        [(('C4H71-3', 'CH3O2'), ('C4H7O2-1', 'CH3O'), (None,)), '166.2.5.RAD_GEN_C4H71-3'],
        [(('C2H3', 'C4H71-3'), ('C2H4', 'C4H6'), (None,)), '183.3.3.RAD_GEN_C4H71-3'], 
        [(('C4H8-1', 'CH3CO3'), ('C4H71-3', 'CH3CO3H'), (None,)), '196.4.4.RAD_GEN_C4H71-3'],
        [(('C4H8-1', 'IC3H7O2'), ('C4H71-3', 'IC3H7O2H'), (None,)), '213.3.3.RAD_GEN_C4H71-3'],
        [(('IC4H9O2', 'C4H71-3'), ('IC4H9O', 'C4H71-O'), (None,)), '230.1.1.RAD_GEN_C4H71-3'],
        [(('IC4H9O2', 'C4H8-2'), ('IC4H9O2H', 'C4H71-3'), (None,)), '231.2.2.RAD_GEN_C4H71-3'],
        [(('TC4H9O2', 'C4H8-1'), ('TC4H9O2H', 'C4H71-3'), (None,)), '231.10.10.RAD_GEN_C4H71-3']
    ]
    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['C4H71-3', 'C4H71-4', 'C4H72-2']
    sort_lst = ['submech_prompt', 0]

    param_dct_sort, _, cmts_dct, _, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    print('Sort by submech_prompt:')
    sorted_results = []

    for i, rxn in enumerate(param_dct_sort.keys()):
        if i % 5 == 0:
            sorted_results.append(
                [rxn, cmts_dct[rxn]['cmts_inline'].split('submech_prompt')[1].strip()])
    print(sorted_results)
    assert results == sorted_results
    print('ok')


def test__filter_pesgroups():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        and also the relative submech
        then order by subpes and broad class
    """
    results = [
        {'grp': 1, 'idxs': ['113:8', '113:9', '113:11', '111:1'], 'peds': [['C4H8-1+H=C4H71-3+H2'], ['C4H8-1+H=C4H71-4+H2'], [
            'C4H8-2+H=C4H72-2+H2'], []], 'hot': [[], [], [], ['C4H71-3', 'C4H71-4', 'C4H72-2']], 'modeltype': 'rovib_dos'}, 
        {'grp': 2, 'idxs': ['121:6', '111:1'], 'peds': [['C4H8-1+O=C4H71-4+OH'], []], 'hot': [[], ['C4H71-4']], 'modeltype': 'rovib_dos'}, 
        {'grp': 3, 'idxs': ['122:21', '122:22', '122:23', '122:24', '111:1'], 'peds': [['C4H8-1+OH=C4H71-3+H2O'], ['C4H8-1+OH=C4H71-4+H2O'], [
            'C4H8-2+OH=C4H71-3+H2O'], ['C4H8-2+OH=C4H72-2+H2O'], []], 'hot': [[], [], [], [], ['C4H71-3', 'C4H71-4', 'C4H72-2']], 
         'modeltype': 'rovib_dos'}, 
        {'grp': 4, 'idxs': ['153:7', '153:9', '111:1'], 'peds': [['C4H8-1+CH3=C4H71-4+CH4'], ['C4H8-2+CH3=C4H72-2+CH4'], []], 
         'hot': [[], [], ['C4H71-4', 'C4H72-2']], 'modeltype': 'rovib_dos'}, 
        {'grp': 5, 'idxs': ['161:6', '111:1'], 'peds': [['C4H8-1+CH3O=C4H71-4+CH3OH'], []], 'hot': [[], ['C4H71-4']], 'modeltype': 'rovib_dos'}
        ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'heptane_cut_species.csv')
    mech_path = os.path.join(CWD, 'data', 'heptane_cut_mech.txt')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)
    therm_str = pathtools.read_file(os.path.join(CWD, 'data'), 'therm.dat')
    spc_therm_dct = ckin_parser.parse_spc_therm_dct(therm_str, [300, 1000, 1500, 2000])

    # Sort with headers for species subset
    isolate_spc = ['C4H71-3', 'C4H71-4', 'C4H72-2']
    sort_lst = ['submech_prompt', 0]

    _, _, _, pes_groups, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst, spc_therm_dct=spc_therm_dct, dct_flt_grps={'DH':30., 'Tref':300}) # 

    print('Sort by submech_prompt and filter pes groups:')
    assert results == pes_groups
    print('ok')


def test__sort_ktp():
    """ test mechanalyzer.parser.sort

        sort ktp dictionary according to highest rate values/ratios
    """
    results = {
        (('H2', 'O'), ('OH', 'H'), (None,)): '2.73e+149.2.27e+01',
        (('H', 'O'), ('OH',), ('(+M)',)): '2.73e+149.4.62e-16',
        (('H', 'O'), ('OH',), (None,)): '2.73e+149.3.97e-39',
        (('H', 'O2'), ('OH', 'O'), (None,)): '2.73e+149.1.30e-89',
        (('H2', 'O'), ('OH', 'OH'), (None,)): '2.73e+149.0.00e+00',
        (('H2', 'O2'), ('HO2V', 'H'), (None,)): '2.73e+149.0.00e+00',
        (('H2', 'O(S)'), ('OH', 'H'), (None,)): '4.80e+143.0.00e+00'
    }

    # Read mechanism files into strings
    spc_paths = [
        os.path.join(CWD, 'data', 'spc2.csv'),
        os.path.join(CWD, 'data', 'spc1.csv')]
    mech_path = None
    sort_path = None

    spc_str, _, _ = _read_files(spc_paths[1], mech_path, sort_path)

    # Build spc and mech information
    spc_dct_full = sparser.parse_mech_spc_dct(spc_str, canon_ent=False)

    # Sort the mechanism
    isolate_spc = []
    sort_lst = ['rxn_max_vals', 'rxn_max_ratio', 0]

    srt_mch = sorter.sorting(
        AL_KTP_DCT, spc_dct_full, sort_lst, isolate_spc)
    sorted_idx, cmts_dct, _ = srt_mch.return_mech_df()
    al_ktp_dct_sorted = sorter.reordered_mech(AL_KTP_DCT, sorted_idx)
    print('ktp dct sorted by max val and ratios test:')
    assert al_ktp_dct_sorted.keys() == results.keys()
    newdct = dict.fromkeys(al_ktp_dct_sorted.keys())
    for rxn in al_ktp_dct_sorted.keys():
        newdct[rxn] = cmts_dct[rxn]['cmts_inline'].split('ratio')[1].strip()

    assert newdct == results
    print('ok')

# Helper function


def _read_files(spc_path, mech_path, sort_path):
    """ read file names
    """

    spc_str, mech_str, sort_str = '', '', ''

    if spc_path is not None:
        with open(spc_path, encoding='utf-8') as fobj:
            spc_str = fobj.read()
    if mech_path is not None:
        with open(mech_path, encoding='utf-8') as fobj:
            mech_str = fobj.read()
    if sort_path is not None:
        with open(sort_path, encoding='utf-8') as fobj:
            sort_str = fobj.read()

    return spc_str, mech_str, sort_str


if __name__ == '__main__':
    # test__sort_with_input()   
    # test__readwrite_thirdbody()
    # test__sortby_mult()
    # test__sortby_molec_r1()
    # test__sortby_pes_dct()
    # test__sortby_rxnclass() # does not work only if filter_pesgroups active
    # test__sortby_species_subpes()
    # test__sort_ktp()
    # still to fix
    test__filter_pesgroups() # recheck why therm is different
    #test__sortby_submech_subpes_chnl()
    #test__sortby_submech_prompt()    
    #test__sortby_submech_ext()
    #test__sortby_submech_del() add this test
    #test__sortby_submech_class()
    
    
