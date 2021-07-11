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
         '! class: N_COH.subpes _NR _rxn type  2004.0 _1 _Decomposition'],
        [(('C2H3', 'H'), ('C2H4',), ('(+M)',)),
         '! class: N_COH.subpes _NR _rxn type  2004.0 _2 _Recombination H'],
        [(('C2H4', 'H'), ('C2H5',), ('(+M)',)),
         '! class: N_COH.subpes _NR _rxn type  2005.0 _2 _Addition H'],
        [(('C2H4', 'H'), ('C2H3', 'H2'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2005.01 _2 _H abstraction'],
        [(('CH2(S)', 'CH3'), ('C2H4', 'H'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2005.02 _2 _Addition-decomposition - propagation'],
        [(('C2H5', 'H'), ('C2H4', 'H2'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2006.0 _2 _Recombination-decomposition - termination'],
        [(('C2H4', 'O'), ('CH3', 'HCO'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2104.0 _2 _Addition-decomposition - branching'],
        [(('C2H4', 'OH'), ('PC2H4OH',), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2105.0 _2 _Addition OH'],
        [(('C2H5', 'OH'), ('C2H4', 'H2O'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2106.0 _2 _Recombination-decomposition - termination'],
        [(('C2H4', 'O2'), ('C2H3', 'HO2'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2204.0 _2 _H abstraction'],
        [(('C2H5O2',), ('C2H4', 'HO2'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  2205.0 _1 _Beta-scission +HO2'],
        [(('C2H4', 'CH3'), ('C2H3', 'CH4'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  3007.0 _2 _H abstraction'],
        [(('C3H4-A', 'O'), ('C2H4', 'CO'), (None,)),
         '! class: N_COH.subpes _NR _rxn type  3104.0 _2 _Addition-decomposition - termination']
    ]

    # Read the mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech.dat')
    sort_path = os.path.join(CWD, 'data', 'sort.dat')

    spc_str, mech_str, sort_str = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism
    isolate_spc, sort_lst = mparser.parse_sort(sort_str)

    param_dct_sort, _, _, cmts_dct, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    index = 0
    for rxn in param_dct_sort.keys():
        assert [rxn, cmts_dct[rxn]['cmts_inline']] == results[index]
        index += 1


def test__readwrite_thirdbody():
    """ test mechanalyzer.parser.sort

        Checks read/write of a small set of rxns involving third bodies
    """

    trd_bdy_dct = {
        (('H', 'OH'), ('H2O',), ('+M',)):
        (([3.5e+22, -2.0, 0.0], None, None, None, None,
          {'H2': 0.73, 'H2O': 3.65, 'CH4': 2.0, 'C2H6': 3.0, 'AR': 0.38}),),
        (('H', 'OH', 'AR'), ('H2O', 'AR'), (None,)):
        (([3.5e+22, -2.0, 0.0], None, None, None, None, None),),
        (('H', 'O2'), ('HO2',), ('(+HE)',)):
        (([4650000000000.0, 0.44, 0.0], [9.19e+18, -1.2, 0.0],
          [0.59, 1e-30, 1e+30, 1e+30], None, None, None),),
        (('CH3', 'IC4H7'), ('AC5H10',), ('(+M)',)):
        (([150000000000000.0, -0.32, -262.3],
          [5.86e+60, -12.81, 6250.0],
            [0.104, 1606.0, 60000.0, 6118.0], None, None,
            {'H2': 2.0, 'H2O': 6.0, 'CH4': 2.0, 'CO': 1.5,
             'CO2': 2.0, 'C2H6': 3.0, 'AR': 0.7}),),
        (('C5H10-2',), ('C4H71-3', 'CH3'), ('(+M)',)):
        (([6.486e+19, -1.367, 76320], [1.53e+104, -24.826, 94800.0],
          [0.005301, 143.7, 16770000000000.0, 3671.0], None, None, None),),
        (('C5H11-1',), ('C2H4', 'NC3H7'), (None,)):
        (([8.06e+20, -2.628, 29232], None, None, None,
          {0.1: [4410.0, 2.192, 18827.0],
           1.0: [8.06e+20, -2.628, 29232.0],
           10.0: [2.17e+28, -4.578, 34864.0],
           100.0: [6.47e+24, -3.383, 34388.0],
           1000.0: [2.34e+17, -1.123, 31176.0]}, None),)}

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

    assert param_dct_sort == trd_bdy_dct


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

    param_dct_sort, _, _, cmts_dct, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    for rxn in param_dct_sort.keys():
        assert cmts_dct[rxn]['cmts_inline'][-1] == results[rxn]


def test__sortby_molec_r1():
    """ test mechanalyzer.parser.sort

        Sort by first (heavier) reactant and molecularity of the reaction
    """
    comments_results = [
        'C2H3_2', 'C2H3_2', 'C2H3_2', 'C2H3_2', 'C2H3_2', 'C2H3_2',
        'C2H3OO_1', 'C2H3OO_1', 'C2H4_1', 'C2H4_2', 'C2H4_2', 'C2H4_2',
        'C2H4_2', 'C2H4_2', 'C2H4_2', 'C2H5_2', 'C2H5_2', 'C2H5_2',
        'C2H5O2_1', 'C2H6_2', 'C3H4-A_2', 'CH3_2', 'CH3_2', 'HOCH2CO_1']

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_C2H4_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by R1-molecularity - No Headers Included
    isolate_spc = []
    sort_lst = ['r1', 'molecularity', 0]  # NO HEADERS INCLUDED
    sort_lst_2 = ['molecularity', 0]  # ONLY BY MOLEC- CHECK HEADERS WITH INT NUMBER

    param_dct_sort, _, _, cmts_dct, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    param_dct_sort_2, _, _, cmts_dct_2, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst_2)
        # NO NEED TO CALL IT AFTERWARDS

    comments = []
    for rxn in param_dct_sort.keys():
        comments.append(''.join(cmts_dct[rxn]['cmts_inline'].split()[-2:]))

    assert comments == comments_results


def test_sortby_pes_dct():
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


def test_sortby_rxnclass():
    """ test mechanalyzer.parser.sort

        sort by reaction class:
            both "broad" (based on multiplicity, reactants/products..)
        and "graph" (based on graph classification - warning, CPU intensive)
        prior to rxn class, the mech is also subdivided into PESs
    """
    results = [
        [(('C2H5', 'H'), ('C2H6',), ('(+M)',)),
         '  addition _Recombination H'],
        [(('C2H3', 'H'), ('C2H4',), ('(+M)',)),
         '  addition _Recombination H'],
        [(('HOCH2CO',), ('CH2OH', 'CO'), (None,)),
         '  beta scission _Decomposition'],
        [(('C2H3OO',), ('CH2O', 'HCO'), (None,)),
         '  elimination _Beta-scission'],
        [(('C2H5O2',), ('C2H4', 'HO2'), (None,)),
         '  elimination _Beta-scission +HO2'],
        [(('C2H4',), ('H2', 'H2CC'), ('(+M)',)),
         '  elimination _Decomposition'],
        [(('C2H4', 'H'), ('C2H3', 'H2'), (None,)),
         '  hydrogen abstraction _H abstraction'],
        [(('C2H5', 'H'), ('C2H4', 'H2'), (None,)),
         '  hydrogen abstraction _Recombination-decomposition - termination'],
        [(('C2H3', 'O2'), ('C2H2', 'HO2'), (None,)),
         '  hydrogen abstraction _Recombination-decomposition - termination'],
        [(('C3H5-A',), ('C3H5-T',), (None,)),
         '  hydrogen migration _Isomerization'],
        [(('C3H5-A',), ('C3H5-S',), (None,)),
         '  hydrogen migration _Isomerization'],
        [(('CH2(S)', 'CH3'), ('C2H4', 'H'), (None,)),
         '  substitution _Addition-decomposition - propagation'],
        [(('CH3', 'CH3'), ('H', 'C2H5'), (None,)),
         '  substitution _Recombination-decomposition - propagation'],
        [(('C3H4-A', 'O'), ('C2H4', 'CO'), (None,)),
         '  unclassified _Addition-decomposition - termination']
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

    param_dct_sort, _, _, cmts_dct, _ = sorter.sorted_mech(
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
         '! subclass: N_COH.subpes  8018.00'],
        [(('IC8', 'O2'), ('IC8-1R', 'HO2'), (None,)),
         '! subclass: N_COH.subpes  8218.00'],
        [(('IC8-1R',), ('IC8-5R',), (None,)),
         '! subclass: N_COH.subpes  8017.00'],
        [(('IC8-1R', 'O2'), ('IC8-1O2R',), (None,)),
         '! subclass: N_COH.subpes  8217.00'],
        [(('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'), (None,)),
         '! subclass: N_COH.subpes  9220.00'],
    ]
    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_submech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort mechanism by subpes with Headers for specues subset
    isolate_spc = ['IC8', 'IC8-1R']
    sort_lst = ['species', 'subpes', 1]

    print('Sort by species subset + subpes test:')

    param_dct_sort, _, _, cmts_dct, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline']])

    assert sorted_results == results

    print('ok')


def test__sortby_submech_class():
    """ test mechanalyzer.parser.sort

        sort by fuel submechanism: extract reactions of
        fuel, fuel radicals, R+O2, R+O4
        then order by subpes and broad class
    """
    results = [
        [(('IC8',), ('NEOC5H11', 'IC3H7'), (None,)),
         '  FUEL _8018.0 _Bond fission'],
        [(('IC8', 'O2'), ('IC8-1R', 'HO2'), (None,)),
         '  FUEL _8218.01 _H abstraction'],
        [(('IC8OOH1',), ('IC8-1OR', 'OH'), (None,)),
         '  FUEL_ADD_O2 _8218.0 _Bond fission +OH'],
        [(('IC8-1O2R', 'HO2'), ('IC8OOH1', 'O2'), (None,)),
         '  FUEL_ADD_O2 _8418.0 _Recombination-decomposition - termination'],
        [(('IC8-1O2R', 'H2O2'), ('IC8OOH1', 'HO2'), (None,)),
         '  FUEL_ADD_O2 _8419.0 _H abstraction'],
        [(('IC8-1R',), ('IC8-5R',), (None,)),
         '  FUEL_RAD _8017.0 _Isomerization'],
        [(('IC8-1R', 'O2'), ('IC8-1O2R',), (None,)),
         '  FUEL_RAD _8217.0 _Recombination O2'],
        [(('IC8-3R', 'O2'), ('IC8D3', 'HO2'), (None,)),
         '  FUEL_RAD _8217.02 _Recombination-decomposition - termination'],
        [(('IC8-1R', 'CH3O2'), ('IC8-1OR', 'CH3O'), (None,)),
         '  FUEL_RAD _9220.0 _Recombination-decomposition - propagation'],
        [(('IC8OOH1-1AR',), ('IC4H7OOH', 'IC4H9'), (None,)),
         '  R_O2 _8217.0 _Beta-scission'],
        [(('IC8OOH1-1AR',), ('IC8O1-1A', 'OH'), (None,)),
         '  R_O2 _8217.0 _Beta-scission +OH'],
        [(('IC8OOH1-1AR',), ('CH2O', 'I24C7D1', 'OH'), (None,)),
         '  R_O2 _8217.0 _Decomposition(lumped)'],
        [(('IC8-1O2R',), ('IC8OOH1-1AR',), (None,)),
         '  R_O2 _8217.0 _Isomerization'],
        [(('IC8-3O2R',), ('IC8D3', 'HO2'), (None,)),
         '  R_O2 _8217.01 _Beta-scission +HO2'],
        [(('IC8OOH1-1AR', 'O2'), ('IC8OOH1-1AO2R',), (None,)),
         '  R_O2 _8417.0 _Recombination O2']
    ]

    # Read mechanism files into strings
    spc_path = os.path.join(CWD, 'data', 'LLNL_species.csv')
    mech_path = os.path.join(CWD, 'data', 'LLNL_IC8_red_mech.dat')
    sort_path = None

    spc_str, mech_str, _ = _read_files(spc_path, mech_path, sort_path)

    # Sort with headers for species subset
    isolate_spc = ['IC8']
    sort_lst = ['submech', 'subpes', 'rxn_class_broad', 0]

    param_dct_sort, _, _, cmts_dct, _ = sorter.sorted_mech(
        spc_str, mech_str, isolate_spc, sort_lst)

    print('Sort by submech-subpes-classbroad test:')
    sorted_results = []
    for rxn in param_dct_sort.keys():
        sorted_results.append(
            [rxn, cmts_dct[rxn]['cmts_inline'].split('type')[1]])

    assert results == sorted_results
    print('ok')


def test_sort_ktp():
    """ test mechanalyzer.parser.sort

        sort ktp dictionary according to highest rate values/ratios
    """
    results = {
        (('H2', 'O'), ('OH', 'H'), (None,)): '2.72670689e+149 _22.692539074854153',
        (('H', 'O'), ('OH',), ('(+M)',)): '2.72670689e+149 _4.619780379447244e-16',
        (('H', 'O'), ('OH',), (None,)): '2.72670689e+149 _3.973593578206074e-39',
        (('H', 'O2'), ('OH', 'O'), (None,)): '2.72670689e+149 _1.3017843126207617e-89',
        (('H2', 'O'), ('OH', 'OH'), (None,)): '2.72670689e+149 _0.0',
        (('H2', 'O2'), ('HO2V', 'H'), (None,)): '2.72670689e+149 _0.0',
        (('H2', 'O(S)'), ('OH', 'H'), (None,)): '4.79926202e+143 _0.0'
    }

    # Read mechanism files into strings
    spc_paths = [
        os.path.join(CWD, 'data', 'spc2.csv'),
        os.path.join(CWD, 'data', 'spc1.csv')]
    mech_path = None
    sort_path = None

    spc_str, _, _ = _read_files(spc_paths[1], mech_path, sort_path)

    # Build spc and mech information
    spc_dct_full = sparser.build_spc_dct(spc_str, SPC_TYPE)
    mech_info = mparser.mech_info(AL_KTP_DCT, spc_dct_full)

    # Sort the mechanism
    isolate_spc = []
    sort_lst = ['rxn_max_vals', 'rxn_max_ratio', 0]

    srt_mch = sorter.sorting(
        mech_info, spc_dct_full, sort_lst, isolate_spc)
    sorted_idx, cmts_dct, _ = srt_mch.return_mech_df()
    al_ktp_dct_sorted = sorter.reordered_mech(AL_KTP_DCT, sorted_idx)
    print('ktp dct sorted by max val and ratios test:')
    assert al_ktp_dct_sorted.keys() == results.keys()
    newdct = dict.fromkeys(al_ktp_dct_sorted.keys())
    for rxn in al_ktp_dct_sorted.keys():
        newdct[rxn] = cmts_dct[rxn]['cmts_inline'].split('ratio')[1].strip()

    assert newdct == results


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


if __name__ == '__main__':
   # __sort_with_input()
   # test__readwrite_thirdbody()
   # test__sortby_mult()
   # test__sortby_molec_r1()
   # test_sortby_pes_dct()
   # test_sortby_rxnclass()
   # test__sortby_species_subpes()
   # test__sortby_submech_class()
   # test_sort_ktp()
