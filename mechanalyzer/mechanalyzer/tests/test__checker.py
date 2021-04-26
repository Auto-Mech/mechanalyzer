""" Test the checker.py functions
"""

from mechanalyzer.builder import checker
import numpy as np

TEMPS = np.array([500,1000,1500])
GOOD_KTS = ([1e9, 1e10, 1e10])
NEGATIVE_KTS = np.array([-1.5, -1.8, -1.9])
LARGE_UNIMOLEC_KTS = np.array([1e12, 1e12, 1e12])
LARGE_BIMOLEC_KTS = np.array([1e16, 1e16, 1e16])
LARGE_TERMOLEC_KTS = np.array([1e23, 1e23, 1e23])
ARR_PARAMS = [0, None, None, None, None, None]
LIND_PARAMS = [0, 0, None, None, None, None]
TROE_PARAMS = [0, 0, 0, None, None, None]
CHEB_PARAMS = [0, None, None, 0, None, None]
PLOG_PARAMS = [None, None, None, None, 0, None]

# For testing sources and sinks and duplicates
RXN_PARAM_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): (ARR_PARAMS, ARR_PARAMS, ARR_PARAMS),
    (('H', 'O2'), ('OH', 'O'), (None,)): (ARR_PARAMS, PLOG_PARAMS),
    (('H2', 'O'), ('OH', 'OH'), (None,)): (LIND_PARAMS, LIND_PARAMS),
    (('H', 'O'), ('OH',), (None,)): (TROE_PARAMS, TROE_PARAMS),
    (('H', 'O'), ('OH',), ('(+M)',)): (CHEB_PARAMS, CHEB_PARAMS),
    (('H', 'O'), ('OH',), ('+O(S)',)): (PLOG_PARAMS, PLOG_PARAMS),
    (('H2', 'O(S)'), ('OH', 'O'), (None,)): (PLOG_PARAMS, PLOG_PARAMS),
    (('H2', 'O2'), ('HO2', 'H'), (None,)): (PLOG_PARAMS, PLOG_PARAMS)
}

# For testing sources and sinks and duplicates
RXN_PARAM_DCT2 = {
    (('H', 'O'), ('OH',), (None,)): (PLOG_PARAMS, PLOG_PARAMS),
    (('OH',), ('H', 'O'), (None,)): (PLOG_PARAMS, PLOG_PARAMS)
}

# For testing negative rates
RXN_KTP_DCT1 = {
    (('H2', 'O'), ('OH', 'H'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
    (('H', 'O2'), ('OH', 'O'), (None,)): {1: (TEMPS, NEGATIVE_KTS), 10: (TEMPS, GOOD_KTS)},
}

# For testing large rates
RXN_KTP_DCT2 = {
    (('OH',), ('H', 'O'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_UNIMOLEC_KTS)},
    (('HO2',), ('HO', 'O'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
    (('H2', 'O'), ('OH', 'H'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_BIMOLEC_KTS)},
    (('H2', 'O'), ('OH', 'OH'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
    (('H', 'O'), ('OH',), ('+O(S)',)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_TERMOLEC_KTS)},
    (('H', 'O'), ('OH',), ('+M',)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
}

# For testing large rates
RXN_KTP_DCT3 = {
    (('OH',), ('H', 'O'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, LARGE_UNIMOLEC_KTS)},
    (('HO2',), ('HO', 'O'), (None,)): {1: (TEMPS, GOOD_KTS), 10: (TEMPS, GOOD_KTS)},
}

CORRECT_NEGATIVE_KTS_STR = '\nNEGATIVE RATE CONSTANTS\n\nH+O2=OH+O\nPressure: 1 atm\n' + \
              '    Temperature (K)\n    500.0       1000.0      1500.0      \n' + \
              '    Rate constant\n    -1.500E+00  -1.800E+00  -1.900E+00  \n\n\n'

CORRECT_LARGE_KTS_STR = '\nLARGE RATE CONSTANTS\n\nUnimolecular threshold: 1.0E+11 s^-1\n' + \
    'Bimolecular threshold: 1.0E+15 cm^3 mol^-1 s^-1\n' + \
    'Termolecular threshold: 1.0E+22 cm^6 mol^-2 s^-1\n\n' + \
    'Unimolecular rate constants that exceed 1.0E+11 s^-1\n\n' + \
    'OH=H+O\nPressure: 10 atm\n' + \
    '    Temperature (K)\n    500.0       1000.0      1500.0      \n' + \
    '    Rate constant\n    1.000E+12   1.000E+12   1.000E+12   \n\n\n' + \
    'No bimolecular reactions exceed 1.0E+15 cm^3 mol^-1 s^-1\n\n' + \
    'No termolecular reactions exceed 1.0E+22 cm^6 mol^-2 s^-1\n\n'

CORRECT_LONE_SPCS_STR = '\nLONE SPECIES\n\nThese species appear in 2 or less reactions\n\n' +\
    'Species  Reactions\nO2       H+O2=OH+O, H2+O2=HO2+H\nO(S)     H2+O(S)=OH+O\nHO2      ' +\
    'H2+O2=HO2+H\n\n\n'

CORRECT_SOURCE_SINK_STR1 = '\nSOURCE AND SINK SPECIES\n\nThese species only appear as ' +\
    'reactants:\nSpecies    Reactions\nH2         H2+O=OH+H, H2+O=OH+OH, H2+O(S)=OH+O, ' +\
    'H2+O2=HO2+H\nO(S)       H2+O(S)=OH+O\nO2         H+O2=OH+O, H2+O2=HO2+H\n\n' +\
    'These species only appear as products:\nSpecies   Reactions\nHO2       H2+O2=HO2+H' +\
    '\nOH        H2+O=OH+H, H+O2=OH+O, H2+O=OH+OH, H+O=OH, H+O(+M)=OH(+M), ' +\
    'H+O+O(S)=OH+O(S), H2+O(S)=OH+O\n\n\n'

CORRECT_SOURCE_SINK_STR2 = '\nSOURCE AND SINK SPECIES\n\nThese species only appear as ' +\
    'reactants:\nNo source species found\n\nThese species only appear as products:\n' +\
    'No sink species found\n\n\n'

CORRECT_DUPLICATES_STR1 = '\nDUPLICATE REACTIONS\n\nThese reactions have more than 2 ' +\
    'rate expressions:\n(Number of rate expressions given in parentheses)\n\n' +\
    'H2+O=OH+H     (3)\n\n\n'

CORRECT_DUPLICATES_STR2 = '\nDUPLICATE REACTIONS\n\nThese reactions have more than 2 ' +\
    'rate expressions:\n(Number of rate expressions given in parentheses)\n\n' +\
    'No reactions with more than 2 expressions found\n\n\n'

CORRECT_MISMATCHES_STR1 = '\nMISMATCHED REACTIONS\n\nThe following reactions have ' +\
    'mismatched rate expressions\nH+O2=OH+O: Arrhenius, PLOG\n\n\n'

CORRECT_MISMATCHES_STR2 = '\nMISMATCHED REACTIONS\n\nNo reactions with mismatching ' +\
    'rate expressions found\n\n\n'


def test__all_checks():
    """ Test the run_all_checks function
    """
    k_thresholds = [1e11, 1e15, 1e22]
    rxn_num_threshold = 2
    _ = checker.run_all_checks(RXN_PARAM_DCT1, RXN_KTP_DCT1, k_thresholds,
                                       rxn_num_threshold,
                                       filename='mech_check.txt')


def test__sources_and_sinks():
    """ Test the get_sources_and_sinks and write_sources_and_sinks functions
    """
    # Test the get_sources_and_sinks function for two different cases
    sources1, sinks1 = checker.get_sources_and_sinks(RXN_PARAM_DCT1)
    sources2, sinks2 = checker.get_sources_and_sinks(RXN_PARAM_DCT2)
    assert list(set(list(sources1.keys())) - set(['O2', 'O(S)', 'H2'])) == []
    assert list(set(list(sinks1.keys())) - set(['OH', 'HO2'])) == []
    assert sources2 == sinks2 == {}

    # Test the write_sources_and_sinks function
    source_sink_str1 = checker.write_sources_and_sinks(sources1, sinks1)
    source_sink_str2 = checker.write_sources_and_sinks(sources2, sinks2)
    assert source_sink_str1 == CORRECT_SOURCE_SINK_STR1
    assert source_sink_str2 == CORRECT_SOURCE_SINK_STR2


def test__negative_rates():
    """ Test the get_negative_kts and write_negative_kts functions
    """
    # Test the get_negative_kts function
    negative_rxn_ktp_dct = checker.get_negative_kts(RXN_KTP_DCT1)
    assert tuple(negative_rxn_ktp_dct.keys()) == ((('H', 'O2'), ('OH', 'O'), (None,)),)
    assert tuple(negative_rxn_ktp_dct[(('H', 'O2'), ('OH', 'O'), (None,))].keys()) == (1,)

    # Test the write_negative_kts function
    negative_kts_str = checker.write_negative_kts(negative_rxn_ktp_dct)
    assert negative_kts_str == CORRECT_NEGATIVE_KTS_STR


def test__large_rates():
    """ Test the get_large_kts and write_large_kts functions
    """
    # Test the get_large_kts function
    thresholds = [1e11, 1e15, 1e22]
    large_rxn_ktp_dcts2 = checker.get_large_kts(RXN_KTP_DCT2, thresholds)
    assert tuple(large_rxn_ktp_dcts2[0].keys()) == ((('OH',), ('H', 'O'), (None,)),)
    assert tuple(large_rxn_ktp_dcts2[1].keys()) == ((('H2', 'O'), ('OH', 'H'), (None,)),)
    assert tuple(large_rxn_ktp_dcts2[2].keys()) == ((('H', 'O'), ('OH',), ('+O(S)',)),)

    # Test the write_large_kts function
    large_rxn_ktp_dcts3 = checker.get_large_kts(RXN_KTP_DCT3, thresholds)
    large_kts_str3 = checker.write_large_kts(large_rxn_ktp_dcts3, thresholds)
    assert large_kts_str3 == CORRECT_LARGE_KTS_STR


def test__lone_species():
    """ Test the get_lone_spcs and write_lone_spcs functions
    """
    # Test the get_lone_spcs function
    threshold = 2
    lone_spcs = checker.get_lone_spcs(RXN_PARAM_DCT1, threshold)
    assert tuple(lone_spcs.keys()) == tuple(['O2', 'O(S)', 'HO2'])

    # Test the write_lone_spcs function
    lone_spcs_str = checker.write_lone_spcs(lone_spcs, threshold)
    assert lone_spcs_str == CORRECT_LONE_SPCS_STR


def test__duplicates():
    """ Test the get_duplicates and write_duplicates functions
    """
    # Test the get_duplicates function
    duplicate_rxns1 = checker.get_duplicates(RXN_PARAM_DCT1)
    duplicate_rxns2 = checker.get_duplicates(RXN_PARAM_DCT2)
    assert tuple(duplicate_rxns1.keys()) == ((('H2', 'O'), ('OH', 'H'), (None,)),)
    assert tuple(duplicate_rxns1.values()) == (3,)
    assert duplicate_rxns2 == {}

    # Test the write_duplicates function
    dup_str1 = checker.write_duplicates(duplicate_rxns1)
    dup_str2 = checker.write_duplicates(duplicate_rxns2)
    assert dup_str1 == CORRECT_DUPLICATES_STR1
    assert dup_str2 == CORRECT_DUPLICATES_STR2


def test__mismatches():
    """ Test the get_mismatches and write_mismatches functions
    """
    # Test the get_mismatches function
    mismatched_rxns1 = checker.get_mismatches(RXN_PARAM_DCT1)
    mismatched_rxns2 = checker.get_mismatches(RXN_PARAM_DCT2)
    assert tuple(mismatched_rxns1.keys()) == ((('H', 'O2'), ('OH', 'O'), (None,)),)
    assert tuple(mismatched_rxns1.values())[0][1] == ['Arrhenius', 'PLOG']
    assert mismatched_rxns2 == {}

    # Test the write_mismatches function
    mismatch_str1 = checker.write_mismatches(mismatched_rxns1)
    mismatch_str2 = checker.write_mismatches(mismatched_rxns2)
    assert mismatch_str1 == CORRECT_MISMATCHES_STR1
    assert mismatch_str2 == CORRECT_MISMATCHES_STR2


if __name__ == '__main__':
    test__all_checks()
    test__sources_and_sinks()
    test__negative_rates()
    test__large_rates()
    test__lone_species()
    test__duplicates()
    test__mismatches()
