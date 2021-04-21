"""
compare thermo
"""

# import os
# import numpy
# from ioformat import remove_whitespace
# from chemkin_io import parser
# from chemkin_io.calculator import combine
#
#
# # Get mechanism information
# def _read_file(file_name):
#     with open(file_name, encoding='utf8', errors='ignore') as file_obj:
#         file_str = file_obj.read()
#     return file_str
#
#
# # Set paths
# PATH = os.path.dirname(os.path.realpath(__file__))
# DATA_PATH = os.path.join(PATH, 'data')
# FAKE1_MECH_NAME = 'fake1_mech.txt'
# FAKE2_MECH_NAME = 'fake2_mech.txt'
# FAKE_CSV_NAME = 'fake_species.csv'
#
# # Read mechanism and csv strings
# FAKE1_MECH_STR = _read_file(
#     os.path.join(DATA_PATH, FAKE1_MECH_NAME))
# FAKE2_MECH_STR = _read_file(
#     os.path.join(DATA_PATH, FAKE2_MECH_NAME))
# FAKE_CSV_STR = _read_file(
#     os.path.join(DATA_PATH, FAKE_CSV_NAME))
#
# # Read species blocks
# FAKE1_THERMO_BLOCK = parser.mechanism.thermo_block(
#     FAKE1_MECH_STR)
# FAKE1_BLOCK_STRS = parser.thermo.data_strings(
#     FAKE1_THERMO_BLOCK)
# FAKE2_THERMO_BLOCK = remove_whitespace(
#     parser.mechanism.thermo_block(FAKE2_MECH_STR))
# FAKE2_BLOCK_STRS = parser.thermo.data_strings(
#     FAKE2_THERMO_BLOCK)
#
# # Temperatures and Pressures to run
# T_REF = 1.0
# TEMPS = numpy.array([500.0, 1000.0, 1500.0])
# PRESSURES = numpy.array([1.0, 4.0, 5.0])
#
#
# def test__compare_rates():
#     """ test chemkin_io.mechparser.compare.rates
#     """
#
#     # Build dictionaries containing:
#     # thermo data strings for each species, k data strings for each reaction
#     # Dictionaries indexed by the given mechanism names or InCHI string
#     # Build name for test but use inchi
#
#     _, mech2_thermo_dct = combine.build_thermo_name_dcts(
#         FAKE1_MECH_STR, FAKE2_MECH_STR, TEMPS)
#     mech1_ktp_dct, mech2_ktp_dct = combine.build_reaction_name_dcts(
#         FAKE1_MECH_STR, FAKE2_MECH_STR,
#         T_REF, TEMPS, PRESSURES)
#
#     _, mech2_thermo_dct = combine.build_thermo_inchi_dcts(
#         FAKE1_MECH_STR, FAKE2_MECH_STR, FAKE_CSV_STR, FAKE_CSV_STR, TEMPS)
#     mech1_ktp_dct, mech2_ktp_dct = combine.build_reaction_inchi_dcts(
#         FAKE1_MECH_STR, FAKE2_MECH_STR, FAKE_CSV_STR, FAKE_CSV_STR,
#         T_REF, TEMPS, PRESSURES)
#
#     print('\nMech1 reaction dct')
#     for key, val in mech1_ktp_dct.items():
#         print(key)
#         print(val)
#     print('\n\nMech2 reaction dct')
#     for key, val in mech2_ktp_dct.items():
#         print(key)
#         print(val)
#
#     # Rate constant
#     ktp_dct = combine.mechanism_rates(
#         mech1_ktp_dct, mech2_ktp_dct,
#         mech2_thermo_dct,
#         TEMPS)
#
#     # print dict
#     print('\n\nktp dct')
#     for rxn, mechs in ktp_dct.items():
#         print(rxn)
#         mech1, mech2 = mechs['mech1'], mechs['mech2']
#         for (pr1, ktp1), (pr2, ktp2) in zip(mech1.items(), mech2.items()):
#             print('Pressure: ', pr1, pr2)
#             print('Mech1 ks: ', ktp1)
#             print('Mech2 ks: ', ktp2)
#             print(' ')
#
#     print('\n\nktp dct')
#     print(ktp_dct)
#
#
# if __name__ == '__main__':
#     test__compare_rates()
