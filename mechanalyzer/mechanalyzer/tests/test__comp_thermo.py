"""
compare thermo
"""

# import os
# from ioformat import remove_whitespace
# from chemkin_io import parser
# from chemkin_io.calculator import combine


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
# HEPTANE_MECH_NAME = 'heptane_mechanism.txt'
# SYNGAS_MECH_NAME = 'syngas_mechanism.txt'
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
#
# # Temperatures to run
# TEMPS = [500.0, 1000.0, 2000.0]
#
#
# def test__compare_thermo():
#     """ test chemkin_io.calculator.combine.mechanism_thermo
#     """
#
#     # Build dictionaries containing the thermo data strings for each species
#     # Dictionaries indexed by the given mechanism names or InCHI string
#     # Build name for test but use inchi
#     mech1_thermo_dct, mech2_thermo_dct = combine.build_thermo_name_dcts(
#         FAKE1_MECH_STR, FAKE2_MECH_STR, TEMPS)
#
#     mech1_thermo_dct, mech2_thermo_dct = combine.build_thermo_inchi_dcts(
#         FAKE1_MECH_STR, FAKE2_MECH_STR,
#         FAKE_CSV_STR, FAKE_CSV_STR, TEMPS)
#
#     print('\nMech1 thermo dct')
#     for key, val in mech1_thermo_dct.items():
#         print(key)
#         print(val)
#     print('\n\nMech2 thermo dct')
#     for key, val in mech2_thermo_dct.items():
#         print(key)
#         print(val)
#
#     # Calculate Enthalpy, Entropy, Gibbs, & Heat Capacity for each species
#     # Mech 1 and 2 thermo data collated for each species
#     thermo_vals_dct = combine.mechanism_thermo(
#         mech1_thermo_dct, mech2_thermo_dct)
#
#     print('\n\n\nCombined thermo vals')
#     for idx in thermo_vals_dct:
#
#         # Print all of the thermo quantities
#         print(idx, '\n')
#         mech1_vals = thermo_vals_dct[idx]['mech1']
#         mech2_vals = thermo_vals_dct[idx]['mech2']
#         print('\nM1 Enthalpy', mech1_vals[0])
#         print('M2 Enthalpy', mech2_vals[0])
#         print('M1 Heat Capacity', mech1_vals[1])
#         print('M2 Heat Capacity', mech2_vals[1])
#         print('M1 Entropy', mech1_vals[2])
#         print('M2 Entropy', mech2_vals[2])
#         print('M1 Gibbs', mech1_vals[3])
#         print('M2 Gibbs', mech2_vals[3])
#
#     print(thermo_vals_dct)
#
#
# if __name__ == '__main__':
#     test__compare_thermo()
