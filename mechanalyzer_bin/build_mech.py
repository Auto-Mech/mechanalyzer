""" Script to generate files with reactions
"""

import sys
import os
import ioformat
import chemkin_io
import mechanalyzer
from mechanalyzer.builder import sorter


# Set up the paths
CWD = os.getcwd()

# Read the build input file and set up info
BLD_STR = ioformat.pathtools.read_file(
    CWD, 'build.dat',
    remove_comments='#', remove_whitespace=True)

FILE_DCT, RSERIES = mechanalyzer.parser.build_input_file(BLD_STR)

# Read the spc_dct and rxn_param_dct
INP_SPC_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_spc'],
    remove_comments='!', remove_whitespace=True)
INP_MECH_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_mech'],
    remove_comments='#', remove_whitespace=True)
SORT_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['sort'],
    remove_comments='#', remove_whitespace=True)

# Check if the input strings exist
if any(string is None for string in (INP_SPC_STR, INP_MECH_STR, SORT_STR)):
    print('ERROR: Input file(s) species.csv, mechanism.dat, sort.dat missing')
    sys.exit()

mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(
    INP_SPC_STR, 'csv')
rxn_param_dct, _, _ = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin', mech_spc_dct)
isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(SORT_STR)

# Generate the requested reactions
mech_spc_dct, rxn_param_dct = mechanalyzer.builder.rxn.build_mechanism(
    mech_spc_dct, rxn_param_dct, rxn_series=RSERIES)

# Write the dictionaries to original strings
csv_str = mechanalyzer.parser.spc.csv_string(
    mech_spc_dct, FILE_DCT['headers'])
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=None,
    mech_spc_dct=mech_spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_param_dct,
    rxn_cmts_dct=None)

# Use strings to generate ordered objects
param_dct_sort, _, mech_spc_dct, cmts_dct, elems = sorter.sorted_mech(
    csv_str, mech_str, isolate_spc, sort_lst)
rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
    rxn_sort_dct=cmts_dct)

# Write the dictionaries to ordered strings
csv_str = mechanalyzer.parser.spc.csv_string(
    mech_spc_dct, FILE_DCT['headers'])
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elems,
    mech_spc_dct=mech_spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_param_dct,
    rxn_cmts_dct=rxn_cmts_dct)

# Write the species and mechanism files
ioformat.pathtools.write_file(csv_str, CWD, FILE_DCT['out_spc'])
ioformat.pathtools.write_file(mech_str, CWD, FILE_DCT['out_mech'])
