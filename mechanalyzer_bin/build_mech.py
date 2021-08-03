""" Script to generate files with reactions
"""

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

spc_dct = mechanalyzer.parser.spc.build_spc_dct(
    INP_SPC_STR, 'csv')
rxn_param_dct, _, _ = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin', spc_dct)
isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(SORT_STR)

# Generate the requested reactions
spc_dct, rxn_param_dct = mechanalyzer.builder.rxn.build_mechanism(
    spc_dct, rxn_param_dct, rxn_series=RSERIES)

# Write the dictionaries to original strings
csv_str = mechanalyzer.parser.spc.csv_str(
    spc_dct, FILE_DCT['headers'])
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=None,
    spc_dct=spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_param_dct,
    comments=None)

# Use strings to generate ordered objects
param_dct_sort, _, spc_dct, cmts_dct, elems = sorter.sorted_mech(
    csv_str, mech_str, isolate_spc, sort_lst)

# Write the dictionaries to ordered strings
csv_str = mechanalyzer.parser.spc.csv_str(
    spc_dct, FILE_DCT['headers'])
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elems,
    spc_dct=spc_dct,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_param_dct,
    comments=cmts_dct)

# Write the species and mechanism files
ioformat.pathtools.read_file(CWD, FILE_DCT['inp_spc'])
ioformat.pathtools.read_file(CWD, FILE_DCT['inp_mech'])
