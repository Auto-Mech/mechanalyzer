""" Script to generate files with reactions
"""

import os
import time
import ioformat
import chemkin_io
import mechanalyzer
from mechanalyzer.builder import sorter


# Initialize the start time for script execution
t0 = time.time()

# Set useful global variables
CWD = os.getcwd()

# READ AND PARSE INPUT
print('\n---- Parsing the species and mechanism files ---\n')

# Read the build input file and set up info
BLD_STR = ioformat.pathtools.read_file(
    CWD, 'build.dat',
    remove_comments='#', remove_whitespace=True)

FILE_DCT, _ = mechanalyzer.parser.build_input_file(BLD_STR)

# Read input species and mechanism files into dictionary
INP_SPC_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_spc'],
    remove_comments='!', remove_whitespace=True)
INP_MECH_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['inp_mech'],
    remove_comments='#', remove_whitespace=True)
SORT_STR = ioformat.pathtools.read_file(
    CWD, FILE_DCT['sort'],
    remove_comments='#', remove_whitespace=True)

# Build the initial dictionaries
spc_dct = mechanalyzer.parser.spc.build_spc_dct(INP_SPC_STR, 'csv')
rxn_param_dct, _, _ = mechanalyzer.parser.mech.parse_mechanism(
    INP_MECH_STR, 'chemkin', spc_dct)
isolate_spc, sort_lst = mechanalyzer.parser.mech.parse_sort(SORT_STR)

# ADD STEREO TO SPC AND REACTIONS
print('\n---- Adding stereochemistry to InChIs of mechanism'
      ' species where needed ---\n')

# Add the stereochemical labels to the species
spc_dct = mechanalyzer.parser.spc.stereochemical_spc_dct(
    spc_dct, nprocs='auto', all_stereo=False)

print('Mechanism species with stereo added')
for name, dct in spc_dct.items():
    print('Name: {0:<25s} InChI: {1}'.format(name, dct['inchi']))

# Expand the reactions in the mechanism to include stereochemical variants
print('\n---- Expanding the list of mechanism reactions to include all'
      ' valid, stereoselective permutations ---\n')
ste_rxn_dct, ste_spc_dct = mechanalyzer.builder.expand_mech_stereo(
    rxn_param_dct, spc_dct)

# OBTAIN SORTED SPECIES AND MECHANISMS

# Write the new (sorted) species dictionary to a string
HEADERS = ('smiles', 'inchi', 'inchikey', 'mult', 'charge')
ste_spc_dct_sort = mechanalyzer.parser.spc.reorder_by_atomcount(ste_spc_dct)
csv_str = mechanalyzer.parser.spc.csv_string(ste_spc_dct_sort, HEADERS)

# Write initial string to call the sorter
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=None,
    spc_dct=ste_spc_dct_sort,
    spc_nasa7_dct=None,
    rxn_param_dct=rxn_param_dct,
    comments=None)

# Use strings to generate ordered objects
param_dct_sort, _, ste_spc_dct_sort, cmts_dct, elems = sorter.sorted_mech(
    csv_str, mech_str, isolate_spc, sort_lst)

# WRITE OUTPUT AND EXIT

# Write the dictionaries to ordered strings
csv_str = mechanalyzer.parser.spc.csv_string(
    ste_spc_dct_sort, FILE_DCT['headers'])
mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elems,
    spc_dct=ste_spc_dct_sort,
    spc_nasa7_dct=None,
    rxn_param_dct=param_dct_sort,
    comments=cmts_dct)

# Write the species and mechanism files
ioformat.pathtools.write_file(csv_str, CWD, FILE_DCT['out_spc'])
ioformat.pathtools.write_file(mech_str, CWD, FILE_DCT['out_mech'])

# Compute script run time and print to screen
tf = time.time()
print('\n\nScript executed successfully.')
print('Time to complete: {:.2f}'.format(tf-t0))
print('Exiting...')