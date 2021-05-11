"""
Sort by selected criteria written in this script
"""

import sys
import os
import argparse
from ioformat import pathtools
import chemkin_io.writer
from mechanalyzer.builder import sorter
from mechanalyzer.parser import mech as mparser


# Set useful global variables
CWD = os.getcwd()

# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-m', '--mech', default='mechanism.dat',
                 help='mechanism file name (mechanism.dat)')
PAR.add_argument('-s', '--spc', default='species.csv',
                 help='species file name (species.csv)')
PAR.add_argument('-i', '--sort', default='sort.dat',
                 help='input file name (sort.dat)')
PAR.add_argument('-o', '--out', default='outmech.dat',
                 help='output file name (outmech.dat)')
# PAR.add_argument('-o2', '--out2', default='outmech2.dat',
#                  help='output file name (outmech2.dat)')
OPTS = vars(PAR.parse_args())

# Read the input files
spc_str = pathtools.read_file(CWD, OPTS['spc'], remove_comments='#')
mech_str = pathtools.read_file(CWD, OPTS['mech'], remove_comments='#')
sort_str = pathtools.read_file(CWD, OPTS['sort'], remove_comments='#')

# Check if the input strings exist
if any(string is None for string in (spc_str, mech_str, sort_str)):
    print('ERROR: Input file missing')
    sys.exit()

# Build sorted mechanism files
isolate_spc, sort_list = mparser.parse_sort(sort_str)
param_dct_sort, _, spc_dct, cmts_dct, elems = sorter.sorted_mech(
    spc_str, mech_str, isolate_spc, sort_str)

# Write the output files (need to make general at some point)
sortd_mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elems, spc_dct=spc_dct,
    rxn_param_dct=param_dct_sort,
    comments=cmts_dct)

sortmech_out_path = os.path.join(CWD, OPTS['out'])
# restmech_out_path = os.path.join(CWD, OPTS['out2'])
with open(sortmech_out_path, 'w') as mech_obj:
    mech_obj.write(sortd_mech_str)
