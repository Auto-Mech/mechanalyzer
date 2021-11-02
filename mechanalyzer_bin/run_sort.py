"""
Sort by selected criteria written in this script
"""

import sys
import os
import argparse
from ioformat import pathtools
import chemkin_io.writer
from mechanalyzer.builder import sorter
from mechanalyzer.parser import spc as sparser
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
PAR.add_argument('-o', '--outmech', default='outmech.dat',
                 help='output file name (outmech.dat)')
PAR.add_argument('-c', '--outspc', default='outspc.csv',
                 help='output file name (outspc.csv)')
OPTS = vars(PAR.parse_args())

# Read the input files
spc_str = pathtools.read_file(CWD, OPTS['spc'], remove_comments='#')
mech_str = pathtools.read_file(CWD, OPTS['mech'], remove_comments='#')
sort_str = pathtools.read_file(CWD, OPTS['sort'], remove_comments='#')

# Check if the input strings exist
if any(string is None for string in (spc_str, mech_str, sort_str)):
    print('ERROR: Input file(s) species.csv, mechanism.dat, sort.dat missing')
    sys.exit()

# Build sorted mechanism files
isolate_spc, sort_lst = mparser.parse_sort(sort_str)
param_dct_sort, _, mech_spc_dct, cmts_dct, elems = sorter.sorted_mech(
    spc_str, mech_str, isolate_spc, sort_lst)
rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
    rxn_sort_dct=cmts_dct)
print(rxn_cmts_dct)

# Write the output files (need to make general at some point)
headers = sparser.csv_headers(mech_spc_dct)
sortd_csv_str = sparser.csv_string(mech_spc_dct, headers)
sortd_mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    elem_tuple=elems,
    mech_spc_dct=mech_spc_dct,
    rxn_param_dct=param_dct_sort,
    rxn_cmts_dct=rxn_cmts_dct)

pathtools.write_file(sortd_csv_str, CWD, OPTS['outspc'])
pathtools.write_file(sortd_mech_str, CWD, OPTS['outmech'])
