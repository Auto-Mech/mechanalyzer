"""
Sort by selected criteria written in this script
"""

import sys
import os
import argparse
import numpy as np
from ioformat import pathtools
import chemkin_io.writer
from mechanalyzer.builder import sorter
from mechanalyzer.builder import sort_fct
from mechanalyzer.parser import spc as sparser
from mechanalyzer.parser import mech as mparser
from mechanalyzer.parser import ckin_ as ckin_parser

# Set useful global variables
CWD = os.getcwd()

# Parse the command line
PAR = argparse.ArgumentParser()
PAR.add_argument('-m', '--mech', default='mechanism.dat',
                 help='mechanism file name (mechanism.dat)')
PAR.add_argument('-s', '--spc', default='species.csv',
                 help='species file name (species.csv)')
PAR.add_argument('-t', '--therm', default='therm.dat',
                 help='thermo file name (therm.dat)')
PAR.add_argument('-i', '--sort', default='sort.dat',
                 help='input file name (sort.dat)')
PAR.add_argument('-o', '--outmech', default='outmech.dat',
                 help='output file name (outmech.dat)')
PAR.add_argument('-c', '--outspc', default='outspc.csv',
                 help='output file name (outspc.csv)')
OPTS = vars(PAR.parse_args())

# Read the input files
spc_str = pathtools.read_file(CWD, OPTS['spc'], remove_comments='!')
mech_str = pathtools.read_file(CWD, OPTS['mech'], remove_comments='!')
therm_str = pathtools.read_file(CWD, OPTS['therm'], remove_comments='!')
sort_str = pathtools.read_file(CWD, OPTS['sort'], remove_comments='#')

# Check if the input strings exist
if any(string is None for string in (spc_str, mech_str, sort_str)):
    print('ERROR: Input file(s) species.csv, mechanism.dat, sort.dat missing')
    sys.exit()

# read thermo and filter pes groups
spc_therm_dct = ckin_parser.parse_spc_therm_dct(therm_str, np.arange(300,2010,10))

# Build sorted mechanism files
isolate_spc, sort_lst = mparser.parse_sort(sort_str)
param_dct_sort, mech_spc_dct, cmts_dct, pes_groups, rxns_filter = sorter.sorted_mech(
    spc_str, mech_str, isolate_spc, sort_lst, spc_therm_dct=spc_therm_dct, thresh_flt_groups=50.)
rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
    rxn_sort_dct=cmts_dct)

# Write the output files (need to make general at some point)
headers = sparser.csv_headers(mech_spc_dct)
sortd_csv_str = sparser.csv_string(mech_spc_dct, headers)
sortd_mech_str = chemkin_io.writer.mechanism.write_chemkin_file(
    mech_spc_dct=mech_spc_dct,
    rxn_param_dct=param_dct_sort,
    rxn_cmts_dct=rxn_cmts_dct)

pes_groups_str = chemkin_io.writer.pesgroups.write_pes_groups(pes_groups)
# replace all '!' with '#' in mech string
sortd_mech_str = sortd_mech_str.replace('!', '#')
pathtools.write_file(sortd_csv_str, CWD, OPTS['outspc'])
pathtools.write_file(sortd_mech_str, CWD, OPTS['outmech'])
pathtools.write_file(pes_groups_str, CWD, 'pes_groups.dat')
np.savetxt(CWD+'/rxns_prompt_dh.out', rxns_filter[2:,:], delimiter='\t\t', header='\t'.join(rxns_filter[1,:]), fmt = list(rxns_filter[0,:]))
# pathtools.write_file(rxns_filter, CWD, 'rxns_prompt_dh.out')

