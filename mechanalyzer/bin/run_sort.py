"""
Sort by selected criteria written in this script
"""
import sys
import os
import tempfile
from ioformat import pathtools
from mechanalyzer.builder import sorter
from mechanalyzer.parser import mech as mparser

# Set Paths to test/data directory and output directory
CWD = os.path.dirname(os.path.realpath(__file__))
TMP_OUT = tempfile.mkdtemp()

# Filenames
try:
    spc_name = os.path.join(CWD, sys.argv[1])
    mech_name = os.path.join(CWD, sys.argv[2])
    sort_inp = os.path.join(CWD, sys.argv[3])
except IndexError:
    print('*ERROR: input files missing - put species, mechanism, and sort.dat files')
    sys.exit()

sort_str = pathtools.read_file(CWD, sort_inp, remove_comments='#')
isolate_species, sort_list = mparser.read_sort_section(sort_str)
sortmech_name = os.path.join(TMP_OUT, 'sorted_mech.txt')
mech_rest_name = os.path.join(TMP_OUT, 'rest_mech.txt')
sorter._sort_main(spc_name, mech_name, sortmech_name,
                  mech_rest_name, isolate_species, sort_list)
