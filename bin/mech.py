"""
Read the mechanism file
"""

import os
import sys
import mechanalyzer


CWD = os.getcwd()
SPC_NAME = sys.argv[1]
MECH_NAME = sys.argv[2]

# Read input species file
with open(os.path.join(CWD, SPC_NAME), 'r') as file_obj:
    SPC_STR = file_obj.read()

# Read input species file
with open(os.path.join(CWD, MECH_NAME), 'r') as file_obj:
    MECH_STR = file_obj.read()

# Build the spc dct
mechanalyzer.parser.spc.write_stereo_csv(
    SPC_STR, outname=OUTNAME, path=CWD, allstereo=ALLSTEREO)

# Build the PES dct
mech_info = mechanism_file(mech_str, mech_type, spc_dct)
# Build basic pes dct unsorted
pes_dct = mechanalyzer.parser.pes.build_pes_dct(*mech_info)
# Build a sorted PES dct
# class restriction
write_mechanism_file(pes_dct, path, outname)

