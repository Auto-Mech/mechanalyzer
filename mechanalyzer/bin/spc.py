"""
Read the mechanism file
"""

import os
import sys


CWD = os.getcwd()
SPC_NAME = sys.argv[1]
MECH_NAME = sys.argv[2]

# Read input species file
with open(os.path.join(CWD, SPC_NAME), 'r') as file_obj:
    SPC_STR = file_obj.read()


# (1) Build ini spc dct
# (2) Modify the spc dct
# (3) Write the spc_dct to a csv str
# (4) Write the file

# OLD
# mechanalyzer.parser.spc.write_stereo_csv(
#     SPC_STR, outname=OUTNAME, path=CWD, allstereo=ALLSTEREO)
