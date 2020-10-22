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
with open(os.path.join(CWD, MECH_NAME), 'r') as file_obj:
    MECH_STR = file_obj.read()

# (1) Build ini pes dct
# (2) Modify the pes dct
# (3) Write the pes_dct to a mech str
# (4) Write the file

