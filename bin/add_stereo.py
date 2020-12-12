"""
  Add stereochemistry to a species lst
"""

import os
import sys
import mechanalyzer


CWD = os.getcwd()
INNAME = sys.argv[1]
OUTNAME = sys.argv[2]

# Set default to only grab one stereoisomer
ALLSTEREO = False

if __name__ == "__main__":
    if len(sys.argv) > 3:
        if sys.argv[3] == 'allstereo':
            ALLSTEREO = True

# Read input species file
    with open(os.path.join(CWD, INNAME), 'r') as file_obj:
        spc_str = file_obj.read()

# Write new string
    mechanalyzer.parser.spc.write_stereo_csv(
        spc_str, outname=OUTNAME, path=CWD, allstereo=ALLSTEREO)
