"""
  Add stereochemistry to a species lst
"""

import os
import sys
import mechanalyzer
import time

CWD = os.getcwd()
INNAME = sys.argv[1]
OUTNAME = sys.argv[2]
PARALLEL = False
if len(sys.argv) > 3:
    if 'p' in sys.argv[3].lower() :
        PARALLEL = True    

t0 = time.time()
# Read input species file
with open(os.path.join(CWD, INNAME), 'r') as file_obj:
    SPC_STR = file_obj.read()

# Write new string
mechanalyzer.parser.spc.write_basis_csv(
    SPC_STR, outname=OUTNAME, path=CWD, parallel=PARALLEL)
tf = time.time()

print('time to complete: {:.2f}'.format(tf-t0))
