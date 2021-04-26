"""
test the sorter
"""

import os
import mechanalyzer


CWD = os.getcwd()

# FUNCTIONS TO PROCESS DIFFERENT MECHANISMS IN data/ WITH DIFFERENT SORTING OPTIONS
def test__pesdct():
    """
    sort by subpes for a subset of species and get the corresponding subpes dictionaries
    """
    SPC_NAME = os.path.join(CWD, 'data', 'LLNL_species.csv')
    MECH_NAME = os.path.join(CWD, 'data', 'LLNL_mech.dat')
    ISOLATE_SPECIES = ['IC8', 'IC8-1R', 'IC8-3R', 'IC8-4R', 'IC8-5R']
    SORT_STR = ['subpes', 0]  # HEADER INDICATING THE SPECIES SUBSET
    sort_main(CWD, SPC_NAME, MECH_NAME, ISOLATE_SPECIES, SORT_STR)

def sort_main(CWD, SPC_NAME, MECH_NAME, ISOLATE_SPECIES, SORT_STR):
    spc_dct_full, rxn_param_dct, elem_tuple = mechanalyzer.parser.mech.readfiles(
        os.path.join(CWD, SPC_NAME), os.path.join(CWD, MECH_NAME))

    # BUILD  MECH INFORMATION
    mech_info = mechanalyzer.parser.mech.build_dct(spc_dct_full, rxn_param_dct)

    # SORTING: sort the mech and build the sorted rxn param dct
    srt_mch = mechanalyzer.parser.mech.sorting(
        mech_info, spc_dct_full, SORT_STR, ISOLATE_SPECIES)

    # get the pes dictionary
    pes_dct = mechanalyzer.parser.mech.sorted_pes_dct(srt_mch)
    print(pes_dct)

if __name__ == '__main__':
    test__pesdct()
